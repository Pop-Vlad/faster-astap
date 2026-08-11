using NINA.Plugin;
using NINA.Plugin.Interfaces;
using NINA.Profile;
using NINA.Profile.Interfaces;
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.ComponentModel.Composition;
using System.Globalization;
using System.IO;
using System.Reflection;
using System.Runtime.CompilerServices;
using System.Threading;
using System.Threading.Tasks;
using System.Windows.Input;

namespace FasterAstap {

    /// <summary>
    /// Makes N.I.N.A. plate solve against a whole-sky quad index instead of
    /// searching the sky.
    ///
    /// N.I.N.A. 3.x has no extension point for plate solvers: PlateSolverEnum is a
    /// closed set and PlateSolverFactory switches on it. What it does offer is that
    /// its ASTAP solver is a CLI solver — it launches whatever executable
    /// PlateSolveSettings.ASTAPLocation names, passes ASTAP's options, and reads the
    /// .ini written beside the image. astap_nina_solve takes those options and
    /// writes that .ini, so pointing that setting at it is the whole integration.
    ///
    /// This class is therefore not a solver, and it is not a supervisor either:
    /// nothing of this plugin runs while a solve does. It writes the configuration
    /// the fixed option set cannot carry, and swaps the ASTAP path over — keeping
    /// the old value so that turning the plugin off puts real ASTAP back exactly as
    /// it was.
    /// </summary>
    [Export(typeof(IPluginManifest))]
    public class FasterAstapPlugin : PluginBase, INotifyPropertyChanged {
        private const string OriginalAstapLocationKey = "OriginalAstapLocation";

        private readonly IProfileService profileService;
        private readonly IPluginOptionsAccessor settings;
        private readonly Solver solver;
        private string activity = "";
        private bool busy;

        [ImportingConstructor]
        public FasterAstapPlugin(IProfileService profileService) {
            this.profileService = profileService;
            settings = new PluginOptionsAccessor(profileService, Guid.Parse(Identifier));

            var here = Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location) ?? ".";
            solver = new Solver(Path.Combine(here, "astap_nina_solve.exe"));

            RefreshCommand = new DelegateCommand(Refresh);
            // Both of these run the solver or delete gigabytes, so neither is
            // allowed to start while the other is still going.
            PrepareCommand = new DelegateCommand(() => Forget(PrepareAsync()), () => !busy);
            CleanUpCommand = new DelegateCommand(() => Forget(CleanUpAsync()), () => !busy);

            profileService.ProfileChanged += OnProfileChanged;

            // Done now, and again on every profile change: the settings are per
            // profile, and so is ASTAPLocation.
            Apply();
        }

        /// <summary>
        /// Runs when N.I.N.A. closes, which is also the last time this code runs
        /// before an uninstall completes.
        ///
        /// The ASTAP path is put back on every exit, not only on uninstall. That is
        /// the property that makes uninstalling safe without a hook for it: the path
        /// points here only while this is loaded, so a plugin folder that disappears
        /// between one run and the next cannot leave the setting standing against an
        /// executable that is no longer there. It is re-adopted on the next start.
        /// </summary>
        public override Task Teardown() {
            profileService.ProfileChanged -= OnProfileChanged;
            RestoreAstapLocation();

            // N.I.N.A. uninstalls a plugin by moving its folder away and deleting it
            // on the next start, so an uninstall during this session leaves the
            // executable beside this assembly already gone by the time this runs.
            // There is no uninstall callback, and that absence is the closest thing
            // to one.
            //
            // It is a guess, and it is only allowed to decide something a wrong
            // guess cannot break. Deleting the plugin folder while N.I.N.A. is shut
            // down never runs this at all, and then the cache is simply left on the
            // disk — which costs space and nothing else. Restoring the ASTAP path is
            // the part that would cost a night, so that happens above, on every
            // exit, without consulting this.
            var uninstalled = !SolverPresent;
            if (uninstalled) {
                solver.DeleteGeneratedFiles();
                if (RemoveCacheOnUninstall) Solver.DeleteCache(CacheDirectory, out _);
            }
            return base.Teardown();
        }

        private void OnProfileChanged(object sender, EventArgs e) {
            RaisePropertyChanged(nameof(DatabaseDirectory));
            RaisePropertyChanged(nameof(Tiers));
            RaisePropertyChanged(nameof(MaxTier));
            RaisePropertyChanged(nameof(Threads));
            RaisePropertyChanged(nameof(UseForPlateSolving));
            Apply();
        }

        // --- settings ----------------------------------------------------------

        /// <summary>
        /// Directory holding the .1476 / .290 star database files. Defaults to
        /// wherever the ASTAP that N.I.N.A. is already configured for lives, since
        /// an ASTAP installation comes with its database beside it.
        /// </summary>
        public string DatabaseDirectory {
            get {
                var stored = settings.GetValueString(nameof(DatabaseDirectory), "");
                return string.IsNullOrWhiteSpace(stored) ? GuessDatabaseDirectory() : stored;
            }
            set {
                settings.SetValueString(nameof(DatabaseDirectory), value ?? "");
                RaisePropertyChanged();
                Apply();
            }
        }

        /// <summary>
        /// Depth ladder to use, empty for the default twelve rungs. Fewer rungs is
        /// less disk and less of the page cache spent on this: the two deepest are
        /// two thirds of the whole ladder.
        /// </summary>
        public string Tiers {
            get => settings.GetValueString(nameof(Tiers), "");
            set {
                settings.SetValueString(nameof(Tiers), value ?? "");
                RaisePropertyChanged();
                Apply();
            }
        }

        /// <summary>
        /// Raises the ceiling of the default ladder, for fields below about a
        /// quarter of a degree. The deep rungs are large — 1800 adds 2.5 GB and
        /// 3600 another 4.8 — so this is off unless asked for.
        /// </summary>
        public double MaxTier {
            get => settings.GetValueDouble(nameof(MaxTier), 0);
            set {
                settings.SetValueDouble(nameof(MaxTier), value);
                RaisePropertyChanged();
                Apply();
            }
        }

        public int Threads {
            get => settings.GetValueInt32(nameof(Threads), 0);
            set {
                settings.SetValueInt32(nameof(Threads), value);
                RaisePropertyChanged();
                Apply();
            }
        }

        /// <summary>
        /// Delete the index cache when the plugin is uninstalled, so that removing
        /// it leaves nothing behind.
        ///
        /// Worth knowing before leaving this on: the cache is keyed to the star
        /// database and the tolerance, not to this plugin, so a command line
        /// astap_index_solve on the same machine is using the very same files. It
        /// costs correctness nothing — the next run rebuilds what it needs — but it
        /// can cost that run several minutes.
        /// </summary>
        public bool RemoveCacheOnUninstall {
            get => settings.GetValueBoolean(nameof(RemoveCacheOnUninstall), true);
            set {
                settings.SetValueBoolean(nameof(RemoveCacheOnUninstall), value);
                RaisePropertyChanged();
            }
        }

        /// <summary>
        /// Whether N.I.N.A.'s ASTAP path points at this plugin's solver. Turning it
        /// off restores whatever it was before.
        ///
        /// On by default, so that installing the plugin is enough and nobody has to
        /// find a checkbox to get what they installed it for.
        /// </summary>
        public bool UseForPlateSolving {
            get => settings.GetValueBoolean(nameof(UseForPlateSolving), true);
            set {
                settings.SetValueBoolean(nameof(UseForPlateSolving), value);
                RaisePropertyChanged();
                if (value) {
                    Apply();
                } else {
                    RestoreAstapLocation();
                    Activity = "";
                }
            }
        }

        // --- what is going on --------------------------------------------------

        public string SolverPath => solver.ExecutablePath;
        public bool SolverPresent => solver.Present;
        public string LogPath => solver.LogPath;
        public string AstapLocationInUse =>
            profileService?.ActiveProfile?.PlateSolveSettings?.ASTAPLocation ?? "";
        public string OriginalAstapLocation => settings.GetValueString(OriginalAstapLocationKey, "");

        public string StatusSummary {
            get {
                if (!SolverPresent)
                    return "astap_nina_solve.exe is missing from the plugin folder.";
                if (!UseForPlateSolving)
                    return "Installed, and not in use: plate solving is left with ASTAP.";
                if (!SamePath(AstapLocationInUse, solver.ExecutablePath))
                    return "Installed, but N.I.N.A.'s ASTAP path points somewhere else.";
                return "In use. N.I.N.A. plate solves against the quad index.";
            }
        }

        /// <summary>
        /// The last lines the solver wrote. N.I.N.A. runs it with no console
        /// anybody can see, so its log is where a solve — or a failure to find the
        /// star database — is reported.
        /// </summary>
        public string RecentLog {
            get {
                var tail = solver.RecentLog(6);
                return string.IsNullOrEmpty(tail) ? "Nothing solved yet." : tail;
            }
        }

        public string Activity {
            get => activity;
            private set {
                activity = value;
                RaisePropertyChanged();
            }
        }

        /// <summary>Where the index rungs are cached.</summary>
        public string CacheDirectory => Solver.DefaultCacheDirectory;

        public string CacheSummary {
            get {
                var bytes = Solver.DirectorySize(CacheDirectory);
                if (bytes == 0)
                    return "No index cache on disk yet — the first solve builds one, which takes " +
                           "minutes. " + CacheDirectory;
                return string.Format(CultureInfo.InvariantCulture, "{0:0.00} GB cached in {1}",
                                     bytes / 1e9, CacheDirectory);
            }
        }

        public ICommand RefreshCommand { get; }
        public ICommand PrepareCommand { get; }
        public ICommand CleanUpCommand { get; }

        /// <summary>
        /// Builds or reads the index without solving anything, so that a ladder
        /// that has never been built is paid for now rather than on the first frame
        /// of a night.
        /// </summary>
        private async Task PrepareAsync() {
            if (busy) return;
            if (!SolverPresent) {
                Activity = "No solver to run.";
                return;
            }
            busy = true;
            Activity = "Reading the index; building any rung that has never been built…";
            try {
                Activity = await solver.PrepareAsync();
            } finally {
                busy = false;
                Refresh();
            }
        }

        /// <summary>
        /// Undoes everything this plugin has done to the machine, without waiting
        /// for an uninstall: the ASTAP path goes back, and the cache and the
        /// generated files are removed. What it cannot do is delete the plugin
        /// itself, so it says what is left to do.
        /// </summary>
        private async Task CleanUpAsync() {
            if (busy) return;
            busy = true;
            try {
                RestoreAstapLocation();
                var directory = CacheDirectory;
                var freed = Solver.DirectorySize(directory);
                solver.DeleteGeneratedFiles();
                // Several gigabytes of small files: off the interface thread, so
                // the options page is not frozen while they go.
                var error = await Task.Run(() => Solver.DeleteCache(directory, out var e) ? "" : e);
                Activity = error.Length == 0
                               ? string.Format(CultureInfo.InvariantCulture,
                                               "Removed {0:0.00} GB and restored the ASTAP path. The plugin itself is removed from Plugins ▸ Installed.",
                                               freed / 1e9)
                               : "Restored the ASTAP path, but the cache could not be deleted: " + error;
            } finally {
                busy = false;
                Refresh();
            }
        }

        // --- doing it ----------------------------------------------------------

        /// <summary>
        /// Writes the settings file and takes the ASTAP path over.
        ///
        /// Nothing has to be warmed up first. The solver reads the index cache
        /// through a memory mapping, so a solve that follows another one finds the
        /// pages it needs already in the operating system's cache and never touches
        /// the disk for them — which is what a resident server used to be for. The
        /// one cost that is not milliseconds is building a ladder that has never
        /// been built, and the options page has a button for paying it deliberately.
        /// </summary>
        private void Apply() {
            try {
                if (!SolverPresent) {
                    Activity = "";
                    Refresh();
                    return;
                }
                // An installation updated over the version that shipped a
                // resident server may still have it, and still be running it.
                // Off this thread: asking it to stop waits on a process that has
                // a solve to finish, and this runs while N.I.N.A. is starting.
                Forget(Task.Run(() => solver.RemoveLegacyServer()));
                solver.WriteConfig(DatabaseDirectory, Tiers, MaxTier, Threads);
                if (UseForPlateSolving) ApplyAstapLocation();
                Activity = "";
                Refresh();
            } catch (Exception e) {
                Activity = "Could not apply the settings: " + e.Message;
            }
        }

        private void Refresh() {
            RaisePropertyChanged(nameof(StatusSummary));
            RaisePropertyChanged(nameof(RecentLog));
            RaisePropertyChanged(nameof(CacheSummary));
            RaisePropertyChanged(nameof(AstapLocationInUse));
        }

        /// <summary>
        /// Points N.I.N.A. at this plugin's solver, remembering what was there. The
        /// old value is only recorded when it is not already ours, so that applying
        /// this twice cannot lose the real ASTAP path.
        /// </summary>
        private void ApplyAstapLocation() {
            var plateSolve = profileService?.ActiveProfile?.PlateSolveSettings;
            if (plateSolve == null || !SolverPresent) return;
            var current = plateSolve.ASTAPLocation ?? "";
            if (SamePath(current, solver.ExecutablePath)) return;
            // A path left pointing at the resident server of an earlier version
            // is one of ours, not a real ASTAP: remembering it as the original
            // would replace the path this is supposed to be able to give back.
            if (!string.IsNullOrWhiteSpace(current) && !SamePath(current, solver.LegacyServerPath))
                settings.SetValueString(OriginalAstapLocationKey, current);
            plateSolve.ASTAPLocation = solver.ExecutablePath;
            RaisePropertyChanged(nameof(AstapLocationInUse));
        }

        private void RestoreAstapLocation() {
            var plateSolve = profileService?.ActiveProfile?.PlateSolveSettings;
            if (plateSolve == null) return;
            var current = plateSolve.ASTAPLocation ?? "";
            // The earlier version's server counts as ours here too, so a profile
            // that was pointed at it gets the real ASTAP back rather than a path
            // to an executable this no longer installs.
            if (!SamePath(current, solver.ExecutablePath) &&
                !SamePath(current, solver.LegacyServerPath)) return;
            var original = settings.GetValueString(OriginalAstapLocationKey, "");
            plateSolve.ASTAPLocation = original;
            RaisePropertyChanged(nameof(AstapLocationInUse));
            RaisePropertyChanged(nameof(StatusSummary));
        }

        /// <summary>
        /// The database that came with the ASTAP N.I.N.A. was already pointed at.
        /// Read from the remembered original first, because by the time this is
        /// asked the live setting may already be ours.
        /// </summary>
        private string GuessDatabaseDirectory() {
            var candidates = new List<string> {
                OriginalAstapLocation,
                profileService?.ActiveProfile?.PlateSolveSettings?.ASTAPLocation ?? "",
                @"C:\Program Files\astap\astap.exe"
            };
            foreach (var candidate in candidates) {
                if (string.IsNullOrWhiteSpace(candidate)) continue;
                if (SamePath(candidate, solver.ExecutablePath) ||
                    SamePath(candidate, solver.LegacyServerPath)) continue;
                var directory = Path.GetDirectoryName(candidate);
                if (!string.IsNullOrEmpty(directory) && Directory.Exists(directory) &&
                    HasStarDatabase(directory))
                    return directory;
            }
            return "";
        }

        private static bool HasStarDatabase(string directory) {
            try {
                foreach (var pattern in new[] { "*.1476", "*.290", "*.001" })
                    if (Directory.EnumerateFiles(directory, pattern).GetEnumerator().MoveNext())
                        return true;
            } catch (Exception) {
                return false;
            }
            return false;
        }

        private static bool SamePath(string a, string b) {
            if (string.IsNullOrWhiteSpace(a) || string.IsNullOrWhiteSpace(b)) return false;
            try {
                return string.Equals(Path.GetFullPath(a), Path.GetFullPath(b),
                                     StringComparison.OrdinalIgnoreCase);
            } catch (Exception) {
                return string.Equals(a, b, StringComparison.OrdinalIgnoreCase);
            }
        }

        /// <summary>
        /// Runs a task without waiting for it, and without letting a failure inside
        /// it reach the finalizer thread as an unobserved exception.
        /// </summary>
        private void Forget(Task task) {
            _ = task.ContinueWith(t => { Activity = "Error: " + t.Exception?.GetBaseException().Message; },
                                  CancellationToken.None, TaskContinuationOptions.OnlyOnFaulted,
                                  TaskScheduler.Default);
        }

        public event PropertyChangedEventHandler PropertyChanged;

        protected void RaisePropertyChanged([CallerMemberName] string propertyName = null) {
            PropertyChanged?.Invoke(this, new PropertyChangedEventArgs(propertyName));
        }
    }
}
