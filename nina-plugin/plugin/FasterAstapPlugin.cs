using NINA.Plugin;
using NINA.Plugin.Interfaces;
using NINA.Profile;
using NINA.Profile.Interfaces;
using NINA.WPF.Base.Interfaces.ViewModel;
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
    /// Makes N.I.N.A. plate solve against a whole-sky quad index held in memory.
    ///
    /// N.I.N.A. 3.x has no extension point for plate solvers: PlateSolverEnum is a
    /// closed set and PlateSolverFactory switches on it. What it does offer is that
    /// its ASTAP solver is a CLI solver — it launches whatever executable
    /// PlateSolveSettings.ASTAPLocation names, passes ASTAP's options, and reads the
    /// .ini written beside the image. astap_index_server takes those options and
    /// writes that .ini, so pointing that setting at it is the whole integration.
    ///
    /// This class is therefore not a solver. It is the part that keeps the index
    /// warm and the settings straight: it starts the server when N.I.N.A. starts,
    /// writes the configuration the fixed option set cannot carry, and swaps the
    /// ASTAP path over — keeping the old value so that turning the plugin off puts
    /// real ASTAP back exactly as it was.
    /// </summary>
    [Export(typeof(IPluginManifest))]
    public class FasterAstapPlugin : PluginBase, INotifyPropertyChanged {
        private const string OriginalAstapLocationKey = "OriginalAstapLocation";

        private readonly IProfileService profileService;
        private readonly IPluginOptionsAccessor settings;
        private readonly IndexServer server;
        private IndexServerStatus status = new IndexServerStatus();
        private string activity = "";

        [ImportingConstructor]
        public FasterAstapPlugin(IProfileService profileService) {
            this.profileService = profileService;
            settings = new PluginOptionsAccessor(profileService, Guid.Parse(Identifier));

            var here = Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location) ?? ".";
            server = new IndexServer(Path.Combine(here, "astap_index_server.exe"));

            RefreshCommand = new DelegateCommand(() => Forget(RefreshStatusAsync()));
            StartServerCommand = new DelegateCommand(() => Forget(StartServerAsync()));
            StopServerCommand = new DelegateCommand(() => Forget(StopServerAsync()));
            CleanUpCommand = new DelegateCommand(() => Forget(CleanUpAsync()));

            profileService.ProfileChanged += OnProfileChanged;

            // The settings are per profile and so is ASTAPLocation, so everything
            // below has to be redone when the profile changes.
            Forget(ApplyAsync(startServer: true));
        }

        /// <summary>
        /// Runs when N.I.N.A. closes, which is also the last time this code runs
        /// before an uninstall completes.
        ///
        /// The ASTAP path is put back on every exit, not only on uninstall. That is
        /// the property that makes uninstalling safe without a hook for it: the path
        /// points here only while this is loaded and serving, so a plugin folder that
        /// disappears between one run and the next cannot stand the setting up
        /// against an executable that is no longer there. It is re-adopted on the
        /// next start, once the index is resident again.
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

            if (StopServerOnExit || uninstalled) {
                try {
                    server.StopAsync().Wait(TimeSpan.FromSeconds(5));
                } catch (Exception) {
                    // Going away regardless. The server also watches this process and
                    // exits when it does, so a missed stop costs seconds, not a
                    // stranded three gigabytes.
                }
            }

            if (uninstalled) {
                server.DeleteGeneratedFiles();
                if (RemoveCacheOnUninstall) IndexServer.DeleteCache(CacheDirectory, out _);
            }
            return base.Teardown();
        }

        private void OnProfileChanged(object sender, EventArgs e) {
            RaisePropertyChanged(nameof(DatabaseDirectory));
            RaisePropertyChanged(nameof(Tiers));
            RaisePropertyChanged(nameof(MaxTier));
            RaisePropertyChanged(nameof(Threads));
            RaisePropertyChanged(nameof(IdleExitMinutes));
            RaisePropertyChanged(nameof(UseForPlateSolving));
            RaisePropertyChanged(nameof(StartWithNina));
            RaisePropertyChanged(nameof(StopServerOnExit));
            Forget(ApplyAsync(startServer: false));
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
                Forget(ApplyAsync(startServer: false));
            }
        }

        /// <summary>
        /// Depth ladder to hold, empty for the default twelve rungs. Fewer rungs
        /// is less memory: the two deepest are two thirds of the whole ladder.
        /// </summary>
        public string Tiers {
            get => settings.GetValueString(nameof(Tiers), "");
            set {
                settings.SetValueString(nameof(Tiers), value ?? "");
                RaisePropertyChanged();
                Forget(ApplyAsync(startServer: false));
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
                Forget(ApplyAsync(startServer: false));
            }
        }

        public int Threads {
            get => settings.GetValueInt32(nameof(Threads), 0);
            set {
                settings.SetValueInt32(nameof(Threads), value);
                RaisePropertyChanged();
                Forget(ApplyAsync(startServer: false));
            }
        }

        /// <summary>
        /// Stop a server that has not been asked to solve anything for this many
        /// minutes. Zero leaves it up, which is the right answer while imaging.
        /// </summary>
        public int IdleExitMinutes {
            get => settings.GetValueInt32(nameof(IdleExitMinutes), 0);
            set {
                settings.SetValueInt32(nameof(IdleExitMinutes), value);
                RaisePropertyChanged();
                Forget(ApplyAsync(startServer: false));
            }
        }

        /// <summary>
        /// Read the index into memory as soon as N.I.N.A. starts, rather than on the
        /// first solve that needs it.
        ///
        /// This is the point of the plugin: the warm-up is seconds, and paying it
        /// while N.I.N.A. is still starting up costs nothing, whereas paying it on
        /// the first solve of the night puts it in the way of something that is
        /// waiting. Separate from taking the ASTAP path over, because holding the
        /// index ready is useful even when something else is doing the solving —
        /// the command line client will use a running server whoever launched it.
        /// </summary>
        public bool StartWithNina {
            get => settings.GetValueBoolean(nameof(StartWithNina), true);
            set {
                settings.SetValueBoolean(nameof(StartWithNina), value);
                RaisePropertyChanged();
                if (value) Forget(ApplyAsync(startServer: true));
            }
        }

        public bool StopServerOnExit {
            get => settings.GetValueBoolean(nameof(StopServerOnExit), true);
            set {
                settings.SetValueBoolean(nameof(StopServerOnExit), value);
                RaisePropertyChanged();
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
        /// find a checkbox to get what they installed it for. What it does not do is
        /// take the path over the moment it is switched on: see AdoptIfReady.
        /// </summary>
        public bool UseForPlateSolving {
            get => settings.GetValueBoolean(nameof(UseForPlateSolving), true);
            set {
                settings.SetValueBoolean(nameof(UseForPlateSolving), value);
                RaisePropertyChanged();
                if (value) {
                    Forget(ApplyAsync(startServer: true));
                } else {
                    RestoreAstapLocation();
                    Activity = "";
                }
            }
        }

        // --- what is going on --------------------------------------------------

        public string SolverPath => server.ExecutablePath;
        public bool SolverPresent => File.Exists(server.ExecutablePath);
        public string AstapLocationInUse =>
            profileService?.ActiveProfile?.PlateSolveSettings?.ASTAPLocation ?? "";
        public string OriginalAstapLocation => settings.GetValueString(OriginalAstapLocationKey, "");

        public IndexServerStatus Status {
            get => status;
            private set {
                status = value;
                RaisePropertyChanged();
                RaisePropertyChanged(nameof(StatusSummary));
                RaisePropertyChanged(nameof(TimingSummary));
                RaisePropertyChanged(nameof(CacheDirectory));
                RaisePropertyChanged(nameof(CacheSummary));
            }
        }

        public string StatusSummary {
            get {
                if (!SolverPresent)
                    return "astap_index_server.exe is missing from the plugin folder.";
                if (!Status.Running) return "Not running. The index is not in memory.";
                return string.Format(CultureInfo.InvariantCulture,
                                     "Ready — {0}, {1} tiers, {2:0.00} GB in memory (pid {3}).",
                                     string.IsNullOrEmpty(Status.Database) ? "no database" : Status.Database,
                                     Status.Tiers, Status.Bytes / 1e9, Status.ProcessId);
            }
        }

        public string TimingSummary {
            get {
                if (!Status.Running || Status.Solves == 0) return "No solves yet this session.";
                return string.Format(CultureInfo.InvariantCulture,
                                     "{0} solves, median {1:0.#} ms, last {2:0.#} ms.",
                                     Status.Solves, Status.MedianMs, Status.LastMs);
            }
        }

        public string Activity {
            get => activity;
            private set {
                activity = value;
                RaisePropertyChanged();
            }
        }

        /// <summary>
        /// Where the index rungs are cached. A running server is the authority,
        /// since it may have been pointed elsewhere; otherwise the default applies.
        /// </summary>
        public string CacheDirectory =>
            !string.IsNullOrWhiteSpace(Status.CachePath) ? Status.CachePath
                                                        : IndexServer.DefaultCacheDirectory;

        public string CacheSummary {
            get {
                var bytes = IndexServer.DirectorySize(CacheDirectory);
                if (bytes == 0) return "No index cache on disk. " + CacheDirectory;
                return string.Format(CultureInfo.InvariantCulture, "{0:0.00} GB cached in {1}",
                                     bytes / 1e9, CacheDirectory);
            }
        }

        public ICommand RefreshCommand { get; }
        public ICommand StartServerCommand { get; }
        public ICommand StopServerCommand { get; }
        public ICommand CleanUpCommand { get; }

        /// <summary>
        /// Undoes everything this plugin has done to the machine, without waiting
        /// for an uninstall: the ASTAP path goes back, the server stops, and the
        /// cache and the generated files are removed. What it cannot do is delete
        /// the plugin itself, so it says what is left to do.
        /// </summary>
        private async Task CleanUpAsync() {
            RestoreAstapLocation();
            await server.StopAsync();
            await Task.Delay(500);
            var directory = CacheDirectory;
            var freed = IndexServer.DirectorySize(directory);
            server.DeleteGeneratedFiles();
            var done = IndexServer.DeleteCache(directory, out var error);
            await RefreshStatusAsync();
            RaisePropertyChanged(nameof(CacheSummary));
            Activity = done
                           ? string.Format(CultureInfo.InvariantCulture,
                                           "Removed {0:0.00} GB and restored the ASTAP path. The plugin itself is removed from Plugins ▸ Installed.",
                                           freed / 1e9)
                           : "Restored the ASTAP path, but the cache could not be deleted: " + error;
        }

        // --- doing it ----------------------------------------------------------

        private async Task ApplyAsync(bool startServer) {
            try {
                if (SolverPresent)
                    server.WriteConfig(DatabaseDirectory, Tiers, MaxTier, Threads, IdleExitMinutes);
                // Wanting the path taken over implies wanting something to hand it
                // to, so either setting is reason enough to start the server.
                if (startServer && SolverPresent && (StartWithNina || UseForPlateSolving))
                    await StartServerAsync();
                else await RefreshStatusAsync();
                AdoptIfReady();
            } catch (Exception e) {
                Activity = "Could not apply the settings: " + e.Message;
            }
        }

        /// <summary>
        /// Takes the ASTAP path over, but only once there is an index in memory to
        /// take it over with.
        ///
        /// The order matters more than it looks. Claiming the path first and warming
        /// up second would leave N.I.N.A. pointed at a solver that cannot answer for
        /// as long as the warm-up takes — seconds when the ladder is cached, minutes
        /// the first time it has to be built, and forever if the star database is
        /// wrong. A plate solve landing in that window is a lost frame, so the path
        /// is left alone until a solve would succeed, and the options page says so.
        /// </summary>
        private void AdoptIfReady() {
            if (!UseForPlateSolving || !SolverPresent) return;
            if (!Status.Running || Status.Tiers == 0) {
                Activity = "Plate solving is left with ASTAP until the index is in memory.";
                return;
            }
            ApplyAstapLocation();
            Activity = "";
        }

        private async Task StartServerAsync() {
            if (!SolverPresent) {
                Activity = "No solver to start.";
                return;
            }
            Activity = "Reading the index into memory…";
            // Reading the whole ladder off disk is seconds; building it the first
            // time, on a machine that has never done it, is minutes.
            var ok = await server.StartAsync(TimeSpan.FromMinutes(30));
            Activity = ok ? "" : "The server did not come up. See " + server.LogPath;
            await RefreshStatusAsync();
            AdoptIfReady();
        }

        private async Task StopServerAsync() {
            Activity = "Stopping…";
            await server.StopAsync();
            await Task.Delay(500);
            Activity = "";
            await RefreshStatusAsync();
        }

        private async Task RefreshStatusAsync() {
            Status = await server.GetStatusAsync();
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
            if (SamePath(current, server.ExecutablePath)) return;
            if (!string.IsNullOrWhiteSpace(current))
                settings.SetValueString(OriginalAstapLocationKey, current);
            plateSolve.ASTAPLocation = server.ExecutablePath;
            RaisePropertyChanged(nameof(AstapLocationInUse));
        }

        private void RestoreAstapLocation() {
            var plateSolve = profileService?.ActiveProfile?.PlateSolveSettings;
            if (plateSolve == null) return;
            if (!SamePath(plateSolve.ASTAPLocation ?? "", server.ExecutablePath)) return;
            var original = settings.GetValueString(OriginalAstapLocationKey, "");
            plateSolve.ASTAPLocation = original;
            RaisePropertyChanged(nameof(AstapLocationInUse));
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
                if (SamePath(candidate, server.ExecutablePath)) continue;
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
