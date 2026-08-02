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

            profileService.ProfileChanged += OnProfileChanged;

            // The settings are per profile and so is ASTAPLocation, so everything
            // below has to be redone when the profile changes.
            Forget(ApplyAsync(startServer: true));
        }

        public override Task Teardown() {
            profileService.ProfileChanged -= OnProfileChanged;
            // Put the ASTAP path back before going away, or a disabled plugin
            // leaves N.I.N.A. pointing at an executable that is about to be
            // deleted with it.
            RestoreAstapLocation();
            if (StopServerOnExit) {
                try {
                    server.StopAsync().Wait(TimeSpan.FromSeconds(5));
                } catch (Exception) {
                    // Going away regardless; a server left running is a held
                    // gigabyte, not a lost frame.
                }
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

        public bool StopServerOnExit {
            get => settings.GetValueBoolean(nameof(StopServerOnExit), true);
            set {
                settings.SetValueBoolean(nameof(StopServerOnExit), value);
                RaisePropertyChanged();
            }
        }

        /// <summary>
        /// Whether N.I.N.A.'s ASTAP path points at this plugin's solver. Turning it
        /// off restores whatever it was before.
        /// </summary>
        public bool UseForPlateSolving {
            get => settings.GetValueBoolean(nameof(UseForPlateSolving), false);
            set {
                settings.SetValueBoolean(nameof(UseForPlateSolving), value);
                RaisePropertyChanged();
                if (value) ApplyAstapLocation();
                else RestoreAstapLocation();
                RaisePropertyChanged(nameof(AstapLocationInUse));
                Forget(ApplyAsync(startServer: value));
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

        public ICommand RefreshCommand { get; }
        public ICommand StartServerCommand { get; }
        public ICommand StopServerCommand { get; }

        // --- doing it ----------------------------------------------------------

        private async Task ApplyAsync(bool startServer) {
            try {
                if (SolverPresent)
                    server.WriteConfig(DatabaseDirectory, Tiers, MaxTier, Threads, IdleExitMinutes);
                if (UseForPlateSolving) ApplyAstapLocation();
                if (startServer && UseForPlateSolving && SolverPresent) await StartServerAsync();
                else await RefreshStatusAsync();
            } catch (Exception e) {
                Activity = "Could not apply the settings: " + e.Message;
            }
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
