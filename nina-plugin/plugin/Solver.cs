using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;

namespace FasterAstap {

    /// <summary>
    /// The solver executable this plugin ships, as far as the plugin is
    /// concerned: where it is, the settings file it reads, the log it writes,
    /// and the index cache it leaves behind.
    ///
    /// N.I.N.A. launches the executable itself, once per plate solve, so there
    /// is no process here to start or keep alive. What the plugin does is write
    /// down the settings that ASTAP's option set has nowhere to carry, and know
    /// where everything lives so it can be taken away again.
    /// </summary>
    public class Solver {

        public Solver(string executablePath) {
            ExecutablePath = executablePath;
            var here = Path.GetDirectoryName(executablePath) ?? ".";
            ConfigPath = Path.Combine(here, "faster-astap.ini");
            LogPath = Path.Combine(here, "faster-astap.log");
            LegacyServerPath = Path.Combine(here, "astap_index_server.exe");
        }

        public string ExecutablePath { get; }
        public string ConfigPath { get; }
        public string LogPath { get; }

        /// <summary>
        /// The resident server earlier versions shipped in this same folder. It
        /// is gone: an index cache read through a memory mapping stays in the
        /// operating system's page cache between runs, which is what the server
        /// was holding it for. An installation updated over one of those still
        /// has the executable, and possibly a copy of it running, so it is
        /// named here to be stopped, removed, and recognised in a saved ASTAP
        /// path that was pointed at it.
        /// </summary>
        public string LegacyServerPath { get; }

        public bool Present => File.Exists(ExecutablePath);

        /// <summary>
        /// Writes the settings a fixed ASTAP option set cannot carry. The solver
        /// reads this file from beside itself on every run, so a solve launched
        /// by N.I.N.A. finds the configuration set here.
        /// </summary>
        public void WriteConfig(string databaseDirectory, string tiers, double maxTier, int threads) {
            var text = new StringBuilder();
            text.AppendLine("# Written by the Faster ASTAP plugin for N.I.N.A.");
            text.AppendLine("# Edited by hand it will be overwritten the next time the plugin starts.");
            text.AppendLine();
            text.AppendLine("database = " + (databaseDirectory ?? ""));
            text.AppendLine("tiers = " + (tiers ?? ""));
            text.AppendLine("maxtier = " + maxTier.ToString(CultureInfo.InvariantCulture));
            text.AppendLine("threads = " + threads.ToString(CultureInfo.InvariantCulture));
            text.AppendLine("logfile = " + LogPath);
            File.WriteAllText(ConfigPath, text.ToString());
        }

        /// <summary>
        /// Reads the index without solving anything, building the rungs that
        /// have never been built.
        ///
        /// Reading a cached ladder is nothing — the files are memory mapped, and
        /// the pages a solve touches stay in the operating system's cache for
        /// the solves after it. Building one that has never been built is
        /// minutes, and this is how those minutes can be spent at a chosen
        /// moment rather than on the first frame of a night.
        /// </summary>
        public async Task<string> PrepareAsync(CancellationToken token = default) {
            if (!Present) return "No solver to run.";
            var info = new ProcessStartInfo {
                FileName = ExecutablePath,
                Arguments = "-prepare -config \"" + ConfigPath + "\"",
                UseShellExecute = false,
                CreateNoWindow = true,
                RedirectStandardOutput = true,
                RedirectStandardError = true,
                WorkingDirectory = Path.GetDirectoryName(ExecutablePath) ?? "."
            };
            using (var process = new Process { StartInfo = info }) {
                process.Start();
                var output = await process.StandardOutput.ReadToEndAsync().ConfigureAwait(false);
                await process.WaitForExitAsync(token).ConfigureAwait(false);
                var lines = output.Split('\n')
                                  .Select(l => l.TrimEnd('\r'))
                                  .Where(l => l.Length > 0)
                                  .ToList();
                if (lines.Count == 0)
                    return process.ExitCode == 0
                               ? "The index is ready."
                               : "The index could not be read. See " + LogPath;
                return lines[lines.Count - 1];
            }
        }

        /// <summary>
        /// The last few lines the solver wrote, which is the only account of a
        /// solve there is: N.I.N.A. launches the executable with no console
        /// anyone can see, so what it prints goes to this file and nowhere else.
        /// </summary>
        public string RecentLog(int lines) {
            try {
                if (!File.Exists(LogPath)) return "";
                var tail = new Queue<string>();
                // Shared read: a solve may be appending to this very file.
                using (var stream = new FileStream(LogPath, FileMode.Open, FileAccess.Read,
                                                   FileShare.ReadWrite))
                using (var reader = new StreamReader(stream)) {
                    string line;
                    while ((line = reader.ReadLine()) != null) {
                        if (line.Trim().Length == 0) continue;
                        tail.Enqueue(line);
                        if (tail.Count > lines) tail.Dequeue();
                    }
                }
                return string.Join(Environment.NewLine, tail);
            } catch (Exception) {
                return "";
            }
        }

        // --- what this leaves on the disk ---------------------------------------
        //
        // Three things outlive the plugin folder: the index cache, and the settings
        // and log written beside the executable. The last two go when the folder
        // does. The cache does not — it is gigabytes, it lives under LocalAppData
        // so that it survives plugin updates, and nothing else will ever tidy it up.

        /// <summary>
        /// Where the solver caches index rungs when it has not been told otherwise.
        /// Matches default_index_cache_path in the C++ side.
        /// </summary>
        public static string DefaultCacheDirectory {
            get {
                var local = Environment.GetFolderPath(Environment.SpecialFolder.LocalApplicationData);
                return Path.Combine(local, "faster-astap", "cache");
            }
        }

        public static long DirectorySize(string directory) {
            try {
                if (File.Exists(directory)) return new FileInfo(directory).Length;
                if (!Directory.Exists(directory)) return 0;
                long total = 0;
                foreach (var file in Directory.EnumerateFiles(directory, "*", SearchOption.AllDirectories)) {
                    try {
                        total += new FileInfo(file).Length;
                    } catch (Exception) {
                        // A file that vanished between listing and measuring is not
                        // worth failing a size report over.
                    }
                }
                return total;
            } catch (Exception) {
                return 0;
            }
        }

        /// <summary>
        /// Deletes the index cache. It is rebuilt on next use, so the cost of being
        /// wrong here is time rather than data.
        /// </summary>
        public static bool DeleteCache(string directory, out string error) {
            error = "";
            try {
                if (!Directory.Exists(directory)) return true;
                Directory.Delete(directory, true);
                // Take the parent too when this emptied it, so nothing is left
                // behind but no unrelated neighbour is removed.
                var parent = Path.GetDirectoryName(directory);
                if (parent != null && Directory.Exists(parent) &&
                    !Directory.EnumerateFileSystemEntries(parent).Any())
                    Directory.Delete(parent);
                return true;
            } catch (Exception e) {
                error = e.Message;
                return false;
            }
        }

        /// <summary>
        /// Stops and removes the resident server an earlier version installed
        /// into this folder, so that updating over one does not leave a process
        /// holding gigabytes for a design that no longer uses it.
        ///
        /// Best effort throughout: a server that will not stop keeps its own
        /// executable locked, and the file then stays until the next start. It
        /// is inert either way — nothing launches it any more.
        /// </summary>
        public void RemoveLegacyServer() {
            if (!File.Exists(LegacyServerPath)) return;
            try {
                foreach (var process in Process.GetProcessesByName(
                             Path.GetFileNameWithoutExtension(LegacyServerPath))) {
                    using (process) {
                        // -stop asks it over its own channel, which lets it
                        // finish the solve it may be on. It also watched the
                        // process that started it and exits when that goes, so
                        // this is belt and braces.
                        try {
                            using (var stop = Process.Start(new ProcessStartInfo {
                                FileName = LegacyServerPath,
                                Arguments = "-stop",
                                UseShellExecute = false,
                                CreateNoWindow = true
                            }))
                                stop?.WaitForExit(5000);
                        } catch (Exception) {
                            // Removed, or not runnable. Either way there is
                            // nothing left to ask.
                        }
                        process.WaitForExit(5000);
                    }
                }
            } catch (Exception) {
                // Enumerating processes can fail on a locked-down machine, and
                // that is not a reason to leave the file behind untried.
            }
            try {
                File.Delete(LegacyServerPath);
            } catch (Exception) {
                // Still running and holding itself open. Next start, then.
            }
        }

        /// <summary>
        /// The settings file and log this plugin writes beside its executable.
        /// </summary>
        public void DeleteGeneratedFiles() {
            foreach (var path in new[] { ConfigPath, LegacyServerPath, LogPath }) {
                try {
                    if (File.Exists(path)) File.Delete(path);
                } catch (Exception) {
                    // Being unable to remove a log is not worth reporting on the way
                    // out; the folder holding it is about to go anyway.
                }
            }
        }
    }
}
