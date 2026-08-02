using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.IO.Pipes;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;

namespace FasterAstap {

    /// <summary>
    /// What a running server reports about itself.
    /// </summary>
    public class IndexServerStatus {
        public bool Running { get; set; }
        public string Database { get; set; } = "";
        public string DatabasePath { get; set; } = "";
        public string CachePath { get; set; } = "";
        public int Tiers { get; set; }
        public long Quads { get; set; }
        public long Bytes { get; set; }
        public string Densities { get; set; } = "";
        public long Solves { get; set; }
        public double MedianMs { get; set; }
        public double LastMs { get; set; }
        public double UptimeSeconds { get; set; }
        public int ProcessId { get; set; }
        public string Error { get; set; } = "";
    }

    /// <summary>
    /// Talks to astap_index_server: writes its settings file, starts and stops it,
    /// and asks it what it is holding.
    ///
    /// The server is a separate process on purpose. It holds gigabytes and it runs
    /// native code, and neither of those belongs in N.I.N.A.'s own process in the
    /// middle of a night's imaging. It also means the index survives this plugin
    /// being reloaded.
    /// </summary>
    public class IndexServer {
        private const string DefaultPipe = "faster-astap-index";

        public IndexServer(string executablePath) {
            ExecutablePath = executablePath;
            ConfigPath = Path.Combine(Path.GetDirectoryName(executablePath) ?? ".", "faster-astap.ini");
            LogPath = Path.Combine(Path.GetDirectoryName(executablePath) ?? ".", "faster-astap.log");
        }

        public string ExecutablePath { get; }
        public string ConfigPath { get; }
        public string LogPath { get; }

        /// <summary>
        /// Writes the settings a fixed ASTAP option set cannot carry. Both the
        /// server and the client program read this, so a solve launched by
        /// N.I.N.A. finds the same configuration the plugin set up.
        /// </summary>
        public void WriteConfig(string databaseDirectory, string tiers, double maxTier, int threads,
                               int idleExitMinutes) {
            var text = new StringBuilder();
            text.AppendLine("# Written by the Faster ASTAP plugin for N.I.N.A.");
            text.AppendLine("# Edited by hand it will be overwritten the next time the plugin starts.");
            text.AppendLine();
            text.AppendLine("database = " + (databaseDirectory ?? ""));
            text.AppendLine("tiers = " + (tiers ?? ""));
            text.AppendLine("maxtier = " + maxTier.ToString(CultureInfo.InvariantCulture));
            text.AppendLine("threads = " + threads.ToString(CultureInfo.InvariantCulture));
            text.AppendLine("idle_exit_minutes = " + idleExitMinutes.ToString(CultureInfo.InvariantCulture));
            text.AppendLine("logfile = " + LogPath);
            // A solve that arrives with no server running starts one rather than
            // failing, which is what keeps a frame from being lost if the server
            // died or N.I.N.A. was started before the plugin got going.
            text.AppendLine("autostart = 1");
            File.WriteAllText(ConfigPath, text.ToString());
        }

        public async Task<IndexServerStatus> GetStatusAsync(CancellationToken token = default) {
            var status = new IndexServerStatus();
            var reply = await RequestAsync("op=status\n", 500, token).ConfigureAwait(false);
            if (reply == null) {
                status.Running = false;
                status.Error = "no server listening";
                return status;
            }
            var fields = Decode(reply);
            status.Running = true;
            status.Database = Get(fields, "database");
            status.DatabasePath = Get(fields, "database_path");
            status.CachePath = Get(fields, "cache");
            status.Tiers = (int)Number(fields, "tiers");
            status.Quads = (long)Number(fields, "quads");
            status.Bytes = (long)Number(fields, "bytes");
            status.Densities = Get(fields, "densities");
            status.Solves = (long)Number(fields, "solves");
            status.MedianMs = Number(fields, "median_ms");
            status.LastMs = Number(fields, "last_ms");
            status.UptimeSeconds = Number(fields, "uptime_seconds");
            status.ProcessId = (int)Number(fields, "pid");
            return status;
        }

        public async Task<bool> IsRunningAsync(CancellationToken token = default) {
            return await RequestAsync("op=ping\n", 300, token).ConfigureAwait(false) != null;
        }

        /// <summary>
        /// Starts a server and waits until it answers. Reading the index takes a
        /// few seconds, and it does not listen until it is ready, so answering at
        /// all is the readiness signal.
        /// </summary>
        public async Task<bool> StartAsync(TimeSpan timeout, CancellationToken token = default) {
            if (await IsRunningAsync(token).ConfigureAwait(false)) return true;
            if (!File.Exists(ExecutablePath)) return false;

            var info = new ProcessStartInfo {
                FileName = ExecutablePath,
                Arguments = "-serve -quiet -config \"" + ConfigPath + "\"",
                UseShellExecute = false,
                CreateNoWindow = true,
                WorkingDirectory = Path.GetDirectoryName(ExecutablePath) ?? "."
            };
            try {
                Process.Start(info);
            } catch (Exception) {
                return false;
            }

            var deadline = DateTime.UtcNow + timeout;
            while (DateTime.UtcNow < deadline) {
                if (token.IsCancellationRequested) return false;
                if (await IsRunningAsync(token).ConfigureAwait(false)) return true;
                await Task.Delay(250, token).ConfigureAwait(false);
            }
            return false;
        }

        public async Task StopAsync(CancellationToken token = default) {
            await RequestAsync("op=shutdown\n", 500, token).ConfigureAwait(false);
        }

        // --- the wire ----------------------------------------------------------
        //
        // One connection, one request, one reply, each a 32 bit little endian
        // length followed by that many bytes of UTF-8. Returns null when nothing
        // was listening, which the caller has to tell apart from a refusal.

        private static async Task<string> RequestAsync(string message, int connectTimeoutMs,
                                                       CancellationToken token) {
            try {
                using (var pipe = new NamedPipeClientStream(".", DefaultPipe, PipeDirection.InOut,
                                                            PipeOptions.Asynchronous)) {
                    try {
                        await pipe.ConnectAsync(connectTimeoutMs, token).ConfigureAwait(false);
                    } catch (TimeoutException) {
                        return null;
                    } catch (IOException) {
                        return null;
                    }

                    var payload = Encoding.UTF8.GetBytes(message);
                    var header = BitConverter.GetBytes(payload.Length);
                    if (!BitConverter.IsLittleEndian) Array.Reverse(header);
                    await pipe.WriteAsync(header, 0, 4, token).ConfigureAwait(false);
                    await pipe.WriteAsync(payload, 0, payload.Length, token).ConfigureAwait(false);
                    await pipe.FlushAsync(token).ConfigureAwait(false);

                    var replyHeader = await ReadExactAsync(pipe, 4, token).ConfigureAwait(false);
                    if (replyHeader == null) return null;
                    if (!BitConverter.IsLittleEndian) Array.Reverse(replyHeader);
                    var length = BitConverter.ToInt32(replyHeader, 0);
                    if (length < 0 || length > 64 * 1024 * 1024) return null;
                    var body = await ReadExactAsync(pipe, length, token).ConfigureAwait(false);
                    return body == null ? null : Encoding.UTF8.GetString(body);
                }
            } catch (OperationCanceledException) {
                return null;
            } catch (Exception) {
                return null;
            }
        }

        private static async Task<byte[]> ReadExactAsync(Stream stream, int count,
                                                         CancellationToken token) {
            var buffer = new byte[count];
            var read = 0;
            while (read < count) {
                var n = await stream.ReadAsync(buffer, read, count - read, token).ConfigureAwait(false);
                if (n <= 0) return null;
                read += n;
            }
            return buffer;
        }

        private static List<KeyValuePair<string, string>> Decode(string text) {
            var result = new List<KeyValuePair<string, string>>();
            foreach (var line in text.Split('\n')) {
                var trimmed = line.TrimEnd('\r');
                if (trimmed.Length == 0) continue;
                var split = trimmed.IndexOf('=');
                if (split < 0) continue;
                result.Add(new KeyValuePair<string, string>(trimmed.Substring(0, split),
                                                            trimmed.Substring(split + 1)));
            }
            return result;
        }

        private static string Get(List<KeyValuePair<string, string>> fields, string key) {
            foreach (var kv in fields)
                if (kv.Key == key) return kv.Value;
            return "";
        }

        private static double Number(List<KeyValuePair<string, string>> fields, string key) {
            var text = Get(fields, key);
            return double.TryParse(text, NumberStyles.Float, CultureInfo.InvariantCulture, out var v)
                       ? v
                       : 0;
        }

        /// <summary>
        /// Server processes started from this executable, so that a plugin
        /// reload does not leave one behind and a user can be told one is up.
        /// </summary>
        public IEnumerable<Process> RunningProcesses() {
            var name = Path.GetFileNameWithoutExtension(ExecutablePath);
            Process[] all;
            try {
                all = Process.GetProcessesByName(name);
            } catch (Exception) {
                return Enumerable.Empty<Process>();
            }
            return all;
        }
    }
}
