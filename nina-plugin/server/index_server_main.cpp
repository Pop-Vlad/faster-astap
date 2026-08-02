// The index solver with the index left in memory between images.
//
// Why this exists: `astap_index_solve` spends 1.65 s reading a 2.7 GB ladder and
// 5 ms solving. For a batch that is the right trade, because the read is paid
// once for all the images. An imaging application solves one frame at a time,
// minutes apart, in a separate process each time — so it pays the read every
// single solve and never sees the 5 ms. This binary splits the two apart.
//
// One executable, two roles:
//
//   astap_index_server -serve      holds the ladder and answers requests
//   astap_index_server -f img.fits solves through it, taking ASTAP's options
//
// The second role is what an imaging application launches, and it takes the
// option set ASTAP takes (-f -z -fov -r -ra -spd -s) and writes the same .ini,
// so anything that can call ASTAP can call this without knowing it has. N.I.N.A.
// in particular has no way to register a new plate solver, but it will launch
// whatever executable its ASTAP setting points at, which is the whole
// integration.
//
// Settings a fixed option set cannot carry — where the star database lives,
// which depth tiers to hold — come from faster-astap.ini next to the executable.
// That file is also how the client knows how to start a server itself when none
// is running, so a solve still succeeds when nobody started one.

#include <algorithm>
#include <atomic>
#include <cctype>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <initializer_list>
#include <iostream>
#include <map>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include "astap/astro_math.h"
#include "astap/image/image_io.h"
#include "ipc.h"
#include "astap/parallel.h"
#include "astap/solve_service.h"
#include "astap/solver.h"

#ifdef _WIN32
#include <windows.h>
#else
#include <cerrno>
#include <csignal>
#include <spawn.h>
#include <sys/wait.h>
#include <unistd.h>
extern char **environ;
#endif

using Clock = std::chrono::steady_clock;

namespace {
  double secs(Clock::time_point a, Clock::time_point b) {
    return std::chrono::duration<double>(b - a).count();
  }

  // --- the wire format ------------------------------------------------------
  //
  // One key=value per line, first '=' splits. Repeated keys are kept in order,
  // which is how a reply carries the several lines a solve would have printed.
  // Values hold file names and numbers, never a newline; any that arrives is
  // dropped rather than allowed to forge a second key.

  std::string sanitise(const std::string &v) {
    std::string out;
    out.reserve(v.size());
    for (char c : v)
      if (c != '\n' && c != '\r') out += c;
    return out;
  }

  class Message {
  public:
    void set(const std::string &k, const std::string &v) { items_.emplace_back(k, sanitise(v)); }
    void set(const std::string &k, long long v) { set(k, std::to_string(v)); }
    void set(const std::string &k, double v, int decimals) {
      set(k, astap::float_to_str(v, decimals));
    }

    std::string encode() const {
      std::string s;
      for (const auto &kv : items_) s += kv.first + "=" + kv.second + "\n";
      return s;
    }

    static Message decode(const std::string &text) {
      Message m;
      size_t p = 0;
      while (p < text.size()) {
        size_t e = text.find('\n', p);
        if (e == std::string::npos) e = text.size();
        const std::string line = text.substr(p, e - p);
        p = e + 1;
        const size_t eq = line.find('=');
        if (eq == std::string::npos) continue;
        m.items_.emplace_back(line.substr(0, eq), sanitise(line.substr(eq + 1)));
      }
      return m;
    }

    std::string get(const std::string &k, const std::string &fallback = "") const {
      for (const auto &kv : items_)
        if (kv.first == k) return kv.second;
      return fallback;
    }

    double number(const std::string &k, double fallback) const {
      const std::string v = get(k);
      return v.empty() ? fallback : std::atof(v.c_str());
    }

    int integer(const std::string &k, int fallback) const {
      const std::string v = get(k);
      return v.empty() ? fallback : std::atoi(v.c_str());
    }

    bool flag(const std::string &k, bool fallback) const {
      const std::string v = get(k);
      return v.empty() ? fallback : (v == "1" || v == "true" || v == "yes");
    }

    std::vector<std::string> all(const std::string &k) const {
      std::vector<std::string> out;
      for (const auto &kv : items_)
        if (kv.first == k) out.push_back(kv.second);
      return out;
    }

  private:
    std::vector<std::pair<std::string, std::string> > items_;
  };

  // --- where the executable is ----------------------------------------------

  std::string executable_path() {
#ifdef _WIN32
    char buf[MAX_PATH * 4];
    const DWORD n = GetModuleFileNameA(nullptr, buf, sizeof(buf));
    return n ? std::string(buf, n) : std::string();
#else
    char buf[4096];
    const ssize_t n = readlink("/proc/self/exe", buf, sizeof(buf) - 1);
    return n > 0 ? std::string(buf, static_cast<size_t>(n)) : std::string();
#endif
  }

  std::string directory_of(const std::string &path) {
    const size_t slash = path.find_last_of("/\\");
    return slash == std::string::npos ? std::string(".") : path.substr(0, slash);
  }

  // --- faster-astap.ini ------------------------------------------------------
  //
  // Read from next to the executable. It exists because the options an imaging
  // application passes are ASTAP's, and ASTAP has nowhere to say "hold these
  // tiers" or "the database is over there". A plugin writes this file; a person
  // can edit it.

  struct Config {
    std::string database_path;
    std::string database = "auto";
    double tolerance = 0.007;
    std::string tiers;    // empty selects the default ladder
    double max_tier = 0;  // raises the ceiling of the default ladder
    std::string cache;
    std::string endpoint;
    std::string logfile;
    int threads = 0;
    bool autostart = true;
    int idle_exit_minutes = 0;  // 0 leaves a server up until it is told to stop
  };

  Config read_config(const std::string &path) {
    Config c;
    std::ifstream f(path);
    if (!f) return c;
    std::string line;
    while (std::getline(f, line)) {
      while (!line.empty() && (line.back() == '\r' || line.back() == ' ')) line.pop_back();
      const size_t start = line.find_first_not_of(" \t");
      if (start == std::string::npos || line[start] == '#' || line[start] == ';') continue;
      const size_t eq = line.find('=', start);
      if (eq == std::string::npos) continue;
      std::string k = line.substr(start, eq - start);
      while (!k.empty() && k.back() == ' ') k.pop_back();
      std::string v = line.substr(eq + 1);
      const size_t vs = v.find_first_not_of(" \t");
      v = vs == std::string::npos ? "" : v.substr(vs);
      if (k == "database") c.database_path = v;
      else if (k == "database_abbreviation") c.database = v.empty() ? "auto" : v;
      else if (k == "tolerance") c.tolerance = std::atof(v.c_str());
      else if (k == "tiers") c.tiers = v;
      else if (k == "maxtier") c.max_tier = std::atof(v.c_str());
      else if (k == "cache") c.cache = v;
      else if (k == "endpoint") c.endpoint = v;
      else if (k == "logfile") c.logfile = v;
      else if (k == "threads") c.threads = std::atoi(v.c_str());
      else if (k == "autostart") c.autostart = v == "1" || v == "true" || v == "yes";
      else if (k == "idle_exit_minutes") c.idle_exit_minutes = std::atoi(v.c_str());
    }
    return c;
  }

  // --- starting a server from a client --------------------------------------
  //
  // Detached and windowless: the client that starts it is about to exit, and the
  // server has to outlive it without ever putting a console on the screen of
  // somebody in the middle of an imaging run.

  bool spawn_server(const std::string &exe, const std::string &config_path, std::string *error) {
#ifdef _WIN32
    std::string cmd = "\"" + exe + "\" -serve";
    if (!config_path.empty()) cmd += " -config \"" + config_path + "\"";
    STARTUPINFOA si = {};
    si.cb = sizeof(si);
    PROCESS_INFORMATION pi = {};
    std::vector<char> mutable_cmd(cmd.begin(), cmd.end());
    mutable_cmd.push_back('\0');
    const BOOL ok = CreateProcessA(nullptr, mutable_cmd.data(), nullptr, nullptr, FALSE,
                                   DETACHED_PROCESS | CREATE_NO_WINDOW, nullptr,
                                   directory_of(exe).c_str(), &si, &pi);
    if (!ok) {
      if (error) *error = "could not start a server (error " + std::to_string(GetLastError()) + ")";
      return false;
    }
    CloseHandle(pi.hThread);
    CloseHandle(pi.hProcess);
    return true;
#else
    std::vector<std::string> args = {exe, "-serve"};
    if (!config_path.empty()) {
      args.push_back("-config");
      args.push_back(config_path);
    }
    std::vector<char *> argv;
    for (std::string &a : args) argv.push_back(const_cast<char *>(a.c_str()));
    argv.push_back(nullptr);
    pid_t pid = 0;
    if (posix_spawn(&pid, exe.c_str(), nullptr, nullptr, argv.data(), environ) != 0) {
      if (error) *error = "could not start a server";
      return false;
    }
    return true;
#endif
  }

  // --- the resident server ---------------------------------------------------

  class ServerApp {
  public:
    ServerApp(astap::SolveService &service, const Config &config)
        : service_(service), config_(config), started_(Clock::now()), last_request_(Clock::now()) {}

    std::string handle(const std::string &raw) {
      const Message req = Message::decode(raw);
      const std::string op = req.get("op", "solve");
      last_request_.store(Clock::now());
      if (op == "ping") {
        Message m;
        m.set("status", "ok");
        return m.encode();
      }
      if (op == "status") return status().encode();
      if (op == "shutdown") {
        Message m;
        m.set("status", "ok");
        m.set("msg", "Shutting down.");
        shutdown_requested_ = true;
        if (server_) server_->stop();
        return m.encode();
      }
      if (op == "solve") return solve(req).encode();
      Message m;
      m.set("status", "error");
      m.set("msg", "Unknown request: " + op);
      return m.encode();
    }

    void set_server(astap::ipc::Server *s) { server_ = s; }
    bool shutdown_requested() const { return shutdown_requested_; }
    Clock::time_point last_request() const { return last_request_.load(); }

  private:
    Message solve(const Message &req) {
      astap::SolveRequest r;
      r.filename = req.get("file");
      r.output_base = req.get("out");
      r.fov = req.number("fov", 0);
      r.max_stars = req.integer("maxstars", 500);
      r.min_star_size = req.number("minstar", 1.5);
      r.downsample = req.integer("downsample", 0);
      r.write_wcs = req.flag("wcs", false);
      r.want_sip = req.flag("sip", false);
      r.refine = req.flag("refine", true);
      r.cmdline = req.get("cmdline");
      const bool progress = req.flag("progress", false);

      Message m;
      if (r.filename.empty()) {
        m.set("status", "error");
        m.set("msg", "No image file in the request.");
        return m;
      }

      // One solve at a time. The star database the second pass reads has a file
      // handle and a tile cache per instance, and the solver already spreads one
      // image across every core, so overlapping two would cost more than it
      // saved even if it were safe.
      const std::lock_guard<std::mutex> lock(solve_mutex_);
      astap::LogFn report;
      if (progress) report = [this](const std::string &s) { log(s); };
      const astap::SolveOutcome out = service_.solve(r, report);

      solves_++;
      recent_ms_.push_back(out.total_seconds * 1000);
      if (recent_ms_.size() > 100) recent_ms_.erase(recent_ms_.begin());

      m.set("status", "ok");
      m.set("solved", out.solved ? 1 : 0);
      m.set("errorlevel", out.errorlevel);
      m.set("solve_seconds", out.solve_seconds, 4);
      m.set("total_seconds", out.total_seconds, 4);
      m.set("stars", out.stars);
      m.set("binning", out.bin);
      if (out.solved) {
        m.set("ra", out.head.ra0, 9);
        m.set("dec", out.head.dec0, 9);
        m.set("scale_arcsec_px", std::fabs(out.head.cdelt2) * 3600, 4);
        m.set("rotation", out.head.crota2, 4);
        m.set("quads", out.nr_inliers);
        m.set("tier", out.tier_density, 1);
        m.set("refined", out.refined ? 1 : 0);
        m.set("sip", out.sip.valid ? 1 : 0);
      }
      for (const std::string &line : out.messages) m.set("msg", line);
      if (!progress)
        for (const std::string &line : out.messages) log(line);
      return m;
    }

    Message status() const {
      Message m;
      m.set("status", "ok");
      m.set("ready", service_.ready() ? 1 : 0);
      m.set("database", service_.database_name());
      m.set("database_path", service_.database_path());
      m.set("cache", service_.cache_path());
      m.set("tiers", static_cast<long long>(service_.tier_count()));
      m.set("tiers_cached", static_cast<long long>(service_.tiers_from_cache()));
      m.set("tiers_built", static_cast<long long>(service_.tiers_built()));
      m.set("quads", static_cast<long long>(service_.quad_count()));
      m.set("bytes", static_cast<long long>(service_.bytes()));
      std::string list;
      for (double d : service_.densities()) {
        if (!list.empty()) list += ",";
        list += astap::float_to_str(d, 1);
      }
      m.set("densities", list);
      m.set("solves", solves_);
      m.set("uptime_seconds", secs(started_, Clock::now()), 1);
      m.set("endpoint", server_ ? server_->endpoint() : std::string());
#ifdef _WIN32
      m.set("pid", static_cast<long long>(GetCurrentProcessId()));
#else
      m.set("pid", static_cast<long long>(getpid()));
#endif
      if (!recent_ms_.empty()) {
        std::vector<double> sorted = recent_ms_;
        std::sort(sorted.begin(), sorted.end());
        m.set("median_ms", sorted[sorted.size() / 2], 1);
        m.set("last_ms", recent_ms_.back(), 1);
      }
      return m;
    }

    void log(const std::string &s) const {
      std::cout << s << std::endl;
      if (config_.logfile.empty()) return;
      std::ofstream f(config_.logfile, std::ios::app);
      if (f) f << s << "\n";
    }

    astap::SolveService &service_;
    Config config_;
    astap::ipc::Server *server_ = nullptr;
    std::mutex solve_mutex_;
    Clock::time_point started_;
    std::atomic<Clock::time_point> last_request_;
    long long solves_ = 0;
    std::vector<double> recent_ms_;
    std::atomic<bool> shutdown_requested_{false};
  };

  // --- outliving nobody ------------------------------------------------------
  //
  // This process holds gigabytes and has no window, so the one thing it must
  // never do is survive the application it was started for. Being asked to stop
  // covers the ordinary exit; it does not cover a crash, a kill, or a machine
  // where the application went away without running its shutdown. So the server
  // also watches the process that started it, and goes when that goes.

  // True only when the process was seen to end. "Could not watch it" has to be
  // told apart from "it ended", or a server that is merely unable to open a
  // handle shuts itself down while the application it belongs to is running
  // perfectly well — the failure would look exactly like the thing it is
  // guarding against.
#ifdef _WIN32
  bool wait_for_process_exit(unsigned long pid, const std::atomic<bool> &give_up,
                             std::string *error) {
    HANDLE h = OpenProcess(SYNCHRONIZE, FALSE, static_cast<DWORD>(pid));
    if (!h) {
      const DWORD e = GetLastError();
      // Gone already is a real answer; anything else means this watch cannot be
      // kept, and the server carries on unwatched rather than quitting.
      if (e == ERROR_INVALID_PARAMETER) return true;
      if (error) *error = "OpenProcess error " + std::to_string(e);
      return false;
    }
    bool exited = false;
    while (!give_up) {
      if (WaitForSingleObject(h, 1000) == WAIT_OBJECT_0) {
        exited = true;
        break;
      }
    }
    CloseHandle(h);
    return exited;
  }
#else
  bool wait_for_process_exit(unsigned long pid, const std::atomic<bool> &give_up,
                             std::string *error) {
    while (!give_up) {
      if (kill(static_cast<pid_t>(pid), 0) != 0) {
        if (errno == ESRCH) return true;
        if (error) *error = std::strerror(errno);
        return false;
      }
      std::this_thread::sleep_for(std::chrono::seconds(1));
    }
    return false;
  }
#endif

  // Ctrl+C and a closing console window, for a server started by hand. Without
  // this the process is torn down where it stands, which loses nothing on
  // Windows but does leave the socket file behind on the platforms that have
  // one.
  astap::ipc::Server *g_server_for_signal = nullptr;

#ifdef _WIN32
  BOOL WINAPI console_handler(DWORD) {
    if (g_server_for_signal) g_server_for_signal->stop();
    return TRUE;
  }
#endif

  int run_server(const Config &config, const std::string &endpoint, bool quiet,
                 unsigned long parent_pid) {
    std::ofstream logstream;
    if (!config.logfile.empty()) logstream.open(config.logfile, std::ios::app);
    auto log = [&](const std::string &s) {
      if (!quiet) std::cout << s << std::endl;
      if (logstream) logstream << s << std::endl;
    };

    if (config.threads > 0) astap::set_thread_count(static_cast<unsigned>(config.threads));

    astap::SolveServiceSettings ss;
    ss.database_path = config.database_path;
    ss.database = config.database;
    ss.quad_tolerance = config.tolerance;
    ss.ladder = astap::resolve_ladder(config.tiers, config.max_tier);
    ss.index_cache = config.cache;

    astap::SolveService service;
    const auto t0 = Clock::now();
    if (!service.load(ss, log)) {
      log("The index could not be loaded, so there is nothing to serve.");
      return astap::kErrNoStarDatabase;
    }
    log("Resident: " + std::to_string(service.tier_count()) + " tiers, " +
        astap::float_to_str(service.bytes() / 1e9, 2) + " GB, ready in " +
        astap::float_to_str(secs(t0, Clock::now()), 2) + " sec.");

    ServerApp app(service, config);
    astap::ipc::Server server;
    std::string err;
    // Listening only now is deliberate: a client that finds the endpoint open
    // knows the answer is milliseconds away, and one that finds it closed can
    // wait on its connect timeout while the ladder is still being read.
    if (!server.start(endpoint, 4, [&app](const std::string &r) { return app.handle(r); }, &err)) {
      log("Cannot listen on " + endpoint + ": " + err);
      return 1;
    }
    app.set_server(&server);
    g_server_for_signal = &server;
#ifdef _WIN32
    SetConsoleCtrlHandler(console_handler, TRUE);
#endif
    log("Listening on " + endpoint + ". Solve with: " + executable_path() + " -f <image>");

    // Tied to the process that asked for this one. It ends when that ends, by
    // whatever means it ended.
    std::atomic<bool> stop_watching{false};
    std::thread parent_watch;
    if (parent_pid) {
      parent_watch = std::thread([&] {
        std::string why;
        const bool exited = wait_for_process_exit(parent_pid, stop_watching, &why);
        if (exited && !stop_watching) {
          log("The process that started this one (" + std::to_string(parent_pid) +
              ") has gone, exiting.");
          server.stop();
        } else if (!exited && !stop_watching && !why.empty()) {
          log("Cannot watch process " + std::to_string(parent_pid) + " (" + why +
              "), carrying on without it.");
        }
      });
    }

    // A server nobody remembers holds gigabytes for nothing. Off by default,
    // because the usual owner of this process is a plugin that stops it itself.
    std::thread idle_watch;
    if (config.idle_exit_minutes > 0) {
      idle_watch = std::thread([&] {
        const auto limit = std::chrono::minutes(config.idle_exit_minutes);
        while (!app.shutdown_requested()) {
          std::this_thread::sleep_for(std::chrono::seconds(10));
          if (app.shutdown_requested()) break;
          if (Clock::now() - app.last_request() > limit) {
            log("Idle for " + std::to_string(config.idle_exit_minutes) + " minutes, exiting.");
            server.stop();
            break;
          }
        }
      });
    }

    server.run();
    stop_watching = true;
    if (idle_watch.joinable()) idle_watch.join();
    // The watcher is parked on the parent's handle for up to a second at a time,
    // so it costs at most that to come back and see that it is finished.
    if (parent_watch.joinable()) parent_watch.join();
    g_server_for_signal = nullptr;
    log("Server stopped, memory released.");
    return 0;
  }

  // --- the client ------------------------------------------------------------

  void print_usage() {
    std::cout
        << "ASTAP astrometric solver, index method, resident (C++)\n"
        "Original algorithm (C) 2018-2026 by Han Kleijn. License MPL 2.0, www.hnsky.org\n"
        "\n"
        "Holds the quad index in memory so that a solve costs milliseconds instead of\n"
        "the 1.6 seconds it takes to read the index back from disk. One executable in\n"
        "two roles: a server that keeps the index, and a client that solves through it\n"
        "while taking the options ASTAP takes.\n"
        "\n"
        "Server:\n"
        "-serve {hold the index and answer requests until stopped}\n"
        "-stop {ask a running server to exit}\n"
        "-status {report what a running server is holding}\n"
        "-config file {settings file, default faster-astap.ini next to this program}\n"
        "-endpoint name {named pipe or socket to use, for more than one server}\n"
        "-parent pid {exit when this process does, however it goes}\n"
        "-quiet {do not write progress to the console}\n"
        "\n"
        "Solving (the options an imaging application passes to ASTAP):\n"
        "-f   filename {image file}\n"
        "-fov diameter_field[degrees] {orders the depth sweep, does not restrict it}\n"
        "-z   downsample_factor[0,1,2,3,4,..] {0 for auto selection}\n"
        "-s   max_number_of_stars {default 500}\n"
        "-m   minimum_star_size[\"] {default 1.5, applied only with -fov}\n"
        "-r, -ra, -spd {accepted and ignored: the index solver needs no start position}\n"
        "-o   file {name the output files with this base path & file name}\n"
        "-wcs {write a .wcs file in the FITS header format}\n"
        "-sip {add SIP distortion coefficients; needs -wcs to be written out}\n"
        "-norefine {skip the second pass, leaving the index solution as it is}\n"
        "-log {write the solver log to a .log text file}\n"
        "-progress {log all progress steps and messages}\n"
        "-noautostart {fail rather than start a server when none is running}\n"
        "-nofallback {fail rather than solve in this process when no server can be had}\n"
        "\n"
        "Settings a fixed option set cannot carry are read from faster-astap.ini next\n"
        "to this program:\n"
        "  database = C:\\Program Files\\astap    directory holding the star database\n"
        "  tiers    = 60,125,250,500             depth ladder, empty for all twelve\n"
        "  maxtier  = 3600                       deepen the default ladder instead\n"
        "  tolerance, cache, endpoint, logfile, threads, autostart, idle_exit_minutes\n"
        "\n"
     << "Image files read by this build: " << astap::supported_image_extensions() << "\n"
     << "\n"
        "Exit status: 0 no errors, 1 no solution, 2 not enough stars detected,\n"
        "16 error reading the image file, 32 no star database found,\n"
        "33 error reading the star database.\n";
  }
} // namespace

int main(int argc, char **argv) {
  std::map<std::string, std::string> opt;
  std::vector<std::string> images;
  std::string cmdline;
  for (int i = 0; i < argc; i++) {
    if (i) cmdline += " ";
    cmdline += argv[i];
  }
  // Which option takes a value is stated rather than guessed, for the reason
  // astap_index_solve states it: guessing "the next token unless it starts with
  // a dash" reads `-progress a.fits b.fits` as -progress=a.fits and silently
  // drops an image.
  //
  // Where this parser has to differ is what it does with an option it does not
  // know. The command line front end can call that a flag and warn, because a
  // person typed it and will read the warning. This program is launched by
  // another program, with an option set that belongs to that program and can
  // gain an entry in any release. Treating an unknown option as a flag would
  // leave its value standing as a bare token, and a bare token here means an
  // image to solve — so a future `-newthing 30` would have this trying to solve
  // a file called "30". An unknown option therefore still swallows a following
  // value, which is the conservative reading when the caller is a machine.
  auto in_list = [](const std::string &k, std::initializer_list<const char *> names) {
    for (const char *n : names)
      if (k == n) return true;
    return false;
  };
  // -r, -ra and -spd are ASTAP's search start position. They are listed because
  // N.I.N.A. passes them and they take a value; the index solver needs no start
  // position, so nothing reads them.
  auto takes_value = [&](const std::string &k) {
    return in_list(k, {"f", "d", "D", "fov", "s", "t", "m", "z", "o", "i", "tiers", "maxtier",
                       "threads", "config", "endpoint", "idle-exit", "parent", "r", "ra", "spd"});
  };
  auto is_flag = [&](const std::string &k) {
    return in_list(k, {"serve", "stop", "status", "quiet", "wcs", "sip", "norefine", "log",
                       "progress", "noautostart", "nofallback", "rebuild", "nocache", "h",
                       "-help"});
  };
  // A token starting with a dash is a value only when it is a negative number.
  auto looks_like_value = [](const char *s) {
    return s[0] != '-' ||
           (s[1] != '\0' && (std::isdigit(static_cast<unsigned char>(s[1])) || s[1] == '.'));
  };

  for (int i = 1; i < argc; i++) {
    std::string a = argv[i];
    if (a.empty() || a[0] != '-') {
      images.push_back(a);
      continue;
    }
    const std::string key = a.substr(1);
    const bool known_flag = is_flag(key);
    if (known_flag) {
      opt[key] = "";
      continue;
    }
    const bool wants = takes_value(key);
    if (i + 1 >= argc || !looks_like_value(argv[i + 1])) {
      if (wants) std::cerr << "Option " << a << " needs a value, ignoring it.\n";
      opt[key] = "";
      continue;
    }
    if (!wants) std::cerr << "Ignoring unknown option " << a << " and its value.\n";
    if (key == "f") images.push_back(argv[++i]);
    else opt[key] = argv[++i];
  }
  auto has = [&](const char *k) { return opt.count(k) != 0; };
  auto val = [&](const char *k) { return opt[k]; };

  if (has("h") || has("-help") || (argc == 1)) {
    print_usage();
    return 0;
  }

  const std::string exe = executable_path();
  const std::string config_path =
      has("config") ? val("config") : directory_of(exe) + "/faster-astap.ini";
  Config config = read_config(config_path);
  // A database directory the settings file does not give is looked for beside
  // the executable, which is where astap_index_solve looks too.
  if (config.database_path.empty()) config.database_path = directory_of(exe);
  if (has("d")) config.database_path = val("d");
  if (has("D")) config.database = val("D");
  if (has("tiers")) config.tiers = val("tiers");
  if (has("t")) config.tolerance = std::atof(val("t").c_str());
  if (has("i")) config.cache = val("i");
  if (has("threads")) config.threads = std::atoi(val("threads").c_str());
  if (has("idle-exit")) config.idle_exit_minutes = std::atoi(val("idle-exit").c_str());
  if (has("noautostart")) config.autostart = false;

  const std::string endpoint = has("endpoint")  ? val("endpoint")
                               : !config.endpoint.empty() ? config.endpoint
                                                          : astap::ipc::default_endpoint();

  if (has("serve"))
    return run_server(config, endpoint, has("quiet"),
                      has("parent") ? std::strtoul(val("parent").c_str(), nullptr, 10) : 0);

  // --- the small requests ----------------------------------------------------
  if (has("stop") || has("status")) {
    Message req;
    req.set("op", has("stop") ? "shutdown" : "status");
    std::string reply, err;
    if (!astap::ipc::request(endpoint, req.encode(), reply, 200, 30000, &err)) {
      std::cout << "No server on " << endpoint << " (" << err << ")." << std::endl;
      return 1;
    }
    const Message m = Message::decode(reply);
    if (has("stop")) {
      std::cout << "Server on " << endpoint << " asked to stop." << std::endl;
      return 0;
    }
    std::cout << "endpoint:  " << endpoint << "\n"
              << "pid:       " << m.get("pid") << "\n"
              << "database:  " << m.get("database") << " in " << m.get("database_path") << "\n"
              << "cache:     " << m.get("cache") << "\n"
              << "resident:  " << m.get("tiers") << " tiers, " << m.get("quads") << " quads, "
              << astap::float_to_str(m.number("bytes", 0) / 1e9, 2) << " GB\n"
              << "tiers:     " << m.get("densities") << "\n"
              << "solves:    " << m.get("solves");
    if (!m.get("median_ms").empty())
      std::cout << ", median " << m.get("median_ms") << " ms, last " << m.get("last_ms") << " ms";
    std::cout << "\nuptime:    " << m.get("uptime_seconds") << " sec" << std::endl;
    return 0;
  }

  if (images.empty()) {
    print_usage();
    return 0;
  }

  // --- solve -----------------------------------------------------------------
  const bool want_log = has("log");
  std::vector<std::string> log_lines;
  auto say = [&](const std::string &s) {
    std::cout << s << std::endl;
    if (want_log) log_lines.push_back(s);
  };

  Message req;
  req.set("op", "solve");
  req.set("out", has("o") ? val("o") : "");
  req.set("fov", has("fov") ? val("fov") : "0");
  req.set("maxstars", has("s") ? val("s") : "500");
  req.set("minstar", has("m") ? val("m") : "1.5");
  req.set("downsample", has("z") ? val("z") : "0");
  req.set("wcs", has("wcs") ? 1 : 0);
  req.set("sip", has("sip") ? 1 : 0);
  req.set("refine", has("norefine") ? 0 : 1);
  req.set("progress", has("progress") ? 1 : 0);
  req.set("cmdline", cmdline);

  int worst = 0;
  bool server_gone = false;
  astap::SolveService fallback;  // loaded at most once, and only if it is needed

  for (const std::string &filename : images) {
    Message one = req;
    one.set("file", filename);

    std::string reply, err;
    bool answered = false;
    if (!server_gone) {
      // A short connect timeout on the first try: either a server is already
      // there, or there is no point waiting for one that was never started.
      answered = astap::ipc::request(endpoint, one.encode(), reply, 200, 600000, &err);
      if (!answered && config.autostart) {
        std::string spawn_err;
        if (spawn_server(exe, config_path, &spawn_err)) {
          // The first solve after a start waits for the ladder to be read, which
          // is the 1.6 s this whole program exists to stop paying repeatedly.
          answered = astap::ipc::request(endpoint, one.encode(), reply, 120000, 600000, &err);
        } else {
          err = spawn_err;
        }
      }
      if (!answered) server_gone = true;
    }

    if (answered) {
      const Message m = Message::decode(reply);
      for (const std::string &line : m.all("msg")) say(line);
      if (m.get("status") != "ok") {
        worst = std::max(worst, 1);
        continue;
      }
      worst = std::max(worst, m.integer("errorlevel", 0));
      continue;
    }

    // No server, and none could be started. Solving here costs the index read
    // this program exists to avoid, but a slow solve beats a lost frame in the
    // middle of a night.
    if (has("nofallback")) {
      say("No index server on " + endpoint + " (" + err + ") and no fallback allowed.");
      astap::Header head;
      astap::write_ini(astap::change_file_ext(has("o") ? val("o") : filename, ".ini"), false, head,
                       cmdline, astap::kErrNoStarDatabase, "");
      worst = std::max(worst, static_cast<int>(astap::kErrNoStarDatabase));
      continue;
    }
    if (!fallback.ready()) {
      say("No index server on " + endpoint + " (" + err + "), solving in this process.");
      if (config.threads > 0) astap::set_thread_count(static_cast<unsigned>(config.threads));
      astap::SolveServiceSettings ss;
      ss.database_path = config.database_path;
      ss.database = config.database;
      ss.quad_tolerance = config.tolerance;
      ss.ladder = astap::resolve_ladder(config.tiers, config.max_tier);
      ss.index_cache = config.cache;
      if (!fallback.load(ss, say)) {
        astap::Header head;
        astap::write_ini(astap::change_file_ext(has("o") ? val("o") : filename, ".ini"), false,
                         head, cmdline, astap::kErrNoStarDatabase, "");
        worst = std::max(worst, static_cast<int>(astap::kErrNoStarDatabase));
        continue;
      }
    }

    astap::SolveRequest r;
    r.filename = filename;
    r.output_base = has("o") ? val("o") : "";
    r.fov = has("fov") ? std::atof(val("fov").c_str()) : 0;
    r.max_stars = has("s") ? std::atoi(val("s").c_str()) : 500;
    r.min_star_size = has("m") ? std::atof(val("m").c_str()) : 1.5;
    r.downsample = has("z") ? std::atoi(val("z").c_str()) : 0;
    r.write_wcs = has("wcs");
    r.want_sip = has("sip");
    r.refine = !has("norefine");
    r.cmdline = cmdline;
    astap::LogFn report;
    if (has("progress")) report = say;
    const astap::SolveOutcome out = fallback.solve(r, report);
    if (!has("progress"))
      for (const std::string &line : out.messages) say(line);
    worst = std::max(worst, out.errorlevel);
  }

  if (want_log) {
    std::ofstream lf(
        astap::change_file_ext(has("o") ? val("o") : images.front(), ".log"));
    lf << cmdline << "\n";
    for (const std::string &l : log_lines) lf << l << "\n";
  }
  return worst;
}
