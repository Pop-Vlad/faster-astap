// The index solver behind an ASTAP command line, for an imaging application.
//
// N.I.N.A. 3.x has no extension point for plate solvers, but its ASTAP solver
// is a CLI solver: it launches whatever executable PlateSolveSettings
// .ASTAPLocation names, passes ASTAP's options, and reads the .ini written
// beside the image. This program takes that option set and writes that .ini, so
// pointing the setting at it is the whole integration.
//
// It is astap_index_solve with two differences, both of them about being
// launched by a program rather than typed by a person. It accepts and ignores
// the options ASTAP takes that an index solve has no use for (-r, -ra, -spd),
// and the settings a fixed option set cannot carry - where the star database
// lives, which depth tiers to use - come from faster-astap.ini next to the
// executable, which the plugin writes and a person can edit.

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <initializer_list>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "astap/astro_math.h"
#include "astap/image/image_io.h"
#include "astap/parallel.h"
#include "astap/solve_service.h"
#include "astap/solver.h"

#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif

namespace {
  // --- where the executable is ----------------------------------------------
  //
  // The settings file sits beside it, not in the working directory: the working
  // directory of a process launched by an imaging application is whatever that
  // application happened to be in.

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

  struct Config {
    std::string database_path;
    std::string database = "auto";
    double tolerance = 0.007;
    std::string tiers; // empty selects the default ladder
    double max_tier = 0; // raises the ceiling of the default ladder
    std::string cache;
    std::string logfile;
    int threads = 0;
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
      else if (k == "logfile") c.logfile = v;
      else if (k == "threads") c.threads = std::atoi(v.c_str());
    }
    return c;
  }

  void print_usage() {
    std::cout
        << "ASTAP astrometric solver, index method, ASTAP compatible front end (C++)\n"
        "Original algorithm (C) 2018-2026 by Han Kleijn. License MPL 2.0, www.hnsky.org\n"
        "\n"
        "Solves against a pre-built whole-sky quad index while taking the options an\n"
        "imaging application passes to ASTAP, and writing the .ini it reads back. The\n"
        "index cache is memory mapped, so the rungs one solve touches stay in the page\n"
        "cache for the solves after it.\n"
        "\n"
        "Solving (the options an imaging application passes to ASTAP):\n"
        "-f   filename {image file. May be repeated, or list files after the options}\n"
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
        "-prepare {build or read the index and exit, solving nothing}\n"
        "\n"
        "Settings a fixed option set cannot carry are read from faster-astap.ini next\n"
        "to this program, and can be overridden on the command line:\n"
        "  database = C:\\Program Files\\astap    directory holding the star database\n"
        "  tiers    = 60,125,250,500             depth ladder, empty for all twelve\n"
        "  maxtier  = 3600                       deepen the default ladder instead\n"
        "  tolerance, cache, logfile, threads\n"
        "-config file {settings file, default faster-astap.ini next to this program}\n"
        "-d path, -D abbreviation, -t tolerance, -i cache, -tiers, -maxtier, -threads\n"
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
  // image to solve - so a future `-newthing 30` would have this trying to solve
  // a file called "30". An unknown option therefore still swallows a following
  // value, which is the conservative reading when the caller is a machine.
  auto in_list = [](const std::string &k, std::initializer_list<const char *> names) {
    for (const char *n: names)
      if (k == n) return true;
    return false;
  };
  // -r, -ra and -spd are ASTAP's search start position. They are listed because
  // N.I.N.A. passes them and they take a value; the index solver needs no start
  // position, so nothing reads them.
  auto takes_value = [&](const std::string &k) {
    return in_list(k, {
                     "f", "d", "D", "fov", "s", "t", "m", "z", "o", "i", "tiers", "maxtier",
                     "threads", "config", "r", "ra", "spd"
                   });
  };
  auto is_flag = [&](const std::string &k) {
    return in_list(k, {
                     "wcs", "sip", "norefine", "log", "progress", "prepare", "rebuild",
                     "nocache", "h", "-help"
                   });
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
    if (is_flag(key)) {
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

  if (has("h") || has("-help") || argc == 1) {
    print_usage();
    return 0;
  }

  const std::string exe = executable_path();
  const std::string config_path =
      has("config") ? val("config") : directory_of(exe) + "/faster-astap.ini";
  Config config = read_config(config_path);
  // A database directory that neither the settings file nor the command line
  // gives is left empty, which sends the service to the usual places: beside
  // this executable, then wherever ASTAP installed one.
  if (has("d")) config.database_path = val("d");
  if (has("D")) config.database = val("D");
  if (has("t")) config.tolerance = std::atof(val("t").c_str());
  if (has("tiers")) config.tiers = val("tiers");
  if (has("maxtier")) config.max_tier = std::atof(val("maxtier").c_str());
  if (has("i")) config.cache = val("i");
  if (has("threads")) config.threads = std::atoi(val("threads").c_str());

  if (!has("prepare") && images.empty()) {
    print_usage();
    return 0;
  }

  // Whatever is said goes to the console, to the .log file when -log asked for
  // one, and to the settings file's logfile when it names one. That last is the
  // only channel that reaches anybody here: a process launched by an imaging
  // application has a console nobody sees.
  const bool want_log = has("log");
  std::vector<std::string> log_lines;
  std::ofstream logstream;
  if (!config.logfile.empty()) logstream.open(config.logfile, std::ios::app);
  astap::LogFn say = [&](const std::string &s) {
    std::cout << s << std::endl;
    if (want_log) log_lines.push_back(s);
    if (logstream) logstream << s << std::endl;
  };

  if (config.threads > 0) astap::set_thread_count(static_cast<unsigned>(config.threads));

  astap::SolveServiceSettings ss;
  ss.database_path = config.database_path;
  ss.database = config.database;
  ss.quad_tolerance = config.tolerance;
  ss.ladder = astap::resolve_ladder(config.tiers, config.max_tier);
  ss.index_cache = config.cache;
  ss.use_cache = !has("nocache");
  ss.rebuild = has("rebuild");

  // --- get an index ----------------------------------------------------------
  astap::SolveService service;
  if (!service.load(ss, say)) {
    // The .ini is written even for this, because the application waiting on it
    // reads a missing file as a solver that never ran at all.
    for (const std::string &f: images) {
      astap::Header head;
      astap::write_ini(astap::change_file_ext(has("o") ? val("o") : f, ".ini"), false, head,
                       cmdline, astap::kErrNoStarDatabase, "");
    }
    return astap::kErrNoStarDatabase;
  }

  // Building a ladder that has never been built takes minutes, and the point of
  // -prepare is to pay that at a chosen moment rather than on the first frame of
  // a night. Reading a cached one takes no time worth choosing a moment for.
  if (has("prepare")) {
    say("Index ready: " + std::to_string(service.tier_count()) + " tiers, " +
        astap::float_to_str(service.bytes() / 1e9, 2) + " GB, cached in " + service.cache_path());
    if (images.empty()) return 0;
  }

  // --- solve ------------------------------------------------------------------
  int worst = 0;
  for (const std::string &filename: images) {
    astap::SolveRequest req;
    req.filename = filename;
    req.output_base = has("o") ? val("o") : "";
    req.params.fov = has("fov") ? std::atof(val("fov").c_str()) : 0;
    req.params.max_stars = has("s") ? std::atoi(val("s").c_str()) : 500;
    req.params.min_star_size = has("m") ? std::atof(val("m").c_str()) : 1.5;
    req.params.downsample = has("z") ? std::atoi(val("z").c_str()) : 0;
    req.write_wcs = has("wcs");
    req.params.want_sip = has("sip");
    req.params.refine = !has("norefine");
    req.cmdline = cmdline;

    // With -progress the service reports as it goes; without it there is nothing
    // to interleave, so the summary is printed once the image is done.
    const astap::SolveOutcome out = service.solve(req, has("progress") ? say : astap::LogFn());
    if (!has("progress"))
      for (const std::string &m: out.messages) say(m);

    worst = std::max(worst, out.errorlevel);
  }

  if (want_log) {
    std::ofstream lf(astap::change_file_ext(has("o") ? val("o") : images.front(), ".log"));
    lf << cmdline << "\n";
    for (const std::string &l: log_lines) lf << l << "\n";
  }
  return worst;
}
