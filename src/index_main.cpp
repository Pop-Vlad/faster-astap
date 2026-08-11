// Command line front end for the index solver, alongside astap_solve.
//
// The two binaries take the same options and write the same .ini and .wcs, so
// one can be swapped for the other. What differs is what happens in between:
// astap_solve walks a squared spiral over the sky, rebuilding database quads at
// every position; this one queries a pre-built whole-sky quad index and votes.
//
// The index is a ladder of depth tiers (see README.md), cached one rung per file
// under ~/.cache/faster-astap. Building a rung takes a few seconds and the result
// is the same for every image and every run, so the rungs a run shares with an
// earlier one are read back and only the new ones are built.
//
// The per image work itself lives in SolveService, which is also what the
// N.I.N.A. front end runs, so the two produce the same output for the same
// image.

// <chrono> went with the timing, which moved into SolveService; <cctype> stays,
// the option parser below needs isdigit.
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <initializer_list>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "astap/image/image_io.h"
#include "astap/parallel.h"
#include "astap/quad_index.h"
#include "astap/solve_service.h"
#include "astap/solver.h"
#include "astap/star_database.h"

namespace {
  std::string dir_of(const std::string &path) {
    size_t slash = path.find_last_of("/\\");
    return slash == std::string::npos ? std::string("./") : path.substr(0, slash + 1);
  }

  void print_usage() {
    std::cout
        << "ASTAP astrometric solver, index method (C++)\n"
        "Original algorithm (C) 2018-2026 by Han Kleijn. License MPL 2.0, www.hnsky.org\n"
        "\n"
        "Solves against a pre-built whole-sky quad index instead of a spiral search.\n"
        "The index is cached under ~/.cache/faster-astap and built on first use.\n"
        "\n"
        "Usage:\n"
        "-f   filename {image file. May be repeated, or list files after the options}\n"
        "-d   path {path to the star database, needed to build an index}\n"
        "-D   abbreviation[d80,d50,...] {select a star database}\n"
        "-fov diameter_field[degrees] {orders the depth sweep, does not restrict it}\n"
        "-s   max_number_of_stars {default 500}\n"
        "-t   quad_tolerance {default 0.007. A cache belongs to one tolerance}\n"
        "-m   minimum_star_size[\"] {default 1.5, applied only with -fov}\n"
        "-z   downsample_factor[0,1,2,3,4,..] {0 for auto selection}\n"
        "-o   file {name the output files with this base path & file name}\n"
        "-sip {add SIP distortion coefficients; needs -wcs to be written out}\n"
        "-norefine {skip the second pass, leaving the index solution as it is}\n"
        "-wcs {write a .wcs file in the FITS header format}\n"
        "-log {write the solver log to a .log text file}\n"
        "-progress {log all progress steps and messages}\n"
        "-threads N {worker threads, 0 or omitted = one per available hardware thread}\n"
        "\n"
        "Index options:\n"
        "-i     file {index cache to use, overriding the default location}\n"
        "-tiers d1,d2,.. {depth ladder in stars/deg^2 for a newly built index}\n"
        "-maxtier d {raise the ceiling of the default ladder to d stars/deg^2,\n"
        "            which is worth doing below 0.25 degrees and not above:\n"
        "            3600 takes the 0.15 degree corpus from 3/16 to 10/16 and\n"
        "            adds rungs of 2.5 and 4.8 GB, built and cached once}\n"
        "-rebuild {rebuild the index even when a usable cache exists}\n"
        "-nocache {build in memory, do not read or write a cache}\n"
        "-cacheinfo {report the cache that would be used, then exit}\n"
        "\n"
        "For an imaging application that launches a solver per frame with ASTAP's\n"
        "options, see astap_nina_solve -h.\n"
        "\n"
        "The solver result is written to filename.ini and, with -wcs, to filename.wcs.\n"
        "\n"
        << "Image files read by this build: " << astap::supported_image_extensions() << "\n"
        << "\n"
        "Exit status: 0 no errors, 1 no solution, 2 not enough stars detected,\n"
        "16 error reading the image file, 32 no star database found,\n"
        "33 error reading the star database.\n";
  }
} // namespace

int main(int argc, char **argv) {
  if (argc == 1 || std::strcmp(argv[1], "-h") == 0 || std::strcmp(argv[1], "--help") == 0) {
    print_usage();
    return 0;
  }

  std::map<std::string, std::string> opt;
  std::vector<std::string> images;
  std::string cmdline;
  for (int i = 0; i < argc; i++) {
    if (i) cmdline += " ";
    cmdline += argv[i];
  }
  // Which option takes a value has to be stated, not guessed. Guessing it as
  // "the next token, unless it starts with a dash" reads
  // `astap_index_solve -progress a.fits b.fits` as -progress=a.fits and then
  // solves only b.fits - an image silently dropped, which is the worst way to
  // get this wrong. An unknown option is a flag and says so, so that a typo like
  // -maxtiers does not quietly leave the default ladder in place.
  auto in_list = [](const std::string &k, std::initializer_list<const char *> names) {
    for (const char *n: names)
      if (k == n) return true;
    return false;
  };
  auto takes_value = [&](const std::string &k) {
    return in_list(k, {
                     "f", "d", "D", "fov", "s", "t", "m", "z", "o", "i", "tiers", "maxtier",
                     "threads"
                   });
  };
  auto is_flag = [&](const std::string &k) {
    return in_list(k, {
                     "sip", "norefine", "wcs", "log", "progress", "rebuild", "nocache",
                     "cacheinfo", "h", "-help"
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
      images.push_back(a); // bare file name, so a shell glob can be passed
      continue;
    }
    const std::string key = a.substr(1);
    if (!takes_value(key)) {
      if (!is_flag(key)) std::cerr << "Ignoring unknown option " << a << "\n";
      opt[key] = ""; // flag without a value
      continue;
    }
    if (i + 1 >= argc || !looks_like_value(argv[i + 1])) {
      std::cerr << "Option " << a << " needs a value, ignoring it.\n";
      continue;
    }
    if (key == "f") images.push_back(argv[++i]);
    else opt[key] = argv[++i];
  }

  auto has = [&](const char *k) { return opt.count(k) != 0; };
  auto val = [&](const char *k) { return opt[k]; };

  const bool want_log = has("log");
  const bool progress = has("progress");
  std::vector<std::string> log_lines;
  astap::LogFn log = [&](const std::string &s) {
    std::cout << s << std::endl;
    if (want_log) log_lines.push_back(s);
  };

  if (has("threads")) astap::set_thread_count((unsigned) std::atoi(val("threads").c_str()));

  astap::SolveServiceSettings ss;
  // Empty means "look in the usual places": beside this executable, then where
  // ASTAP installs its databases. -d overrides that outright.
  ss.database_path = has("d") ? val("d") : "";
  ss.database = has("D") ? val("D") : "auto";
  ss.quad_tolerance = has("t") ? std::atof(val("t").c_str()) : 0.007;
  // -tiers names a ladder outright; -maxtier just raises the ceiling of the
  // default one, which is the common case: a narrow field needs deeper rungs and
  // nothing else about the ladder changes.
  ss.ladder = astap::resolve_ladder(has("tiers") ? val("tiers") : "",
                                    has("maxtier") ? std::atof(val("maxtier").c_str()) : 0);
  ss.index_cache = has("i") ? val("i") : "";
  ss.use_cache = !has("nocache");
  ss.rebuild = has("rebuild");

  // --- report the cache and stop ---------------------------------------------
  if (has("cacheinfo")) {
    astap::StarDatabase db;
    // The same search the solve would do, so the report names the database a
    // solve would actually use.
    bool have_db = false;
    if (!ss.database_path.empty()) {
      have_db = db.select(astap::with_separator(ss.database_path), ss.database, 1.0);
    } else {
      for (const std::string &dir: astap::default_database_directories())
        if ((have_db = db.select(astap::with_separator(dir), ss.database, 1.0))) break;
    }
    const std::string db_name = have_db ? db.name() : "unknown";
    const int db_type = have_db ? db.database_type() : 0;
    const std::string ladder_file =
        ss.index_cache.empty()
          ? astap::default_index_cache_path(db_name, db_type, ss.quad_tolerance)
          : ss.index_cache;

    std::cout << "star database: "
        << (have_db ? db.name() + " in " + db.path() : std::string("(none found)"))
        << "\nladder wanted:";
    for (double d: ss.ladder) std::cout << " " << d;
    std::cout << "\n";
    astap::QuadIndexFile info;
    std::string err;
    double total = 0;
    for (double d: ss.ladder) {
      const std::string path =
          ss.index_cache.empty()
            ? astap::index_tier_cache_path(db_name, db_type, ss.quad_tolerance, d)
            : ladder_file;
      if (astap::read_index_file_header(path, info, &err)) {
        std::cout << "  tier " << d << ": " << info.bytes / 1e9 << " GB  " << path << "\n";
        total += info.bytes / 1e9;
      } else {
        std::cout << "  tier " << d << ": not built\n";
      }
      if (!ss.index_cache.empty()) break;
    }
    // The pre-split layout, still read when it is there.
    if (ss.index_cache.empty() && astap::read_index_file_header(ladder_file, info, &err)) {
      std::cout << "  a ladder file from before the tiers were split is also present ("
          << info.densities.size() << " tiers, " << info.bytes / 1e9
          << " GB): " << ladder_file << "\n";
    }
    std::cout << "  total " << total << " GB\n";
    return 0;
  }

  if (images.empty()) {
    print_usage();
    return 0;
  }

  // --- get an index ----------------------------------------------------------
  astap::SolveService service;
  if (!service.load(ss, log)) {
    for (const std::string &f: images) {
      astap::Header head;
      astap::write_ini(astap::change_file_ext(has("o") ? val("o") : f, ".ini"), false, head,
                       cmdline, astap::kErrNoStarDatabase, "");
    }
    return astap::kErrNoStarDatabase;
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
    const astap::SolveOutcome out = service.solve(req, progress ? log : astap::LogFn());
    if (!progress)
      for (const std::string &m: out.messages) log(m);

    worst = std::max(worst, out.errorlevel);
  }

  if (want_log) {
    std::ofstream lf(astap::change_file_ext(has("o") ? val("o") : images.front(), ".log"));
    lf << cmdline << "\n";
    for (const std::string &l: log_lines) lf << l << "\n";
  }
  return worst;
}
