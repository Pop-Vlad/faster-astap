// Command line front end, modelled on astap_command_line.lpr.

#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "astap/astro_math.h"
#include "astap/image/fits.h"
#include "astap/image/image_io.h"
#include "astap/parallel.h"
#include "astap/solver.h"

namespace {
  std::string change_file_ext(const std::string &path, const std::string &ext) {
    size_t slash = path.find_last_of("/\\");
    size_t dot = path.find_last_of('.');
    if (dot == std::string::npos || (slash != std::string::npos && dot < slash)) return path + ext;
    return path.substr(0, dot) + ext;
  }

  std::string dir_of(const std::string &path) {
    size_t slash = path.find_last_of("/\\");
    return slash == std::string::npos ? std::string("./") : path.substr(0, slash + 1);
  }

  // The database path is concatenated with a bare file name, so it has to end in
  // a separator. Windows takes '/' too, and a path the user already ended with
  // '\' is left as it is.
  std::string with_separator(std::string dir) {
    if (!dir.empty() && dir.back() != '/' && dir.back() != '\\') dir += '/';
    return dir;
  }

  void print_usage() {
    std::cout
        << "ASTAP astrometric solver, C++ port\n"
        "Original algorithm (C) 2018-2026 by Han Kleijn. License MPL 2.0, www.hnsky.org\n"
        "Usage:\n"
        "-f   filename {image file, see the list of extensions below}\n"
        "-r   radius_area_to_search[degrees]\n"
        "-fov diameter_field[degrees] {enter zero for auto}\n"
        "-ra  right_ascension[hours]\n"
        "-spd south_pole_distance[degrees]\n"
        "-s   max_number_of_stars {default 500}\n"
        "-t   quad_tolerance {default 0.007}\n"
        "-m   minimum_star_size[\"] {default 1.5}\n"
        "-z   downsample_factor[0,1,2,3,4,..] {0 for auto selection}\n"
        "-check {apply the check pattern filter, raw OSC images with binning 1x1 only}\n"
        "-d   path {path to the star database}\n"
        "-D   abbreviation[d80,d50,...] {select a star database}\n"
        "-o   file {name the output files with this base path & file name}\n"
        "-sip {add SIP (Simple Image Polynomial) coefficients}\n"
        "-speed mode[auto/slow] {slow forces more area overlap while searching}\n"
        "-wcs {write a .wcs file in the FITS header format}\n"
        "-log {write the solver log to a .log text file}\n"
        "-progress {log all progress steps and messages}\n"
        "-threads N {worker threads, 0 or omitted = one per available hardware thread}\n"
        "\n"
        "Preference is given to the command line values. The solver result is written to\n"
        "filename.ini and, with -wcs, to filename.wcs.\n"
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
  std::string cmdline;
  for (int i = 0; i < argc; i++) {
    if (i) cmdline += " ";
    cmdline += argv[i];
  }
  for (int i = 1; i < argc; i++) {
    std::string a = argv[i];
    if (a.empty() || a[0] != '-') continue;
    std::string key = a.substr(1);
    if (i + 1 < argc && argv[i + 1][0] != '-') {
      opt[key] = argv[i + 1];
      i++;
    } else {
      opt[key] = ""; // flag without a value
    }
  }

  auto has = [&](const char *k) { return opt.count(k) != 0; };
  auto val = [&](const char *k) { return opt[k]; };

  if (!has("f")) {
    print_usage();
    return 0;
  }

  astap::SolverSettings settings;
  // Without -d, look where a database is likely to already be: beside this
  // executable, then wherever ASTAP installed one. The solver takes a single
  // directory, so the choice is made here; the first that holds a readable
  // database wins, and if none does the original path is left in place so the
  // solver reports the failure it would have reported anyway.
  if (has("d")) {
    settings.database_path = with_separator(val("d"));
  } else {
    settings.database_path = with_separator(dir_of(argv[0]));
    astap::StarDatabase probe;
    for (const std::string &dir: astap::default_database_directories())
      if (probe.select(with_separator(dir), has("D") ? val("D") : "auto", 1.0)) {
        settings.database_path = with_separator(dir);
        break;
      }
  }

  const std::string filename = val("f");

  astap::Header head;
  astap::ImageArray img;
  astap::ImageLoadResult loaded = astap::load_image(filename, head, img);

  std::vector<std::string> log_lines;
  const bool want_log = has("log");
  const bool progress = has("progress");

  astap::LogFn log = [&](const std::string &s) {
    std::cout << s << std::endl;
    if (want_log) log_lines.push_back(s);
  };

  if (!loaded.ok) {
    log(loaded.error);
    const std::string out = has("o") ? val("o") : filename;
    astap::write_ini(change_file_ext(out, ".ini"), false, head, cmdline, astap::kErrImageRead, "");
    return astap::kErrImageRead;
  }
  if (!loaded.warning.empty()) log(loaded.warning);

  // Command line values override the header.
  if (has("fov")) {
    settings.fov_specified = true; // do not calculate it from the header
    settings.search_fov = std::atof(val("fov").c_str());
  }
  if (has("r")) settings.radius_search = std::atof(val("r").c_str());
  if (has("ra")) head.ra0 = std::atof(val("ra").c_str()) * astap::kPi / 12;
  // South pole distance, because negative values cannot be passed on the command line.
  if (has("spd")) head.dec0 = (std::atof(val("spd").c_str()) - 90) * astap::kPi / 180;
  if (has("z")) settings.downsample = std::atoi(val("z").c_str());
  if (has("s")) settings.max_stars = std::atoi(val("s").c_str());
  if (has("t")) settings.quad_tolerance = std::atof(val("t").c_str());
  if (has("m")) settings.min_star_size = std::atof(val("m").c_str());
  if (has("sip")) settings.add_sip = val("sip") != "n";
  if (has("speed")) settings.force_oversize = val("speed").find("slow") != std::string::npos;
  if (has("check")) settings.check_pattern_filter = val("check") != "n";
  if (has("D")) settings.star_database = val("D");
  settings.show_log = progress;
  if (has("threads")) astap::set_thread_count((unsigned) std::atoi(val("threads").c_str()));

  const std::string filename_output = has("o") ? val("o") : filename;

  astap::Solver solver(settings);
  solver.set_log(log);
  solver.set_mount(loaded.ra_mount, loaded.dec_mount);

  const bool solution = solver.solve(std::move(img), head);

  {
    const astap::Solver::Timing &t = solver.timing();
    char buf[320];
    std::snprintf(buf, sizeof(buf),
                  "Timing: image stages %.2f s wall, spiral search %.2f s wall "
                  "(database read %.2f s, quad build %.2f s, match %.2f s)",
                  t.image_wall, t.spiral_wall, t.read_stars_cpu, t.quads_cpu, t.match_cpu);
    log(buf);
  }

  astap::write_ini(change_file_ext(filename_output, ".ini"), solution, head, cmdline,
                   solver.errorlevel(), solver.warning());

  if (solution && has("wcs")) {
    std::string comment = "Solved in " + astap::float_to_str(solver.solved_seconds(), 1) +
                          " sec. Offset was " +
                          astap::distance_to_string(solver.search_offset(), solver.search_offset()) +
                          ".";
    astap::update_solution_cards(head.cards, head, solver.sip(), true, comment);
    // A minimal header: the .wcs file carries no image data.
    astap::remove_key(head.cards, "NAXIS1  =");
    astap::remove_key(head.cards, "NAXIS2  =");
    astap::update_integer(head.cards, "NAXIS   =",
                          " / Minimal header                                 ", 0);
    astap::update_integer(head.cards, "BITPIX  =",
                          " /                                                ", 8);
    astap::write_fits_header_file(change_file_ext(filename_output, ".wcs"), head.cards);
  }

  if (want_log) {
    std::ofstream lf(change_file_ext(filename_output, ".log"));
    lf << cmdline << "\n";
    for (const std::string &l: log_lines) lf << l << "\n";
  }

  return solution ? 0 : solver.errorlevel();
}
