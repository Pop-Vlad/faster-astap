// Command line front end for the index solver, alongside astap_solve.
//
// The two binaries take the same options and write the same .ini and .wcs, so
// one can be swapped for the other. What differs is what happens in between:
// astap_solve walks a squared spiral over the sky, rebuilding database quads at
// every position; this one queries a pre-built whole-sky quad index and votes.
//
// The index is a ladder of depth tiers (see README.md). Building it
// takes a few seconds and the result is the same for every image and every run,
// so it is cached under ~/.cache/faster-astap and reused. The first run on a new
// star database pays for the build; later runs do not.

#include <cctype>
#include <chrono>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <fstream>
#include <initializer_list>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "astap/astro_math.h"
#include "astap/image/fits.h"
#include "astap/image/image_io.h"
#include "astap/index_solver.h"
#include "astap/parallel.h"
#include "astap/quad_index.h"
#include "astap/quads.h"
#include "astap/solver.h"
#include "astap/star_detection.h"

using Clock = std::chrono::steady_clock;

namespace {
  double secs(Clock::time_point a, Clock::time_point b) {
    return std::chrono::duration<double>(b - a).count();
  }

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

  // The default depth ladder. One tier reaches about a factor of two in image
  // star density either side of itself, so the 2x steps below overlap and cover
  // 0.2 to 2000 stars/deg^2 continuously. Measured over the corpus, an image
  // whose density falls between two rungs is solved by one of them.
  const std::vector<double> kDefaultLadder = {0.5, 1, 2, 4, 8, 16, 32, 60, 125, 250, 500, 900};

  // Rungs past the default ceiling, for -maxtier. They are not in the default
  // ladder because they are most of the cache — 900 is 1.2 GB on disk, 1800 adds
  // 2.5 and 3600 another 4.8 — and because the density matching in
  // solve_stars_with_tiers reaches the same depths for free down to about 0.25
  // degrees, which is where the default ladder now stops rather than 0.5. Only
  // smaller fields need these. A deeper rung than 3600 buys nothing measurable:
  // the database runs out of stars and the images that still fail are short of
  // stars, not of depth. See the small fields section of README.md.
  const std::vector<double> kDeepLadder = {1800, 3600};

  std::vector<double> parse_densities(const std::string &v) {
    std::vector<double> out;
    size_t p = 0;
    while (p <= v.size()) {
      const size_t c = v.find(',', p);
      const std::string tok = v.substr(p, c == std::string::npos ? std::string::npos : c - p);
      if (!tok.empty()) out.push_back(std::atof(tok.c_str()));
      if (c == std::string::npos) break;
      p = c + 1;
    }
    return out;
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
        "The solver result is written to filename.ini and, with -wcs, to filename.wcs.\n"
        "\n"
     << "Image files read by this build: " << astap::supported_image_extensions() << "\n"
     << "\n"
        "Exit status: 0 no errors, 1 no solution, 2 not enough stars detected,\n"
        "16 error reading the image file, 32 no star database found,\n"
        "33 error reading the star database.\n";
  }

  // Mono, binned copy of the image, plus the binning factor that produced it.
  //
  // Binning past about 1200 px on the long side throws away the stars the solve
  // depends on, so the automatic factor stops there. The index solver works on
  // this copy, and its solution is scaled back to original pixels afterwards.
  int bin_mono(const astap::ImageArray &img, int width, int height, int forced,
               astap::ImageArray &out) {
    int bin = forced;
    if (bin <= 0) {
      bin = 1;
      while (std::max(width, height) / (bin + 1) >= 1200) bin++;
    }
    out = astap::ImageArray(1, height / bin, width / bin);
    for (int y = 0; y < out.height(); y++)
      for (int x = 0; x < out.width(); x++) {
        double v = 0;
        for (int c = 0; c < img.colours(); c++)
          for (int i = 0; i < bin; i++)
            for (int j = 0; j < bin; j++) v += img.at(c, y * bin + i, x * bin + j);
        out.at(0, y, x) = static_cast<float>(v / (img.colours() * bin * bin));
      }
    return bin;
  }

  // The index solver returns a WCS in binned pixels. Binned pixel xb covers
  // original pixels [xb*bin, xb*bin+bin), so its centre sits at
  // xb*bin + (bin-1)/2 in original zero-based coordinates. Rewriting the
  // solution in those coordinates is a scale of the pixel size and a shift of
  // the reference pixel; the rotation is unchanged.
  void scale_solution_to_original(const astap::IndexSolveResult &r, int bin, astap::Header &head) {
    head.ra0 = r.ra0;
    head.dec0 = r.dec0;
    head.crpix1 = (r.crpix1 - 1) * bin + (bin - 1) / 2.0 + 1;
    head.crpix2 = (r.crpix2 - 1) * bin + (bin - 1) / 2.0 + 1;
    head.cdelt1 = r.cdelt1 / bin;
    head.cdelt2 = r.cdelt2 / bin;
    head.crota1 = r.crota1;
    head.crota2 = r.crota2;
    head.cd1_1 = r.cd1_1 / bin;
    head.cd1_2 = r.cd1_2 / bin;
    head.cd2_1 = r.cd2_1 / bin;
    head.cd2_2 = r.cd2_2 / bin;
  }
  // SIP coefficients from the binned frame, rewritten for original pixels.
  //
  // The polynomial takes an offset from the reference pixel and returns a
  // correction, both in pixels, and offsets scale exactly: u_orig = bin * u_bin,
  // because the reference pixel is mapped by the same relation as every other.
  // A term c*u^p*v^q therefore becomes c * bin^(1-p-q) in the original frame.
  void scale_sip_to_original(astap::SipCoefficients &sip, int bin) {
    if (!sip.valid || bin == 1) return;
    double (*tables[4])[4] = {sip.a, sip.b, sip.ap, sip.bp};
    for (double (*t)[4] : tables)
      for (int p = 0; p < 4; p++)
        for (int q = 0; q < 4; q++)
          if (t[p][q] != 0) t[p][q] *= std::pow(static_cast<double>(bin), 1 - p - q);
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
  // solves only b.fits — an image silently dropped, which is the worst way to
  // get this wrong. An unknown option is a flag and says so, so that a typo like
  // -maxtiers does not quietly leave the default ladder in place.
  auto in_list = [](const std::string &k, std::initializer_list<const char *> names) {
    for (const char *n : names)
      if (k == n) return true;
    return false;
  };
  auto takes_value = [&](const std::string &k) {
    return in_list(k, {"f", "d", "D", "fov", "s", "t", "m", "z", "o", "i", "tiers", "maxtier",
                       "threads"});
  };
  auto is_flag = [&](const std::string &k) {
    return in_list(k, {"sip", "norefine", "wcs", "log", "progress", "rebuild", "nocache",
                       "cacheinfo", "h", "-help"});
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

  const double tolerance = has("t") ? std::atof(val("t").c_str()) : 0.007;
  const int max_stars = has("s") ? std::atoi(val("s").c_str()) : 500;
  const double min_star_size = has("m") ? std::atof(val("m").c_str()) : 1.5;
  const int forced_bin = has("z") ? std::atoi(val("z").c_str()) : 0;
  const double fov = has("fov") ? std::atof(val("fov").c_str()) : 0;
  // -tiers names a ladder outright; -maxtier just raises the ceiling of the
  // default one, which is the common case: a narrow field needs deeper rungs and
  // nothing else about the ladder changes.
  const double maxtier = has("maxtier") ? std::atof(val("maxtier").c_str()) : 0;
  std::vector<double> ladder = has("tiers") ? parse_densities(val("tiers")) : kDefaultLadder;
  if (!has("tiers"))
    for (double d : kDeepLadder)
      if (d <= maxtier) ladder.push_back(d);

  // --- locate the star database ---------------------------------------------
  const std::string dbpath = has("d") ? with_separator(val("d")) : dir_of(argv[0]);
  astap::StarDatabase db;
  const bool have_db = db.select(dbpath, has("D") ? val("D") : "auto", 1.0);

  // --- resolve the index cache ----------------------------------------------
  const bool nocache = has("nocache");
  std::string cache = has("i") ? val("i")
                               : astap::default_index_cache_path(have_db ? db.name() : "unknown",
                                                                 have_db ? db.database_type() : 0,
                                                                 tolerance);
  if (has("cacheinfo")) {
    std::cout << "star database: " << (have_db ? db.name() : "(none found in " + dbpath + ")")
              << "\nladder wanted:";
    for (double d : ladder) std::cout << " " << d;
    std::cout << "\n";
    astap::QuadIndexFile info;
    std::string err;
    double total = 0;
    for (double d : ladder) {
      const std::string path = has("i") ? cache
                                        : astap::index_tier_cache_path(
                                              have_db ? db.name() : "unknown",
                                              have_db ? db.database_type() : 0, tolerance, d);
      if (astap::read_index_file_header(path, info, &err)) {
        std::cout << "  tier " << d << ": " << info.bytes / 1e9 << " GB  " << path << "\n";
        total += info.bytes / 1e9;
      } else {
        std::cout << "  tier " << d << ": not built\n";
      }
      if (has("i")) break;
    }
    // The pre-split layout, still read when it is there.
    if (!has("i") && astap::read_index_file_header(cache, info, &err)) {
      std::cout << "  a ladder file from before the tiers were split is also present ("
                << info.densities.size() << " tiers, " << info.bytes / 1e9 << " GB): " << cache
                << "\n";
    }
    std::cout << "  total " << total << " GB\n";
    return 0;
  }

  if (images.empty()) {
    print_usage();
    return 0;
  }

  // --- get an index ----------------------------------------------------------
  //
  // A rung at a time. Each is cached in its own file, so the rungs this run
  // shares with a previous one are read back and only the new ones are built:
  // raising the ceiling costs the deep rung, not the whole ladder. -i overrides
  // that with one file holding everything, which is also the shape older caches
  // have, so an existing ladder file is still read for the rungs it covers
  // instead of being thrown away.
  std::vector<astap::QuadIndex> tiers(ladder.size());
  std::vector<bool> have(ladder.size(), false);
  auto t0 = Clock::now();
  std::string err;
  size_t from_cache = 0;

  auto take_from = [&](const std::string &path, bool quiet) {
    astap::QuadIndexFile info;
    if (!astap::read_index_file_header(path, info, &err)) return;
    if (std::fabs(info.quad_tolerance - tolerance) > 1e-12) {
      if (!quiet)
        log("Index cache " + path + " was built at tolerance " +
            astap::float_to_str(info.quad_tolerance, 4) + ", rebuilding at " +
            astap::float_to_str(tolerance, 4) + ".");
      return;
    }
    std::vector<astap::QuadIndex> in;
    if (!astap::load_index_file(path, in, 0, 0, &err)) {
      if (!quiet) log("Index cache " + path + " unusable (" + err + "), rebuilding it.");
      return;
    }
    for (astap::QuadIndex &ix : in)
      for (size_t i = 0; i < ladder.size(); i++)
        if (!have[i] && std::fabs(ix.settings().star_density - ladder[i]) < 1e-9) {
          tiers[i] = std::move(ix);
          have[i] = true;
          from_cache++;
          break;
        }
  };

  if (!nocache && !has("rebuild")) {
    if (has("i")) {
      take_from(cache, false);
    } else {
      for (size_t i = 0; i < ladder.size(); i++)
        take_from(astap::index_tier_cache_path(have_db ? db.name() : "unknown",
                                               have_db ? db.database_type() : 0, tolerance,
                                               ladder[i]),
                  true);
      // Whatever is still missing may be in a ladder file written before the
      // rungs were split apart.
      if (from_cache < ladder.size()) take_from(cache, true);
    }
  }

  std::vector<double> missing;
  for (size_t i = 0; i < ladder.size(); i++)
    if (!have[i]) missing.push_back(ladder[i]);

  if (missing.empty()) {
    log("Index: " + std::to_string(tiers.size()) + " tiers read from cache in " +
        astap::float_to_str(secs(t0, Clock::now()), 2) + " sec.");
  } else {
    if (!have_db) {
      log("No star database found in " + dbpath);
      for (const std::string &f : images) {
        astap::Header head;
        astap::write_ini(change_file_ext(has("o") ? val("o") : f, ".ini"), false, head, cmdline,
                         astap::kErrNoStarDatabase, "");
      }
      return astap::kErrNoStarDatabase;
    }
    astap::QuadIndexSettings qs;
    qs.quad_tolerance = tolerance;
    if (progress)
      log("Building " + std::to_string(missing.size()) + " index tier(s) from " + db.name() +
          "; each one happens once.");
    // One pass over the database covers every missing rung: a tile is read once
    // at the deepest of them and the shallower ones take a prefix.
    std::vector<astap::QuadIndex> built;
    if (!astap::build_tiers(db, qs, missing, built)) {
      log("Could not read the star database in " + dbpath);
      return astap::kErrStarDatabaseRead;
    }
    for (astap::QuadIndex &ix : built) {
      const double d = ix.settings().star_density;
      if (!nocache && !has("i")) {
        const std::string path = astap::index_tier_cache_path(db.name(), db.database_type(),
                                                              tolerance, d);
        std::vector<astap::QuadIndex> one;
        one.push_back(ix);
        if (!astap::ensure_parent_directory(path) || !astap::save_index_file(path, one))
          log("Warning: could not write the index cache to " + path);
      }
      for (size_t i = 0; i < ladder.size(); i++)
        if (!have[i] && std::fabs(d - ladder[i]) < 1e-9) {
          tiers[i] = std::move(ix);
          have[i] = true;
          break;
        }
    }
    if (!nocache && has("i")) {
      if (astap::ensure_parent_directory(cache) && astap::save_index_file(cache, tiers))
        log("Cached to " + cache);
      else
        log("Warning: could not write the index cache to " + cache);
    }
    size_t nq = 0;
    for (const astap::QuadIndex &ix : tiers) nq += ix.size();
    log("Index: " + std::to_string(tiers.size()) + " tiers (" + std::to_string(from_cache) +
        " cached, " + std::to_string(missing.size()) + " built), " + std::to_string(nq) +
        " quads, ready in " + astap::float_to_str(secs(t0, Clock::now()), 2) + " sec.");
  }

  // A rung that could not be built at all would otherwise sit in the sweep as an
  // empty index; drop it so the tier count reported is the tier count used.
  for (size_t i = tiers.size(); i-- > 0;)
    if (!have[i]) tiers.erase(tiers.begin() + static_cast<long>(i));
  if (tiers.empty()) {
    log("No usable index tiers.");
    return astap::kErrStarDatabaseRead;
  }

  // --- solve ------------------------------------------------------------------
  int worst = 0;
  for (const std::string &filename : images) {
    astap::Header head;
    astap::ImageArray img;
    const astap::ImageLoadResult r = astap::load_image(filename, head, img);
    const std::string out_base = has("o") ? val("o") : filename;
    if (!r.ok) {
      log(r.error);
      astap::write_ini(change_file_ext(out_base, ".ini"), false, head, cmdline,
                       astap::kErrImageRead, "");
      worst = std::max(worst, static_cast<int>(astap::kErrImageRead));
      continue;
    }
    if (!r.warning.empty()) log(r.warning);

    astap::ImageArray small;
    const int bin = bin_mono(img, head.width, head.height, forced_bin, small);

    // A star size in arcseconds only means something once the plate scale is
    // known, which a blind solve does not have. With -fov it does.
    double hfd_min = 0.8;
    if (fov > 0) {
      const double arcsec_per_px = fov * 3600 / head.height;
      hfd_min = std::max(0.8, min_star_size / (bin * arcsec_per_px));
    }

    astap::Header bhead = head;
    astap::Histogram hist;
    astap::get_background(0, small, bhead, true, true, max_stars, hist);
    astap::RowList stars;
    double mean_hfd = 0;
    astap::find_stars(small, bhead, hfd_min, max_stars, stars, mean_hfd,
                      progress ? log : astap::LogFn());
    if (stars.count() < 4) {
      log("Not enough stars detected in " + filename);
      astap::write_ini(change_file_ext(out_base, ".ini"), false, head, cmdline,
                       astap::kErrNotEnoughStars, "");
      worst = std::max(worst, static_cast<int>(astap::kErrNotEnoughStars));
      continue;
    }

    if (progress)
      log("Image: " + std::to_string(stars.count()) + " stars, binning " + std::to_string(bin));

    // The tier sweep is ordered by the image's star density when the field size
    // is known; without -fov every tier is tried, cheapest first.
    double density_hint = 0;
    if (fov > 0) {
      const double w_deg = fov * head.width / head.height;
      density_hint = stars.count() / std::max(1e-9, fov * w_deg);
    }

    const auto s0 = Clock::now();
    astap::IndexSolveResult res = astap::solve_stars_with_tiers(
        tiers, stars, small.width(), small.height(), {}, density_hint);
    const double elapsed = secs(s0, Clock::now());

    if (!res.solved) {
      log("No solution for " + filename + " (" + res.reason + ", " +
          std::to_string(res.nr_matches) + " candidate quads)");
      astap::write_ini(change_file_ext(out_base, ".ini"), false, head, cmdline, astap::kErrNone,
                       "");
      worst = std::max(worst, 1);
      continue;
    }

    // Second pass: with the position known, read the database once at it and
    // redo the match there. This is what lifts the accuracy and produces enough
    // matched quads for a distortion fit.
    astap::SipCoefficients sip;
    if (has("norefine") && has("sip"))
      log("Note: -sip needs the second pass, which -norefine disables. No SIP written.");
    if (!has("norefine")) {
      const auto p0 = Clock::now();
      astap::IndexRefineResult ref = astap::refine_with_database(
          db, stars, small.width(), small.height(), res, {}, has("sip"));
      const double refine_secs = secs(p0, Clock::now());
      if (ref.ok && progress)
        log("Second pass: " + std::to_string(ref.nr_quads) + " quads matched against " +
            std::to_string(ref.nr_candidates) + " database quads in " +
            astap::float_to_str(refine_secs * 1000, 1) + " ms" +
            (ref.residual_before >= 0
                 ? ", residual " + astap::float_to_str(ref.residual_before, 3) + "\" -> " +
                       astap::float_to_str(ref.residual_after, 3) + "\""
                 : "") +
            (ref.kept ? "." : ", discarded: " + ref.reason));
      if (!ref.ok && progress) log("Second pass skipped: " + ref.reason);
      if (ref.sip_valid) {
        sip = ref.sip;
        scale_sip_to_original(sip, bin);
      } else if (has("sip")) {
        log("No SIP coefficients: " + ref.reason);
      }
    }

    scale_solution_to_original(res, bin, head);
    log("Solution found: " + astap::prepare_ra(head.ra0, ": ") + " " +
        astap::prepare_dec(head.dec0, "d "));
    log("Solved in " + astap::float_to_str(elapsed, 3) + " sec. Scale " +
        astap::float_to_str(std::fabs(head.cdelt2) * 3600, 4) + "\"/px, rotation " +
        astap::float_to_str(head.crota2, 2) + "d, " + std::to_string(res.nr_inliers) +
        " quads at depth tier " + astap::float_to_str(res.tier_density, 1) + " stars/deg^2 (" +
        std::to_string(res.tiers_tried) + " tiers tried" +
        (res.many_quads_pass ? ", larger quad groups" : "") +
        (res.stars_used < res.stars_detected
             ? ", brightest " + std::to_string(res.stars_used) + " of " +
                   std::to_string(res.stars_detected) + " stars"
             : "") +
        (res.refined ? ", refined" : "") + (sip.valid ? ", SIP" : "") + ").");

    astap::write_ini(change_file_ext(out_base, ".ini"), true, head, cmdline, astap::kErrNone, "");

    if (has("wcs")) {
      const std::string comment = "Solved in " + astap::float_to_str(elapsed, 3) +
                                  " sec by the index method.";
      astap::update_solution_cards(head.cards, head, sip, true, comment);
      astap::remove_key(head.cards, "NAXIS1  =");
      astap::remove_key(head.cards, "NAXIS2  =");
      astap::update_integer(head.cards, "NAXIS   =",
                            " / Minimal header                                 ", 0);
      astap::update_integer(head.cards, "BITPIX  =",
                            " /                                                ", 8);
      astap::write_fits_header_file(change_file_ext(out_base, ".wcs"), head.cards);
    }
  }

  if (want_log) {
    std::ofstream lf(change_file_ext(has("o") ? val("o") : images.front(), ".log"));
    lf << cmdline << "\n";
    for (const std::string &l : log_lines) lf << l << "\n";
  }
  return worst;
}
