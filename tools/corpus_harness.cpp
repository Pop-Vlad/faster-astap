// Runs both solvers over a corpus and reports solve rate and position agreement.
//
// This is the capability gate. The index solver is allowed to be no better than
// the port at finding hard fields, but it must not be worse: every image the
// port solves has to solve through the index too. One test image cannot tell a
// 300x speedup from one that quietly stops solving five percent of images, so no
// claim about the index solver should be made without this table.
//
// Ground truth comes from each file's own WCS (SkyView writes CRVAL1/2 and
// CDELT1/2 when it makes the cutout), so a solve is checked against the sky, not
// against the other solver.
//
//   corpus_harness <corpus dir> <database dir> [options]
//     --density 100,200,300   index depth tiers to try, in order
//     --hint                  give both solvers the true position and scale
//     --csv <file>            also write the table as CSV
//     --limit N               stop after N images
//     --filter <substring>    only images whose name contains this
//     --tol <x>               quad matching tolerance, default 0.007

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <dirent.h>
#include <fstream>
#include <string>
#include <vector>

#include "astap/astro_math.h"
#include "astap/fits.h"
#include "astap/index_solver.h"
#include "astap/quad_index.h"
#include "astap/quads.h"
#include "astap/solver.h"
#include "astap/star_detection.h"

using namespace astap;
using Clock = std::chrono::steady_clock;

namespace {
  double secs(Clock::time_point a, Clock::time_point b) {
    return std::chrono::duration<double>(b - a).count();
  }

  double arcsec_between(double ra1, double dec1, double ra2, double dec2) {
    double sep;
    ang_sep(ra1, dec1, ra2, dec2, sep);
    return sep * 3600 * 180 / kPi;
  }

  // A solve counts only when it lands on the right field. "Right" has to scale
  // with the image: 60 arcsec is a coarse error on a half degree frame and two
  // pixels on a ten degree one, so the threshold is the larger of a fixed floor
  // and two percent of the field.
  double correct_within_arcsec(double fov_deg) { return std::max(60.0, 0.02 * fov_deg * 3600); }

  // The index solver works on a mono, binned copy, the same one quad_index_bench
  // uses. Binning a small image hard throws away the stars the solve depends on,
  // so the factor is chosen to land the long side near 1200 px and no lower.
  int binning_for(int width, int height) {
    int bin = 1;
    while (std::max(width, height) / (bin + 1) >= 1200) bin++;
    return bin;
  }

  struct ImageCase {
    std::string file;
    int width = 0, height = 0;         // original
    int bwidth = 0, bheight = 0;       // binned, what the index solver sees
    double true_ra = 0, true_dec = 0;  // radians, from the file's own WCS
    double true_scale = 0;             // arcsec per pixel, original pixels
    double fov_deg = 0;
    size_t nstars = 0;
    double density = 0;  // detected stars per square degree
    RowList stars;

    // Port
    bool port_solved = false;
    bool port_good = false;
    double port_error = 0;  // arcsec from truth
    double port_secs = 0;
    double port_scale = 0;

    // Index solver, sweeping the whole tier ladder
    bool idx_solved = false;
    double idx_error = 0;
    double idx_secs = 0;
    int idx_inliers = 0;
    double idx_tier = 0;
    int idx_tiers_tried = 0;
    bool idx_many = false;
    std::string idx_reason;
  };

  std::vector<std::string> list_fits(const std::string &dir, const std::string &filter) {
    std::vector<std::string> out;
    DIR *d = opendir(dir.c_str());
    if (!d) return out;
    while (dirent *e = readdir(d)) {
      const std::string n = e->d_name;
      if (n.size() < 6 || n.compare(n.size() - 5, 5, ".fits") != 0) continue;
      if (!filter.empty() && n.find(filter) == std::string::npos) continue;
      out.push_back(dir + "/" + n);
    }
    closedir(d);
    std::sort(out.begin(), out.end());
    return out;
  }

  std::string basename_of(const std::string &p) {
    const size_t s = p.find_last_of('/');
    return s == std::string::npos ? p : p.substr(s + 1);
  }
} // namespace

int main(int argc, char **argv) {
  if (argc < 3) {
    std::printf(
        "usage: corpus_harness <corpus dir> <database dir> [--density 100,200,300] [--hint]\n"
        "                      [--csv out.csv] [--limit N] [--filter substring] [--tol x]\n");
    return 2;
  }
  const std::string corpus = argv[1];
  std::string dbpath = argv[2];
  if (dbpath.back() != '/') dbpath += '/';

  std::vector<double> densities;
  bool hint = false;
  std::string csv, filter;
  int limit = 0;
  double tol = 0.007;
  for (int i = 3; i < argc; i++) {
    const std::string a = argv[i];
    if (a == "--hint") hint = true;
    else if (a == "--csv" && i + 1 < argc) csv = argv[++i];
    else if (a == "--limit" && i + 1 < argc) limit = std::atoi(argv[++i]);
    else if (a == "--filter" && i + 1 < argc) filter = argv[++i];
    else if (a == "--tol" && i + 1 < argc) tol = std::atof(argv[++i]);
    else if (a == "--density" && i + 1 < argc) {
      std::string v = argv[++i];
      size_t p = 0;
      while (p <= v.size()) {
        const size_t c = v.find(',', p);
        const std::string tok = v.substr(p, c == std::string::npos ? std::string::npos : c - p);
        if (!tok.empty()) densities.push_back(std::atof(tok.c_str()));
        if (c == std::string::npos) break;
        p = c + 1;
      }
    }
  }
  // A geometric ladder. One tier reaches about a factor of two in image density
  // either side of itself, so a ratio of 2.5 between rungs covers the range
  // continuously, and the corpus spans 1 to 1300 detected stars/deg^2.
  if (densities.empty()) densities = {0.5, 1, 2, 4, 8, 16, 32, 60, 125, 250, 500, 900};

  std::vector<std::string> files = list_fits(corpus, filter);
  if (limit > 0 && static_cast<int>(files.size()) > limit) files.resize(limit);
  if (files.empty()) {
    std::printf("no .fits files in %s\n", corpus.c_str());
    return 1;
  }
  std::printf("%zu images, densities:", files.size());
  for (double d : densities) std::printf(" %.0f", d);
  std::printf(", %s\n\n", hint ? "hinted" : "blind");

  // --- pass 1: image stages and the port -------------------------------------
  std::vector<ImageCase> cases;
  for (const std::string &f : files) {
    ImageCase c;
    c.file = f;

    Header head;
    ImageArray img;
    FitsLoadResult r = load_fits(f, head, img);
    if (!r.ok) {
      std::printf("  %-42s cannot read: %s\n", basename_of(f).c_str(), r.error.c_str());
      continue;
    }
    c.width = head.width;
    c.height = head.height;
    c.true_ra = head.ra0;
    c.true_dec = head.dec0;
    c.true_scale = std::fabs(head.cdelt2) * 3600;
    c.fov_deg = c.true_scale * std::max(head.width, head.height) / 3600;
    if (c.true_scale <= 0) {
      std::printf("  %-42s no WCS in the file, skipped\n", basename_of(f).c_str());
      continue;
    }

    // Stars and quads for the index solver, on the binned mono copy.
    const int bin = binning_for(head.width, head.height);
    ImageArray small(1, head.height / bin, head.width / bin);
    for (int y = 0; y < small.height(); y++)
      for (int x = 0; x < small.width(); x++) {
        double v = 0;
        for (int col = 0; col < img.colours(); col++)
          for (int i = 0; i < bin; i++)
            for (int j = 0; j < bin; j++) v += img.at(col, y * bin + i, x * bin + j);
        small.at(0, y, x) = static_cast<float>(v / (img.colours() * bin * bin));
      }
    c.bwidth = small.width();
    c.bheight = small.height();

    Header bhead = head;
    Histogram hist;
    get_background(0, small, bhead, true, true, 500, hist);
    RowList stars;
    double mean_hfd = 0;
    find_stars(small, bhead, 0.8, 500, stars, mean_hfd);
    c.nstars = stars.count();
    c.density = c.nstars / std::max(1e-9, static_cast<double>(c.fov_deg) * c.fov_deg *
                                              (std::min(c.width, c.height) /
                                               static_cast<double>(std::max(c.width, c.height))));
    c.stars = stars;

    // The port, on the original image.
    SolverSettings ss;
    ss.database_path = dbpath;
    ss.radius_search = hint ? 5 : 180;
    if (hint) {
      ss.fov_specified = true;
      ss.search_fov = c.true_scale * head.height / 3600;
    }
    Header phead = head;
    if (!hint) {
      phead.ra0 = 0;  // blind: the solver may not start from the answer
      phead.dec0 = 0;
    }
    Solver solver(ss);
    solver.set_mount(99999, 99999);
    auto t0 = Clock::now();
    c.port_solved = solver.solve(img, phead);
    c.port_secs = secs(t0, Clock::now());
    if (c.port_solved) {
      c.port_error = arcsec_between(phead.ra0, phead.dec0, c.true_ra, c.true_dec);
      c.port_scale = std::fabs(phead.cdelt2) * 3600;
    }
    c.port_good = c.port_solved && c.port_error < correct_within_arcsec(c.fov_deg);

    std::printf("  %-42s %5.2f deg  %4zu stars (%6.0f/deg^2)  port %s\n", basename_of(f).c_str(),
                c.fov_deg, c.nstars, c.density,
                c.port_solved ? (c.port_good ? "ok" : "WRONG POSITION") : "no solution");
    cases.push_back(std::move(c));
  }
  if (cases.empty()) return 1;

  // --- pass 2: build the whole ladder once, sweep it per image ---------------
  StarDatabase db;
  if (!db.select(dbpath, "auto", 1.0)) {
    std::printf("no star database in %s\n", dbpath.c_str());
    return 1;
  }
  std::printf("\ndatabase %s (type %d)\n", db.name().c_str(), db.database_type());

  QuadIndexSettings qs;
  qs.quad_tolerance = tol;
  std::vector<QuadIndex> tiers;
  auto t0 = Clock::now();
  if (!build_tiers(db, qs, densities, tiers)) {
    std::printf("index build failed\n");
    return 1;
  }
  const double build_secs = secs(t0, Clock::now());
  size_t total_quads = 0, total_bytes = 0;
  for (const QuadIndex &ix : tiers) {
    total_quads += ix.size();
    total_bytes += ix.bytes();
    std::printf("  tier @%-6.0f %9zu quads %6.0f MB\n", ix.settings().star_density, ix.size(),
                ix.bytes() / 1e6);
  }
  std::printf("ladder: %zu tiers, %zu quads, %.1f GB, built in %.2f s\n\n", tiers.size(),
              total_quads, total_bytes / 1e9, build_secs);

  for (ImageCase &c : cases) {
    IndexSolveSettings is;
    auto s0 = Clock::now();
    // No density hint: this is the blind case, where the field size is unknown
    // and the ladder has to be swept. A hint only reorders the sweep.
    IndexSolveResult r =
        solve_stars_with_tiers(tiers, c.stars, c.bwidth, c.bheight, is, hint ? c.density : 0.0);
    c.idx_secs = secs(s0, Clock::now());
    c.idx_error = r.solved ? arcsec_between(r.ra0, r.dec0, c.true_ra, c.true_dec) : 0;
    // A solve that lands somewhere else is a failure, not a solve.
    c.idx_solved = r.solved && c.idx_error < correct_within_arcsec(c.fov_deg);
    c.idx_inliers = r.nr_inliers;
    c.idx_tier = r.tier_density;
    c.idx_tiers_tried = r.tiers_tried;
    c.idx_many = r.many_quads_pass;
    c.idx_reason = r.solved ? (c.idx_solved ? "" : "wrong position") : r.reason;
  }

  // --- report ----------------------------------------------------------------
  std::printf("%-42s %6s %6s %7s | %-15s | %-30s\n", "image", "fov", "stars", "dens",
              "port (err, time)", "index (ladder sweep)");
  std::printf("%s\n", std::string(120, '-').c_str());

  int port_ok = 0, idx_ok = 0, gate_fail = 0;
  double port_time = 0, idx_time = 0;
  for (const ImageCase &c : cases) {
    std::printf("%-42s %6.2f %6zu %7.0f | ", basename_of(c.file).c_str(), c.fov_deg, c.nstars,
                c.density);
    if (c.port_good) {
      std::printf("%5.2fpx %6.2fs |", c.port_error / c.true_scale, c.port_secs);
      port_ok++;
      port_time += c.port_secs;
      idx_time += c.idx_secs;
    } else {
      std::printf("%-14s |", c.port_solved ? "WRONG" : "no solution");
    }
    if (c.idx_solved) {
      std::printf(" %5.2fpx %6.3fs  tier %-6.0f %2d inliers%s", c.idx_error / c.true_scale,
                  c.idx_secs, c.idx_tier, c.idx_inliers, c.idx_many ? " [many]" : "");
      idx_ok++;
    } else {
      std::printf(" %-30s", c.idx_reason.c_str());
    }
    std::printf("%s\n", c.port_good && !c.idx_solved ? "   <== MISSED" : "");
    if (c.port_good && !c.idx_solved) gate_fail++;
  }

  std::printf("\nport  solved %d/%zu\n", port_ok, cases.size());
  std::printf("index solved %d/%zu\n", idx_ok, cases.size());
  if (port_ok > 0)
    std::printf("on the %d images the port solves: index takes %.2f s against the port's %.2f s "
                "(%.0fx)\n",
                port_ok, idx_time, port_time, port_time / std::max(1e-9, idx_time));
  std::printf("\nGATE: images the port solves that the index misses: %d\n", gate_fail);
  if (gate_fail == 0 && port_ok > 0) std::printf("PASS\n");

  if (!csv.empty()) {
    std::ofstream f(csv);
    f << "image,fov_deg,stars,density,port_solved,port_error_arcsec,port_secs,"
         "idx_solved,idx_error_arcsec,idx_secs,idx_tier,idx_inliers,idx_tiers_tried,"
         "idx_many_quads,idx_reason\n";
    for (const ImageCase &c : cases)
      f << basename_of(c.file) << "," << c.fov_deg << "," << c.nstars << "," << c.density << ","
        << c.port_good << "," << c.port_error << "," << c.port_secs << ","
        << c.idx_solved << "," << c.idx_error << "," << c.idx_secs << "," << c.idx_tier << ","
        << c.idx_inliers << "," << c.idx_tiers_tried << "," << c.idx_many << ",\""
        << c.idx_reason << "\"\n";
    std::printf("wrote %s\n", csv.c_str());
  }
  return gate_fail == 0 ? 0 : 1;
}
