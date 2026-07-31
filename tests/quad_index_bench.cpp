// Exercises the quad index against a real image: build the index, look up the
// image quads, and check that the matches concentrate on the true position.
// This is the first end-to-end evidence that design B can work; it is not yet
// a solver. See docs/index_solver.md.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <map>
#include <string>
#include <vector>

#include "astap/astro_math.h"
#include "astap/fits.h"
#include "astap/index_solver.h"
#include "astap/quad_index.h"
#include "astap/quads.h"
#include "astap/star_detection.h"

using namespace astap;
using Clock = std::chrono::steady_clock;

static double secs(Clock::time_point a, Clock::time_point b) {
  return std::chrono::duration<double>(b - a).count();
}

int main(int argc, char** argv) {
  if (argc < 3) {
    std::printf("usage: quad_index_bench <image.fits> <database dir> "
                "[density] [true_ra_deg true_dec_deg]\n");
    return 2;
  }
  const std::string image = argv[1];
  std::string dbpath = argv[2];
  if (dbpath.back() != '/') dbpath += '/';
  const double density = argc > 3 ? std::atof(argv[3]) : 300.0;
  const bool have_truth = argc > 5;
  const double true_ra = have_truth ? std::atof(argv[4]) * kPi / 180 : 0;
  const double true_dec = have_truth ? std::atof(argv[5]) * kPi / 180 : 0;

  StarDatabase db;
  if (!db.select(dbpath, "auto", 1.0)) {
    std::printf("no star database in %s\n", dbpath.c_str());
    return 1;
  }
  std::printf("database %s (type %d)\n", db.name().c_str(), db.database_type());

  // --- build the index -----------------------------------------------------
  QuadIndexSettings s;
  s.star_density = density;
  QuadIndex index;
  auto t0 = Clock::now();
  if (!index.build(db, s)) {
    std::printf("index build failed\n");
    return 1;
  }
  auto t1 = Clock::now();
  std::printf("index: %zu quads, %.0f MB, built in %.2f s (density %.0f stars/deg^2)\n",
              index.size(), index.bytes() / 1e6, secs(t0, t1), density);

  // --- image side ----------------------------------------------------------
  Header head;
  ImageArray img;
  if (!load_fits(image, head, img).ok) {
    std::printf("cannot read %s\n", image.c_str());
    return 1;
  }
  // Mono, binned so the working image lands near 1500 px on its long side.
  // Binning a small image hard, as an earlier fixed factor of 3 did, throws
  // away the stars the solve depends on.
  int bin = 1;
  while (std::max(head.width, head.height) / (bin + 1) >= 1200) bin++;
  ImageArray small(1, head.height / bin, head.width / bin);
  for (int y = 0; y < small.height(); y++)
    for (int x = 0; x < small.width(); x++) {
      double v = 0;
      for (int c = 0; c < img.colours(); c++)
        for (int i = 0; i < bin; i++)
          for (int j = 0; j < bin; j++) v += img.at(c, y * bin + i, x * bin + j);
      small.at(0, y, x) = static_cast<float>(v / (img.colours() * bin * bin));
    }

  Histogram hist;
  get_background(0, small, head, true, true, 500, hist);
  RowList stars;
  double mean_hfd = 0;
  find_stars(small, head, 0.8, 500, stars, mean_hfd);

  RowList iquads;
  find_quads(static_cast<int>(stars.count()), stars, iquads);
  std::printf("image: %zu stars, %zu quads\n", stars.count(), iquads.count());

  // --- query ---------------------------------------------------------------
  t0 = Clock::now();
  std::vector<uint32_t> hits;
  struct Pair { uint32_t iq, dq; double scale; };
  std::vector<Pair> pairs;
  for (size_t q = 0; q < iquads.count(); q++) {
    float r[5];
    for (int k = 0; k < 5; k++) r[k] = static_cast<float>(iquads(k + 1, q));
    hits.clear();
    index.query(r, hits);
    for (uint32_t h : hits)
      pairs.push_back({static_cast<uint32_t>(q), h, index.d1(h) / iquads(0, q)});
  }
  t1 = Clock::now();
  std::printf("query: %zu candidate pairs in %.3f s (%.1f per image quad)\n", pairs.size(),
              secs(t0, t1), pairs.empty() ? 0.0 : double(pairs.size()) / iquads.count());
  if (pairs.empty()) return 1;

  // --- solve ---------------------------------------------------------------
  IndexSolveSettings ss;
  t0 = Clock::now();
  IndexSolveResult res = solve_with_index(index, iquads, small.width(), small.height(), ss);
  t1 = Clock::now();
  if (!res.solved) {
    std::printf("NOT SOLVED (%s), %d candidate pairs\n", res.reason.c_str(), res.nr_matches);
    return 1;
  }
  std::printf("solved in %.3f s: %s %s, %.4f arcsec/px, %.3f deg rotation, %d quads\n",
              secs(t0, t1), prepare_ra(res.ra0, ": ").c_str(),
              prepare_dec(res.dec0, "d ").c_str(), res.scale_arcsec_px, res.crota2,
              res.nr_references);

  if (have_truth) {
    double sep;
    ang_sep(res.ra0, res.dec0, true_ra, true_dec, sep);
    std::printf("distance from the known solution: %.2f arcsec\n", sep * 3600 * 180 / kPi);
    if (sep * 180 / kPi < 0.5) {
      std::printf("\nPASS\n");
      return 0;
    }
    std::printf("\nFAIL: solved, but not at the true position\n");
    return 1;
  }
  return 0;
}
