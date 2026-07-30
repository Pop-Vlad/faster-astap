// Self contained checks for the ported routines. No external test framework.

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <string>
#include <vector>

#include "astap/astro_math.h"
#include "astap/matching.h"
#include "astap/quads.h"
#include "astap/star_database.h"
#include "astap/star_detection.h"

using namespace astap;

namespace {
  int failures = 0;

  void check(bool cond, const std::string &what) {
    if (!cond) {
      std::printf("FAIL: %s\n", what.c_str());
      failures++;
    } else {
      std::printf("ok  : %s\n", what.c_str());
    }
  }

  void check_near(double got, double want, double tol, const std::string &what) {
    if (!(std::fabs(got - want) <= tol)) {
      std::printf("FAIL: %s (got %.10g, want %.10g, tol %g)\n", what.c_str(), got, want, tol);
      failures++;
    } else {
      std::printf("ok  : %s\n", what.c_str());
    }
  }

  // ---------------------------------------------------------------------------

  void test_math() {
    check_near(fnmodulo(-0.5, 2 * kPi), 2 * kPi - 0.5, 1e-12, "fnmodulo wraps negative values");
    check_near(fnmodulo(7.0, 2 * kPi), 7.0 - 2 * kPi, 1e-12, "fnmodulo wraps large values");

    std::vector<double> v = {5, 1, 3};
    check_near(smedian(v, 3), 3, 1e-12, "smedian of three values");
    std::vector<double> v2 = {1, 2, 3, 4};
    check_near(smedian(v2, 4), 2.5, 1e-12, "smedian of four values");
    // For an odd count above three the original averages the three middle values.
    std::vector<double> v3 = {1, 2, 3, 4, 100};
    check_near(smedian(v3, 5), (2 + 3 + 4) / 3.0, 1e-12, "smedian averages three middle values");

    // equatorial_standard and standard_equatorial must be each other's inverse.
    const double ra0 = 1.234, dec0 = 0.567;
    const double ra = ra0 + 0.01, dec = dec0 - 0.008;
    double x, y, ra_back, dec_back;
    equatorial_standard(ra0, dec0, ra, dec, 1.0, x, y);
    standard_equatorial(ra0, dec0, x, y, 1.0, ra_back, dec_back);
    check_near(ra_back, ra, 1e-12, "standard_equatorial inverts equatorial_standard (ra)");
    check_near(dec_back, dec, 1e-12, "standard_equatorial inverts equatorial_standard (dec)");

    // A star due north of the reference has a position angle of zero.
    check_near(position_angle(0.5, 0.31, 0.5, 0.3), 0.0, 1e-9, "position angle north is zero");
    // A star at the same declination is very close to due east. It is not exactly
    // 90 degrees: a line of constant declination is not a great circle, so the
    // angle falls slightly short of a right angle.
    check_near(position_angle(0.51, 0.3, 0.5, 0.3), kPi / 2, 5e-3, "position angle east is ~+90d");

    double ra_txt;
    check(ra_text_to_radians("12 30 00", ra_txt), "ra_text_to_radians parses 12 30 00");
    check_near(ra_txt, 12.5 * kPi / 12, 1e-12, "ra_text_to_radians value");
    double dec_txt;
    check(dec_text_to_radians("-05 30 00", dec_txt), "dec_text_to_radians parses -05 30 00");
    check_near(dec_txt, -5.5 * kPi / 180, 1e-12, "dec_text_to_radians value");
  }

  // ---------------------------------------------------------------------------

  void test_lsq_fit() {
    // Six equations for X := 2x + 0.5y + 10.
    RowList a(3, 6);
    std::vector<double> b(6);
    const double px[6] = {0, 100, 200, 50, 150, 250};
    const double py[6] = {0, 30, 60, 200, 10, 90};
    for (size_t i = 0; i < 6; i++) {
      a(0, i) = px[i];
      a(1, i) = py[i];
      a(2, i) = 1;
      b[i] = 2 * px[i] + 0.5 * py[i] + 10;
    }
    SolutionVector s;
    check(lsq_fit(a, b, s), "lsq_fit succeeds");
    check_near(s[0], 2.0, 1e-9, "lsq_fit recovers a");
    check_near(s[1], 0.5, 1e-9, "lsq_fit recovers b");
    check_near(s[2], 10.0, 1e-9, "lsq_fit recovers c");

    // The input matrix must not be modified.
    check_near(a(0, 1), 100.0, 1e-12, "lsq_fit leaves the A matrix untouched");
  }

  // ---------------------------------------------------------------------------

  // Builds a random star field, transforms it with a known affine transformation
  // and checks that the quad matching recovers that transformation.
  void test_quad_matching() {
    std::mt19937 rng(12345);
    std::uniform_real_distribution<double> ux(0, 2000), uy(0, 1500);

    const size_t n = 300;
    RowList db(3, n); // the "database" list
    for (size_t i = 0; i < n; i++) {
      db(0, i) = ux(rng);
      db(1, i) = uy(rng);
      db(2, i) = 100;
    }

    // The image is the same field rotated by 20 degrees, scaled by 0.5 and
    // shifted. The solver solves for image -> database, so build the image by
    // applying the inverse of what we expect to recover.
    const double rot = 20 * kPi / 180;
    const double scale = 0.5;
    RowList image(3, n);
    for (size_t i = 0; i < n; i++) {
      const double x = db(0, i), y = db(1, i);
      image(0, i) = (std::cos(rot) * x - std::sin(rot) * y) * scale + 30;
      image(1, i) = (std::sin(rot) * x + std::cos(rot) * y) * scale - 15;
      image(2, i) = 100;
    }

    // find_quads sorts the star list in X in place when there are 150 stars or
    // more, so keep copies for the verification below.
    RowList image_orig = image;
    RowList db_orig = db;

    MatchState st;
    find_quads(static_cast<int>(n), image, st.quad_star_distances2);
    find_quads(static_cast<int>(n), db, st.quad_star_distances1);

    check(st.quad_star_distances2.count() > 100, "image quads found");
    check(st.quad_star_distances1.count() > 100, "database quads found");
    check(image(0, 0) <= image(0, n - 1), "find_quads sorted the star list in X in place");

    check(find_offset_and_rotation(st, 3, 0.007), "find_offset_and_rotation solves the field");
    check(st.nr_references >= 3, "at least three matching quads survive the outlier filter");

    // The recovered solution must map an image position back onto the database
    // position for every star.
    double worst = 0;
    for (size_t i = 0; i < n; i++) {
      const double xi = image_orig(0, i), yi = image_orig(1, i);
      const double X = st.solution_vector_x[0] * xi + st.solution_vector_x[1] * yi +
                       st.solution_vector_x[2];
      const double Y = st.solution_vector_y[0] * xi + st.solution_vector_y[1] * yi +
                       st.solution_vector_y[2];
      worst = std::max(worst, std::max(std::fabs(X - db_orig(0, i)), std::fabs(Y - db_orig(1, i))));
    }
    check_near(worst, 0.0, 1e-6, "solution maps every image star onto its database star");

    // The solution is the inverse transformation, so its scale is 1/scale and,
    // written as [[cos, sin],[-sin, cos]]/scale, its rotation is +rot.
    const double recovered_scale =
        std::sqrt(sqr(st.solution_vector_x[0]) + sqr(st.solution_vector_x[1]));
    check_near(recovered_scale, 1 / scale, 1e-6, "solution recovers the scale");
    const double recovered_rot = std::atan2(st.solution_vector_x[1], st.solution_vector_x[0]);
    check_near(recovered_rot, rot, 1e-6, "solution recovers the rotation");
  }

  // A field with few stars takes the find_many_quads path (groups of 5, 6 or 7).
  void test_many_quads_path() {
    std::mt19937 rng(999);
    std::uniform_real_distribution<double> ux(0, 500), uy(0, 400);

    const size_t n = 25; // below 30, so groups of six stars are used
    RowList db(3, n);
    for (size_t i = 0; i < n; i++) {
      db(0, i) = ux(rng);
      db(1, i) = uy(rng);
      db(2, i) = 100;
    }
    RowList image(3, n);
    for (size_t i = 0; i < n; i++) {
      image(0, i) = db(0, i) * 2 + 7;
      image(1, i) = db(1, i) * 2 - 3;
      image(2, i) = 100;
    }

    MatchState st;
    find_quads(static_cast<int>(n), image, st.quad_star_distances2);
    find_quads(static_cast<int>(n), db, st.quad_star_distances1);
    // C(6,4) = 15 quads per star, minus the duplicates.
    check(st.quad_star_distances2.count() > n, "find_many_quads produces several quads per star");
    check(find_offset_and_rotation(st, 3, 0.007), "few-star field solves through find_many_quads");
    check_near(st.solution_vector_x[0], 0.5, 1e-6, "few-star solution recovers the scale");
    check_near(st.solution_vector_x[2], -3.5, 1e-6, "few-star solution recovers the x offset");
  }

  // ---------------------------------------------------------------------------

  // Renders gaussian stars into an image and checks that find_stars locates them.
  void test_find_stars() {
    const int w = 400, h = 300;
    ImageArray img(1, h, w);
    const float background = 1000.0f;

    std::mt19937 rng(7);
    std::normal_distribution<float> noise(0.0f, 8.0f);
    for (int y = 0; y < h; y++)
      for (int x = 0; x < w; x++) img.at(0, y, x) = background + noise(rng);

    struct P {
      double x, y, amp;
    };
    const std::vector<P> stars = {
      {50.4, 60.2, 6000}, {120.0, 200.0, 9000}, {300.7, 100.3, 12000},
      {200.0, 150.0, 4000}, {350.0, 250.0, 15000}
    };
    const double sigma = 1.6;
    for (const P &s: stars)
      for (int y = static_cast<int>(s.y) - 10; y <= static_cast<int>(s.y) + 10; y++)
        for (int x = static_cast<int>(s.x) - 10; x <= static_cast<int>(s.x) + 10; x++) {
          if (x < 0 || y < 0 || x >= w || y >= h) continue;
          const double r2 = sqr(x - s.x) + sqr(y - s.y);
          img.at(0, y, x) += static_cast<float>(s.amp * std::exp(-r2 / (2 * sigma * sigma)));
        }

    Header head;
    head.bitpix = 16;
    head.width = w;
    head.height = h;
    Histogram hist;
    get_background(0, img, head, true, true, 500, hist);

    check_near(head.backgr, background, 30, "get_background finds the background level");
    check(head.noise_level > 2 && head.noise_level < 20, "get_background finds a sane noise level");

    RowList found;
    double mean_hfd = 0;
    find_stars(img, head, 0.8, 500, found, mean_hfd);

    check(found.count() >= stars.size(), "find_stars finds all planted stars");

    // Every planted star must be matched within half a pixel.
    for (const P &s: stars) {
      double best = 1e9;
      for (size_t i = 0; i < found.count(); i++)
        best = std::min(best, std::sqrt(sqr(found(0, i) - s.x) + sqr(found(1, i) - s.y)));
      check_near(best, 0.0, 0.5, "star at " + float_to_str(s.x, 1) + "," + float_to_str(s.y, 1) +
                                 " located within half a pixel");
    }
    check(mean_hfd > 1.0 && mean_hfd < 6.0, "mean HFD is plausible for sigma 1.6 stars");
  }

  // A single star measured directly by hfd().
  void test_hfd() {
    const int w = 100, h = 100;
    ImageArray img(1, h, w);
    for (int y = 0; y < h; y++)
      for (int x = 0; x < w; x++) img.at(0, y, x) = 500.0f;

    const double cx = 50.0, cy = 50.0, sigma = 2.0;
    for (int y = 30; y < 70; y++)
      for (int x = 30; x < 70; x++)
        img.at(0, y, x) += static_cast<float>(
          20000 * std::exp(-(sqr(x - cx) + sqr(y - cy)) / (2 * sigma * sigma)));

    double hfd1, fwhm, snr, flux, xc, yc;
    hfd(img, 50, 50, 14, hfd1, fwhm, snr, flux, xc, yc);
    check_near(xc, cx, 0.05, "hfd centre of gravity x");
    check_near(yc, cy, 0.05, "hfd centre of gravity y");
    // For a gaussian the HFD is close to 2*sqrt(2*ln2)*sigma * 1.0 ~= FWHM.
    check_near(hfd1, 2.3548 * sigma, 0.7, "hfd of a gaussian approximates the FWHM");
    check(snr > 50, "hfd reports a high SNR for a bright star");
  }

  // ---------------------------------------------------------------------------

  // The area numbering and file naming of the star database.
  void test_database_areas() {
    StarDatabase db290;
    StarDatabase db1476;
    // select() needs files on disk, so exercise find_areas through a helper that
    // does not depend on them: the default type is 1476, and a 290 instance is
    // created by selecting a non-existing database (which leaves the default).
    // Instead the public surface is checked for the default .1476 layout.

    int a1, a2, a3, a4;
    double f1, f2, f3, f4;

    // A small field in the middle of a tile must resolve to exactly one area.
    db1476.find_areas(1.0, 0.3, 0.2 * kPi / 180, a1, a2, a3, a4, f1, f2, f3, f4);
    check(a1 >= 1 && a1 <= 1476, "1476 area number is in range");
    check(a2 == 0 && a3 == 0 && a4 == 0, "a small field needs only one 1476 area");
    check_near(f1, 1.0, 1e-9, "the single area covers the whole field");

    // The celestial poles.
    db1476.find_areas(0.0, kPi / 2 - 0.001, 0.2 * kPi / 180, a1, a2, a3, a4, f1, f2, f3, f4);
    check(a1 == 1476, "the north pole maps to area 1476");
    db1476.find_areas(0.0, -kPi / 2 + 0.001, 0.2 * kPi / 180, a1, a2, a3, a4, f1, f2, f3, f4);
    check(a1 == 1, "the south pole maps to area 1");

    // A field straddling a tile boundary needs more than one area and the
    // fractions must add up to one.
    db1476.find_areas(0.0, 0.0, 4.0 * kPi / 180, a1, a2, a3, a4, f1, f2, f3, f4);
    check(a2 != 0 || a3 != 0 || a4 != 0, "a field on a boundary needs several areas");
    check_near(f1 + f2 + f3 + f4, 1.0, 1e-6, "the area fractions add up to one");

    (void) db290;
  }
} // namespace

int main() {
  std::printf("-- math --\n");
  test_math();
  std::printf("-- lsq_fit --\n");
  test_lsq_fit();
  std::printf("-- hfd --\n");
  test_hfd();
  std::printf("-- find_stars --\n");
  test_find_stars();
  std::printf("-- quad matching --\n");
  test_quad_matching();
  std::printf("-- find_many_quads --\n");
  test_many_quads_path();
  std::printf("-- star database areas --\n");
  test_database_areas();

  if (failures) {
    std::printf("\n%d check(s) FAILED\n", failures);
    return 1;
  }
  std::printf("\nall checks passed\n");
  return 0;
}
