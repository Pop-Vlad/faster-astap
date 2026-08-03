// Checks that solving from an array in memory gives the same answer as solving
// the file it came from, for both solvers.
//
// That equivalence is the point of the split: the file front ends are supposed
// to be a read and two writes wrapped around a pipeline that never touches a
// disk, and the only way that stays true is to compare the two paths.
//
// Needs a star database and one corpus image. Without either it reports itself
// skipped and passes; fetch images with tools/fetch_skyview_corpus.py.

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <filesystem>
#include <string>
#include <vector>

#include "astap/image/image_io.h"
#include "astap/solve_service.h"
#include "astap/spiral_service.h"
#include "astap/star_database.h"

using namespace astap;

static int failures = 0;

static void check(bool ok, const std::string &what) {
  if (!ok) {
    std::printf("FAIL: %s\n", what.c_str());
    failures++;
  }
}

// Both solvers report the centre in radians; a tenth of an arcsecond is far
// tighter than either solver's own scatter and far looser than the bit
// identical agreement the two paths should actually show.
static void same_position(const Header &a, const Header &b, const std::string &what) {
  const double arcsec = 180 * 3600 / 3.14159265358979323846;
  check(std::fabs(a.ra0 - b.ra0) * arcsec < 0.1, what + ": same right ascension");
  check(std::fabs(a.dec0 - b.dec0) * arcsec < 0.1, what + ": same declination");
  check(std::fabs(a.cdelt2 - b.cdelt2) * 3600 < 1e-4, what + ": same scale");
}

int main() {
  bool have_db = false;
  for (const std::string &dir : default_database_directories()) {
    StarDatabase probe;
    if (probe.select(with_separator(dir), "auto", 1.0)) {
      have_db = true;
      break;
    }
  }
  if (!have_db) {
    std::printf("SKIPPED: no star database found\n");
    return 0;
  }

  // The corpus is a capability grid, so a good part of it is images nothing can
  // solve on purpose — the 0.05 degree cutouts are down to a handful of stars
  // and are there to fix where the floor is, and at 10 degrees the port's own
  // field size sweep gives up. This test is not about which images solve, so it
  // wants one both solvers are comfortable with: around a degree, which is the
  // middle of the range for both and quick either way. The corpus tool names
  // them <field>_<survey>_<fov>deg.fits; anything else keeps its place.
  std::vector<std::pair<double, std::string>> candidates;
  std::error_code ec;
  for (const auto &e : std::filesystem::directory_iterator(ASTAP_TEST_CORPUS_DIR, ec)) {
    const std::string n = e.path().filename().string();
    if (n.size() <= 5 || n.compare(n.size() - 5, 5, ".fits") != 0) continue;
    double fov = 0;
    const size_t deg = n.rfind("deg.fits");
    if (deg != std::string::npos) {
      size_t start = n.find_last_of('_', deg);
      if (start != std::string::npos) fov = std::atof(n.substr(start + 1, deg - start - 1).c_str());
    }
    candidates.emplace_back(std::fabs(fov - 1.0), e.path().string()); // nearest a degree first
  }
  std::sort(candidates.begin(), candidates.end());
  if (candidates.empty()) {
    std::printf("SKIPPED: no .fits image in %s\n", ASTAP_TEST_CORPUS_DIR);
    return 0;
  }

  const std::string out_base =
      (std::filesystem::path(ASTAP_TEST_TMP_DIR) / "service_api_tests").string();

  SolveService svc;
  if (!svc.load(SolveServiceSettings{})) {
    std::printf("SKIPPED: no index could be loaded or built\n");
    return 0;
  }

  // Settle on an image by solving it, rather than by trusting its name.
  std::string image;
  SolveOutcome from_file;
  for (size_t i = 0; i < candidates.size() && i < 4; i++) {
    SolveRequest req;
    req.filename = candidates[i].second;
    req.output_base = out_base;
    from_file = svc.solve(req);
    if (from_file.solved) {
      image = req.filename;
      break;
    }
  }
  if (image.empty()) {
    std::printf("SKIPPED: none of the first images in %s solved\n", ASTAP_TEST_CORPUS_DIR);
    return 0;
  }
  std::printf("image: %s\n", image.c_str());

  Header file_head;
  ImageArray img;
  const ImageLoadResult lr = load_image(image, file_head, img);
  if (!lr.ok) {
    std::printf("FAIL: could not read %s: %s\n", image.c_str(), lr.error.c_str());
    return 1;
  }

  // --- the index solver ------------------------------------------------------
  {
    check(from_file.solved, "index: the file path solves");

    // The same pixels, with a header that never saw a file.
    Header head = header_for_image(img);
    check(head.width == img.width() && head.height == img.height(),
          "header_for_image carries the dimensions");
    const SolveOutcome from_array = svc.solve_image(img, head, SolveParams{});
    check(from_array.solved, "index: the array path solves");

    if (from_file.solved && from_array.solved) {
      same_position(from_file.head, from_array.head, "index");
      check(from_file.stars == from_array.stars, "index: same stars detected");
      check(from_file.nr_inliers == from_array.nr_inliers, "index: same consensus");
      check(from_file.tier_density == from_array.tier_density, "index: same depth tier");
    }
    // The array path is not allowed to have written anything of its own.
    check(!std::filesystem::exists(out_base + ".wcs"),
          "index: no .wcs unless one was asked for");
  }

  // --- the spiral port -------------------------------------------------------
  {
    SpiralService port;
    if (!port.load(SpiralServiceSettings{})) {
      std::printf("FAIL: the port found no database the index solver had just used\n");
      return 1;
    }

    SpiralRequest req;
    req.filename = image;
    req.output_base = out_base;
    const SolveOutcome port_file = port.solve(req);
    check(port_file.solved, "port: the file path solves");

    // The port's array path gets the header the file path had. That is not a
    // convenience: unlike the index solver, which takes the field size from its
    // parameters, the port reads the plate scale out of the header to decide
    // what field size to search for. Handing it header_for_image() instead is a
    // different question — a blind solve — and would rightly give a different
    // answer, so it would not test what this is testing.
    Header head = file_head;
    const SolveOutcome port_array = port.solve_image(img, head, SpiralParams{});
    check(port_array.solved, "port: the array path solves");

    if (port_file.solved && port_array.solved)
      same_position(port_file.head, port_array.head, "port");

    // And the two solvers have to agree with each other, which is the check
    // that neither service is quietly solving a different image.
    if (port_file.solved && from_file.solved)
      same_position(port_file.head, from_file.head, "port against index");
  }

  std::filesystem::remove(out_base + ".ini");
  std::printf(failures ? "%d failure(s)\n" : "service api tests passed\n", failures);
  return failures ? 1 : 0;
}
