// Verifies that whatever QuadBatchBuilder is available reproduces find_quads()
// bit for bit. This is the gate the SYCL backend has to pass.
#include <cmath>
#include <cstdio>
#include <cstring>
#include <random>
#include <vector>

#include "astap/quad_batch.h"
#include "astap/quads.h"

using namespace astap;

static int failures = 0;

static bool bit_equal(double a, double b) {
  return std::memcmp(&a, &b, sizeof(double)) == 0;
}

// Compares a builder against per-position find_quads for a batch of lists.
static void compare(QuadBatchBuilder& b, int nrstars_image, size_t npos, size_t nstars,
                    unsigned seed) {
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> u(-3000, 3000);

  std::vector<RowList> a_lists(npos), b_lists(npos);
  for (size_t p = 0; p < npos; p++) {
    a_lists[p].resize(3, nstars);
    for (size_t i = 0; i < nstars; i++) {
      a_lists[p](0, i) = u(rng);
      a_lists[p](1, i) = u(rng);
      a_lists[p](2, i) = 100;
    }
    b_lists[p] = a_lists[p];
  }

  std::vector<RowList> ref(npos), got;
  for (size_t p = 0; p < npos; p++) find_quads(nrstars_image, a_lists[p], ref[p]);
  b.build(nrstars_image, b_lists, got);

  size_t total = 0, mismatched = 0;
  for (size_t p = 0; p < npos; p++) {
    if (ref[p].count() != got[p].count()) {
      std::printf("FAIL: position %zu quad count %zu vs %zu\n", p, ref[p].count(), got[p].count());
      failures++;
      return;
    }
    for (size_t i = 0; i < ref[p].count(); i++)
      for (int r = 0; r < 8; r++) {
        total++;
        if (!bit_equal(ref[p](r, i), got[p](r, i))) mismatched++;
      }
  }
  if (mismatched) {
    std::printf("FAIL: %zu of %zu values differ (image stars %d, %zu stars/position)\n", mismatched,
                total, nrstars_image, nstars);
    failures++;
  } else {
    std::printf("ok  : %-22s image stars %3d, %3zu positions x %3zu stars -> %zu values identical\n",
                b.backend().c_str(), nrstars_image, npos, nstars, total);
  }
}

int main() {
  struct Case { int nrstars_image; size_t nstars; const char* what; };
  const Case cases[] = {
      {7, 20, "mode 7"},   {7, 60, "mode 7"},  {25, 40, "mode 6"},
      {40, 80, "mode 5"},  {124, 94, "find_quads"}, {200, 300, "find_quads sorted"},
  };

  auto cpu = make_cpu_quad_builder();
  for (const Case& c : cases) compare(*cpu, c.nrstars_image, 24, c.nstars, 7);

  auto gpu = make_sycl_quad_builder();
  if (!gpu) {
    std::printf("\n(no SYCL device available, backend not exercised)\n");
  } else {
    std::printf("\n");
    for (const Case& c : cases) compare(*gpu, c.nrstars_image, 24, c.nstars, 7);
    // A batch with mixed star counts must still come out right.
    std::printf("\n");
    compare(*gpu, 7, 40, 12, 99);
  }

  if (failures) {
    std::printf("\n%d check(s) FAILED\n", failures);
    return 1;
  }
  std::printf("\nall quad batch checks passed\n");
  return 0;
}
