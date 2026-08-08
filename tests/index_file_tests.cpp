// Verifies the on-disk index layout: every array starts on a boundary a mapped
// file could address it at, and the values survive a save/load round trip.
//
// The offsets are re-derived here from the header rather than taken from the
// writer, so that a writer and a reader agreeing on a layout no mapping could
// use is still a failure. Version 2 packed the arrays end to end, which put
// d1/ra/dec at a 4 byte offset whenever a tier held an odd number of quads.
//
// Needs a star database to build an index from; without one the test reports
// that it was skipped and passes, as there is nothing it can check.

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "astap/quad_index.h"
#include "astap/star_database.h"

using namespace astap;

static int failures = 0;

static void check(bool ok, const std::string &what) {
  if (!ok) {
    std::printf("FAIL: %s\n", what.c_str());
    failures++;
  }
}

namespace {
  // The layout rule of the format, written out independently of quad_index.cpp.
  constexpr uint64_t kAlign = 64;

  uint64_t align_up(uint64_t x) { return (x + kAlign - 1) / kAlign * kAlign; }

  struct Tier {
    double density;
    uint64_t nquads;
    uint64_t offset;
  };

  template<typename T>
  T read_pod(std::istream &f) {
    T v{};
    f.read(reinterpret_cast<char *>(&v), sizeof(T));
    return v;
  }

  // Parses the header and returns the tier table, or an empty vector.
  std::vector<Tier> parse_header(std::istream &f, uint32_t &version, uint32_t &nbins) {
    char magic[8];
    f.read(magic, sizeof(magic));
    if (!f || std::memcmp(magic, "ASTAPQIX", 8) != 0) return {};
    version = read_pod<uint32_t>(f);
    read_pod<uint32_t>(f); // byte order mark
    read_pod<uint32_t>(f); // database type
    const uint32_t ntiers = read_pod<uint32_t>(f);
    nbins = read_pod<uint32_t>(f);
    read_pod<uint32_t>(f); // reserved
    for (int i = 0; i < 4; i++) read_pod<double>(f); // tolerance, centre, radius
    std::vector<Tier> tiers(ntiers);
    for (uint32_t k = 0; k < ntiers; k++) {
      tiers[k].density = read_pod<double>(f);
      tiers[k].nquads = read_pod<uint64_t>(f);
      tiers[k].offset = read_pod<uint64_t>(f);
    }
    return f ? tiers : std::vector<Tier>();
  }
} // namespace

int main() {
  // Any database will do: the layout does not depend on which one.
  StarDatabase db;
  bool have_db = false;
  for (const std::string &dir: default_database_directories()) {
    std::string d = dir;
    if (!d.empty() && d.back() != '/' && d.back() != '\\') d += '/';
    if (db.select(d, "auto", 1.0)) {
      have_db = true;
      std::printf("database: %s in %s\n", db.name().c_str(), db.path().c_str());
      break;
    }
  }
  if (!have_db) {
    std::printf("SKIPPED: no star database found, nothing to build an index from\n");
    return 0;
  }

  // A small cap keeps the build to a second or two. Two sparse tiers are enough
  // to exercise the tier-to-tier stride as well as the arrays within one tier.
  QuadIndexSettings s;
  s.quad_tolerance = 0.007;
  s.centre_ra = 0;
  s.centre_dec = 0;
  s.radius_deg = 6;
  std::vector<QuadIndex> built;
  if (!build_tiers(db, s, {0.5, 2.0}, built) || built.empty()) {
    std::printf("FAIL: could not build an index from %s\n", db.path().c_str());
    return 1;
  }
  for (const QuadIndex &ix: built)
    std::printf("built tier %g stars/deg^2: %zu quads\n", ix.settings().star_density, ix.size());

  const std::string path =
      (std::filesystem::path(ASTAP_TEST_TMP_DIR) / "index_file_tests.qix").string();
  if (!save_index_file(path, built)) {
    std::printf("FAIL: could not write %s\n", path.c_str());
    return 1;
  }

  const uint64_t file_size = std::filesystem::file_size(path);

  // --- the layout, re-derived from the header ------------------------------
  {
    std::ifstream f(path, std::ios::binary);
    uint32_t version = 0, nbins = 0;
    const std::vector<Tier> tiers = parse_header(f, version, nbins);
    check(!tiers.empty(), "header parses");
    check(version >= 4, "version is at least 4 (padded arrays, stored grid)");
    check(nbins >= 2, "the header carries a plausible bin count");
    check(tiers.size() == built.size(), "tier count matches");
    const uint64_t ncells = static_cast<uint64_t>(nbins) * nbins * nbins;

    for (size_t k = 0; k < tiers.size() && k < built.size(); k++) {
      const uint64_t n = tiers[k].nquads;
      const std::string at = "tier " + std::to_string(k);
      check(n == built[k].size(), at + " quad count matches");

      const uint64_t ratio = tiers[k].offset;
      const uint64_t d1 = align_up(ratio + n * 5 * sizeof(float));
      const uint64_t ra = align_up(d1 + n * sizeof(double));
      const uint64_t dec = align_up(ra + n * sizeof(double));
      const uint64_t cell_start = align_up(dec + n * sizeof(double));
      const uint64_t items = align_up(cell_start + (ncells + 1) * sizeof(uint32_t));

      // The invariant that makes the file mappable.
      check(ratio % kAlign == 0, at + " ratio array is aligned");
      check(d1 % kAlign == 0, at + " d1 array is aligned");
      check(ra % kAlign == 0, at + " ra array is aligned");
      check(dec % kAlign == 0, at + " dec array is aligned");
      check(cell_start % kAlign == 0, at + " cell_start array is aligned");
      check(items % kAlign == 0, at + " items array is aligned");
      check(items + n * sizeof(uint32_t) <= file_size, at + " lies inside the file");

      // The offsets are only right if the data is actually there. Reading one
      // double at the derived position of the last element ties the two.
      if (n > 0 && items + n * sizeof(uint32_t) <= file_size) {
        f.seekg(static_cast<std::streamoff>(d1 + (n - 1) * sizeof(double)));
        check(read_pod<double>(f) == built[k].d1(n - 1), at + " d1 lands where derived");
        f.seekg(static_cast<std::streamoff>(ra + (n - 1) * sizeof(double)));
        check(read_pod<double>(f) == built[k].centre_ra(n - 1), at + " ra lands where derived");
        f.seekg(static_cast<std::streamoff>(dec + (n - 1) * sizeof(double)));
        check(read_pod<double>(f) == built[k].centre_dec(n - 1), at + " dec lands where derived");

        // The CSR grid is only consistent if its last offset is the quad count:
        // every quad has to sit in exactly one cell.
        f.seekg(static_cast<std::streamoff>(cell_start + ncells * sizeof(uint32_t)));
        check(read_pod<uint32_t>(f) == n, at + " the stored grid indexes every quad");
        f.seekg(static_cast<std::streamoff>(cell_start));
        check(read_pod<uint32_t>(f) == 0, at + " the stored grid starts at zero");
      }
    }
  }

  // --- the round trip -------------------------------------------------------
  {
    std::vector<QuadIndex> loaded;
    std::string err;
    if (!load_index_file(path, loaded, 0, 0, &err)) {
      std::printf("FAIL: could not read it back: %s\n", err.c_str());
      failures++;
    } else {
      check(loaded.size() == built.size(), "round trip tier count");
      for (size_t k = 0; k < loaded.size() && k < built.size(); k++) {
        check(loaded[k].mapped(), "the loaded tier reads out of the mapping");
        check(loaded[k].size() == built[k].size(), "round trip quad count");

        // The grid came off disk rather than being rebuilt, so what matters is
        // that it still answers exactly as the one built in memory does.
        bool same_hits = true;
        for (size_t i = 0; i < loaded[k].size() && same_hits; i += 97) {
          std::vector<uint32_t> a, b;
          built[k].query(built[k].ratios(i), a);
          loaded[k].query(loaded[k].ratios(i), b);
          same_hits = a == b && !a.empty();
        }
        check(same_hits, "the stored grid answers queries as the built one does");
        check(loaded[k].settings().star_density == built[k].settings().star_density,
              "round trip density");
        bool same = true;
        for (size_t i = 0; i < loaded[k].size() && same; i++) {
          same = loaded[k].d1(i) == built[k].d1(i) &&
                 loaded[k].centre_ra(i) == built[k].centre_ra(i) &&
                 loaded[k].centre_dec(i) == built[k].centre_dec(i) &&
                 std::memcmp(loaded[k].ratios(i), built[k].ratios(i), 5 * sizeof(float)) == 0;
        }
        check(same, "round trip quad values");
      }
    }
  }

  // A file written by an older layout has to be refused, not misread.
  {
    const std::string stale =
        (std::filesystem::path(ASTAP_TEST_TMP_DIR) / "index_file_tests_v2.qix").string();
    std::filesystem::copy_file(path, stale, std::filesystem::copy_options::overwrite_existing);
    std::fstream f(stale, std::ios::binary | std::ios::in | std::ios::out);
    f.seekp(8);
    const uint32_t old_version = 2;
    f.write(reinterpret_cast<const char *>(&old_version), sizeof(old_version));
    f.close();

    std::vector<QuadIndex> loaded;
    std::string err;
    check(!load_index_file(stale, loaded, 0, 0, &err), "an older version is refused");
    std::filesystem::remove(stale);
  }

  std::filesystem::remove(path);
  std::printf(failures ? "%d failure(s)\n" : "index file tests passed\n", failures);
  return failures ? 1 : 0;
}
