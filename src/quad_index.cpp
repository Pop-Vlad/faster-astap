#include "astap/quad_index.h"

#include <algorithm>
#include <atomic>
#include <cmath>

#include "astap/astro_math.h"
#include "astap/parallel.h"
#include "astap/quads.h"

namespace astap {
  namespace {
    // Tile geometry, mirroring the ring tables of the database reader.
    const int kCells290[18] = {1, 4, 8, 12, 16, 20, 24, 28, 32, 32, 28, 24, 20, 16, 12, 8, 4, 1};
    const int kCells1476[36] = {
      1, 3, 9, 15, 21, 27, 33, 38, 43, 48, 52, 56,
      60, 63, 65, 67, 68, 69, 69, 68, 67, 65, 63, 60,
      56, 52, 48, 43, 38, 33, 27, 21, 15, 9, 3, 1
    };

    int tile_count(int database_type) { return database_type == kDatabase290 ? 290 : 1476; }

    // Centre of tile `area` (1 based) and the tile's angular size, so the tile
    // can be projected into its own tangent plane.
    void tile_geometry(int database_type, int area, double &ra, double &dec, double &size_deg) {
      const int *cells = database_type == kDatabase290 ? kCells290 : kCells1476;
      const int nrings = database_type == kDatabase290 ? 18 : 36;
      const double ring_h = 180.0 / nrings;

      int remaining = area;
      for (int ring = 0; ring < nrings; ring++) {
        if (remaining <= cells[ring]) {
          const double dec_lo = -90.0 + ring * ring_h;
          dec = (dec_lo + ring_h * 0.5) * kPi / 180;
          ra = (remaining - 0.5) * 2 * kPi / cells[ring];
          // The wider of the two extents, used only to size the search field.
          const double w = 360.0 / cells[ring] * std::cos(dec);
          size_deg = std::max(ring_h, w);
          return;
        }
        remaining -= cells[ring];
      }
      ra = 0;
      dec = 0;
      size_deg = ring_h;
    }
  } // namespace

  int QuadIndex::bin_of(float v) const {
    int b = static_cast<int>(v / static_cast<float>(settings_.quad_tolerance));
    if (b < 0) b = 0;
    if (b >= nbins_) b = nbins_ - 1;
    return b;
  }

  uint32_t QuadIndex::cell_of(int b0, int b1, int b2) const {
    return static_cast<uint32_t>((b0 * nbins_ + b1) * nbins_ + b2);
  }

  size_t QuadIndex::bytes() const {
    return ratio_.size() * sizeof(float) + (d1_.size() + ra_.size() + dec_.size()) * sizeof(double) +
           (cell_start_.size() + items_.size()) * sizeof(uint32_t);
  }

  bool QuadIndex::build(StarDatabase &db, const QuadIndexSettings &s,
                        const std::function<void(double)> &progress) {
    settings_ = s;
    ratio_.clear();
    d1_.clear();
    ra_.clear();
    dec_.clear();

    const int ntiles = tile_count(db.database_type());
    const double radius_rad = s.radius_deg * kPi / 180;

    // Tiles are independent, so they are built in parallel and merged in tile
    // order afterwards, which keeps the index deterministic.
    std::vector<std::vector<float>> t_ratio(ntiles);
    std::vector<std::vector<double>> t_d1(ntiles), t_ra(ntiles), t_dec(ntiles);
    std::vector<char> t_ok(ntiles, 1);

    std::atomic<int> done{0};
    parallel_for(0, static_cast<size_t>(ntiles), [&](size_t ti, unsigned) {
      const int area = static_cast<int>(ti) + 1;
      double tra, tdec, tsize;
      tile_geometry(db.database_type(), area, tra, tdec, tsize);

      if (s.radius_deg < 180) {
        double sep;
        ang_sep(tra, tdec, s.centre_ra, s.centre_dec, sep);
        if (sep > radius_rad + tsize * kPi / 180) return; // tile outside the cap
      }

      // How many stars this tile should contribute at the requested density.
      const double tile_area_deg2 = 41253.0 / ntiles;
      const int want = std::max(8, static_cast<int>(s.star_density * tile_area_deg2));

      // Each thread needs its own reader: a StarDatabase owns a file handle and
      // a read cursor, so it cannot be shared or copied.
      StarDatabase reader;
      reader.configure(db.path(), db.name(), db.database_type());
      if (!reader.open_area(tdec, area)) {
        t_ok[ti] = 0;
        return;
      }

      // A field a little larger than the tile so the whole tile is covered.
      const double field = tsize * 1.5 * kPi / 180;
      std::vector<double> sra, sdec;
      sra.reserve(want);
      sdec.reserve(want);
      {
        double ra2 = 0, dec2 = 0, mag = 0, bv = 0;
        while (static_cast<int>(sra.size()) < want &&
               reader.read_star(tra, tdec, field, ra2, dec2, mag, bv)) {
          sra.push_back(ra2);
          sdec.push_back(dec2);
        }
      }
      if (sra.size() < 8) return;

      // Project into the tile's own tangent plane, in arcseconds.
      RowList stars(3, sra.size());
      for (size_t i = 0; i < sra.size(); i++) {
        double x, y;
        equatorial_standard(tra, tdec, sra[i], sdec[i], 1, x, y);
        stars(0, i) = x;
        stars(1, i) = y;
        stars(2, i) = 100;
      }

      // The image side uses find_quads with its own star count; here the tile
      // is dense, so the regular three-nearest path applies.
      RowList quads;
      find_quads(1000 /* forces the regular path */, stars, quads);

      std::vector<float> &r = t_ratio[ti];
      std::vector<double> &o1 = t_d1[ti];
      std::vector<double> &ora = t_ra[ti];
      std::vector<double> &odec = t_dec[ti];
      r.reserve(quads.count() * 5);
      o1.reserve(quads.count());
      ora.reserve(quads.count());
      odec.reserve(quads.count());

      for (size_t q = 0; q < quads.count(); q++) {
        for (int k = 1; k <= 5; k++) r.push_back(static_cast<float>(quads(k, q)));
        o1.push_back(quads(0, q));
        // Quad centre back to absolute coordinates.
        double qra, qdec;
        standard_equatorial(tra, tdec, quads(6, q), quads(7, q), 1, qra, qdec);
        ora.push_back(qra);
        odec.push_back(qdec);
      }

      const int d = ++done;
      if (progress && (d % 32 == 0)) progress(static_cast<double>(d) / ntiles);
    });

    size_t total = 0;
    for (int t = 0; t < ntiles; t++) total += t_d1[t].size();
    ratio_.reserve(total * 5);
    d1_.reserve(total);
    ra_.reserve(total);
    dec_.reserve(total);
    for (int t = 0; t < ntiles; t++) {
      ratio_.insert(ratio_.end(), t_ratio[t].begin(), t_ratio[t].end());
      d1_.insert(d1_.end(), t_d1[t].begin(), t_d1[t].end());
      ra_.insert(ra_.end(), t_ra[t].begin(), t_ra[t].end());
      dec_.insert(dec_.end(), t_dec[t].begin(), t_dec[t].end());
    }
    if (d1_.empty()) return false;

    // Bin on the first three ratios. They lie in (0,1], so the bin count is
    // 1/tolerance plus a slot for the value 1.
    nbins_ = static_cast<int>(1.0 / s.quad_tolerance) + 2;
    const size_t ncells = static_cast<size_t>(nbins_) * nbins_ * nbins_;
    cell_start_.assign(ncells + 1, 0);
    std::vector<uint32_t> cell_of_quad(d1_.size());
    for (size_t i = 0; i < d1_.size(); i++) {
      const float *r = ratios(i);
      const uint32_t c = cell_of(bin_of(r[0]), bin_of(r[1]), bin_of(r[2]));
      cell_of_quad[i] = c;
      cell_start_[c + 1]++;
    }
    for (size_t c = 0; c < ncells; c++) cell_start_[c + 1] += cell_start_[c];
    items_.resize(d1_.size());
    std::vector<uint32_t> fill(cell_start_.begin(), cell_start_.end() - 1);
    for (size_t i = 0; i < d1_.size(); i++) items_[fill[cell_of_quad[i]]++] = static_cast<uint32_t>(i);

    if (progress) progress(1.0);
    return true;
  }

  void QuadIndex::query(const float *r, std::vector<uint32_t> &hits) const {
    if (items_.empty()) return;
    const float tol = static_cast<float>(settings_.quad_tolerance);
    const int b0 = bin_of(r[0]), b1 = bin_of(r[1]), b2 = bin_of(r[2]);

    for (int i = -1; i <= 1; i++) {
      const int c0 = b0 + i;
      if (c0 < 0 || c0 >= nbins_) continue;
      for (int j = -1; j <= 1; j++) {
        const int c1 = b1 + j;
        if (c1 < 0 || c1 >= nbins_) continue;
        for (int k = -1; k <= 1; k++) {
          const int c2 = b2 + k;
          if (c2 < 0 || c2 >= nbins_) continue;
          const uint32_t c = cell_of(c0, c1, c2);
          for (uint32_t p = cell_start_[c]; p < cell_start_[c + 1]; p++) {
            const uint32_t q = items_[p];
            const float *s = ratios(q);
            if (std::fabs(s[0] - r[0]) <= tol && std::fabs(s[1] - r[1]) <= tol &&
                std::fabs(s[2] - r[2]) <= tol && std::fabs(s[3] - r[3]) <= tol &&
                std::fabs(s[4] - r[4]) <= tol)
              hits.push_back(q);
          }
        }
      }
    }
  }
} // namespace astap
