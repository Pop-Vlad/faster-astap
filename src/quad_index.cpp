#include "astap/quad_index.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <istream>
#include <ostream>
#include <system_error>

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

    // One tile's stars, brightest first, already projected into the tile's own
    // tangent plane in arcseconds. Tile-local projection is fine: over a 5.14
    // degree tile the differential tangent-plane distortion is about 6e-4
    // relative, against a quad tolerance of 7e-3.
    struct TileStars {
      double ra = 0, dec = 0;  // tile centre
      RowList projected;       // x, y, magnitude placeholder
    };

    // Reads at most `want` stars of one tile. Returns false when the tile could
    // not be opened; an empty result is not an error, only a sparse tile.
    bool read_tile(const StarDatabase &db, int area, int want, TileStars &out, bool &io_error) {
      double tsize;
      tile_geometry(db.database_type(), area, out.ra, out.dec, tsize);

      // Each thread needs its own reader: a StarDatabase owns a file handle and
      // a read cursor, so it cannot be shared or copied.
      StarDatabase reader;
      reader.configure(db.path(), db.name(), db.database_type());
      if (!reader.open_area(out.dec, area)) {
        io_error = true;
        return false;
      }

      // A field a little larger than the tile so the whole tile is covered.
      const double field = tsize * 1.5 * kPi / 180;
      std::vector<double> sra, sdec;
      sra.reserve(want);
      sdec.reserve(want);
      double ra2 = 0, dec2 = 0, mag = 0, bv = 0;
      while (static_cast<int>(sra.size()) < want &&
             reader.read_star(out.ra, out.dec, field, ra2, dec2, mag, bv)) {
        sra.push_back(ra2);
        sdec.push_back(dec2);
      }
      if (sra.size() < 4) return false;

      out.projected.resize(3, sra.size());
      for (size_t i = 0; i < sra.size(); i++) {
        double x, y;
        equatorial_standard(out.ra, out.dec, sra[i], sdec[i], 1, x, y);
        out.projected(0, i) = x;
        out.projected(1, i) = y;
        out.projected(2, i) = 100;
      }
      return true;
    }

    // Quads from the brightest `take` stars of a tile, appended to the per-tile
    // output buffers with their centres back in absolute coordinates.
    void quads_from_tile(const TileStars &t, size_t take, bool many, std::vector<float> &ratio,
                         std::vector<double> &d1, std::vector<double> &ra,
                         std::vector<double> &dec) {
      take = std::min(take, t.projected.count());
      if (take < 8) return;

      // find_quads sorts its input in place, so each tier works on its own copy.
      RowList stars(3, take);
      for (size_t i = 0; i < take; i++) {
        stars(0, i) = t.projected(0, i);
        stars(1, i) = t.projected(1, i);
        stars(2, i) = t.projected(2, i);
      }

      RowList quads;
      if (many)
        find_many_quads(stars, quads, 6);  // 15 per star, see many_quads_below_density
      else
        find_quads(1000 /* forces the regular three-nearest path */, stars, quads);

      ratio.reserve(ratio.size() + quads.count() * 5);
      d1.reserve(d1.size() + quads.count());
      ra.reserve(ra.size() + quads.count());
      dec.reserve(dec.size() + quads.count());
      for (size_t q = 0; q < quads.count(); q++) {
        for (int k = 1; k <= 5; k++) ratio.push_back(static_cast<float>(quads(k, q)));
        d1.push_back(quads(0, q));
        double qra, qdec;
        standard_equatorial(t.ra, t.dec, quads(6, q), quads(7, q), 1, qra, qdec);
        ra.push_back(qra);
        dec.push_back(qdec);
      }
    }

    // How many stars a tile should contribute at a given density.
    int stars_wanted(int ntiles, double density) {
      const double tile_area_deg2 = 41253.0 / ntiles;
      return std::max(8, static_cast<int>(density * tile_area_deg2));
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
    // What the index addresses, whether that memory is this process's or a
    // mapped file's. Reported as the index size, so it should not change with
    // where the quads happen to live.
    const uint64_t ncells = static_cast<uint64_t>(nbins_) * nbins_ * nbins_;
    return nquads_ * (5 * sizeof(float) + 3 * sizeof(double) + sizeof(uint32_t)) +
           static_cast<size_t>(nbins_ > 0 ? (ncells + 1) * sizeof(uint32_t) : 0);
  }

  void QuadIndex::repoint() {
    // A mapped index has every view aimed into the mapping by the loader, and
    // holds no vectors that could have moved.
    if (map_) return;
    nquads_ = d1_.size();
    v_ratio_ = ratio_.data();
    v_d1_ = d1_.data();
    v_ra_ = ra_.data();
    v_dec_ = dec_.data();
    v_cell_start_ = cell_start_.data();
    v_items_ = items_.data();
    nitems_ = items_.size();
  }

  QuadIndex &QuadIndex::operator=(const QuadIndex &o) {
    if (this != &o) {
      settings_ = o.settings_;
      map_ = o.map_;
      ratio_ = o.ratio_;
      d1_ = o.d1_;
      ra_ = o.ra_;
      dec_ = o.dec_;
      cell_start_ = o.cell_start_;
      items_ = o.items_;
      nbins_ = o.nbins_;
      nquads_ = o.nquads_;
      // Carried over as they are for a mapped index, which now shares the same
      // mapping; overwritten by repoint() for an owning one, whose vectors this
      // has just copied.
      v_ratio_ = o.v_ratio_;
      v_d1_ = o.v_d1_;
      v_ra_ = o.v_ra_;
      v_dec_ = o.v_dec_;
      v_cell_start_ = o.v_cell_start_;
      v_items_ = o.v_items_;
      nitems_ = o.nitems_;
      repoint();
    }
    return *this;
  }

  QuadIndex &QuadIndex::operator=(QuadIndex &&o) noexcept {
    if (this != &o) {
      settings_ = o.settings_;
      map_ = std::move(o.map_);
      ratio_ = std::move(o.ratio_);
      d1_ = std::move(o.d1_);
      ra_ = std::move(o.ra_);
      dec_ = std::move(o.dec_);
      cell_start_ = std::move(o.cell_start_);
      items_ = std::move(o.items_);
      nbins_ = o.nbins_;
      nquads_ = o.nquads_;
      v_ratio_ = o.v_ratio_;
      v_d1_ = o.v_d1_;
      v_ra_ = o.v_ra_;
      v_dec_ = o.v_dec_;
      v_cell_start_ = o.v_cell_start_;
      v_items_ = o.v_items_;
      nitems_ = o.nitems_;
      repoint();
    }
    return *this;
  }

  void QuadIndex::finalise() {
    // So that size() and ratios() below read from whichever storage this index
    // holds: the loader has mapped the quads, the build path has just filled the
    // vectors, and the grid is built the same way from either.
    repoint();

    // Bin on the first three ratios. They lie in (0,1], so the bin count is
    // 1/tolerance plus a slot for the value 1.
    nbins_ = static_cast<int>(1.0 / settings_.quad_tolerance) + 2;
    const size_t ncells = static_cast<size_t>(nbins_) * nbins_ * nbins_;
    cell_start_.assign(ncells + 1, 0);
    std::vector<uint32_t> cell_of_quad(nquads_);
    for (size_t i = 0; i < nquads_; i++) {
      const float *r = ratios(i);
      const uint32_t c = cell_of(bin_of(r[0]), bin_of(r[1]), bin_of(r[2]));
      cell_of_quad[i] = c;
      cell_start_[c + 1]++;
    }
    for (size_t c = 0; c < ncells; c++) cell_start_[c + 1] += cell_start_[c];
    items_.resize(nquads_);
    std::vector<uint32_t> fill(cell_start_.begin(), cell_start_.end() - 1);
    for (size_t i = 0; i < nquads_; i++)
      items_[fill[cell_of_quad[i]]++] = static_cast<uint32_t>(i);

    repoint(); // the grid vectors have just been rebuilt
  }

  bool QuadIndex::build(StarDatabase &db, const QuadIndexSettings &s,
                        const std::function<void(double)> &progress) {
    std::vector<QuadIndex> one;
    if (!build_tiers(db, s, {s.star_density}, one, progress)) return false;
    *this = std::move(one[0]);
    return true;
  }

  bool build_tiers(StarDatabase &db, const QuadIndexSettings &base,
                   const std::vector<double> &densities, std::vector<QuadIndex> &out,
                   const std::function<void(double)> &progress) {
    out.clear();
    if (densities.empty()) return false;

    const int ntiles = tile_count(db.database_type());
    const int ntiers = static_cast<int>(densities.size());
    const double radius_rad = base.radius_deg * kPi / 180;
    double deepest = 0;
    for (double d : densities) deepest = std::max(deepest, d);
    const int want_deepest = stars_wanted(ntiles, deepest);

    // Tiles are independent, so they are built in parallel and merged in tile
    // order afterwards, which keeps the index deterministic. Every tier gets its
    // own per-tile buffer.
    std::vector<std::vector<std::vector<float>>> t_ratio(
        ntiers, std::vector<std::vector<float>>(ntiles));
    std::vector<std::vector<std::vector<double>>> t_d1(
        ntiers, std::vector<std::vector<double>>(ntiles));
    std::vector<std::vector<std::vector<double>>> t_ra(t_d1), t_dec(t_d1);

    std::atomic<int> done{0};
    std::atomic<int> io_errors{0};
    parallel_for(0, static_cast<size_t>(ntiles), [&](size_t ti, unsigned) {
      const int area = static_cast<int>(ti) + 1;
      double tra, tdec, tsize;
      tile_geometry(db.database_type(), area, tra, tdec, tsize);

      if (base.radius_deg < 180) {
        double sep;
        ang_sep(tra, tdec, base.centre_ra, base.centre_dec, sep);
        if (sep > radius_rad + tsize * kPi / 180) return; // tile outside the cap
      }

      // Read and project once, at the deepest tier; the others take a prefix.
      TileStars tile;
      bool io_error = false;
      const bool have = read_tile(db, area, want_deepest, tile, io_error);
      if (io_error) io_errors++;
      if (have) {
        for (int k = 0; k < ntiers; k++)
          quads_from_tile(tile, static_cast<size_t>(stars_wanted(ntiles, densities[k])),
                          densities[k] < base.many_quads_below_density, t_ratio[k][ti],
                          t_d1[k][ti], t_ra[k][ti], t_dec[k][ti]);
      }

      const int d = ++done;
      if (progress && (d % 32 == 0)) progress(static_cast<double>(d) / ntiles);
    });

    out.resize(ntiers);
    for (int k = 0; k < ntiers; k++) {
      QuadIndex &ix = out[k];
      ix.settings_ = base;
      ix.settings_.star_density = densities[k];
      size_t total = 0;
      for (int t = 0; t < ntiles; t++) total += t_d1[k][t].size();
      ix.ratio_.reserve(total * 5);
      ix.d1_.reserve(total);
      ix.ra_.reserve(total);
      ix.dec_.reserve(total);
      for (int t = 0; t < ntiles; t++) {
        ix.ratio_.insert(ix.ratio_.end(), t_ratio[k][t].begin(), t_ratio[k][t].end());
        ix.d1_.insert(ix.d1_.end(), t_d1[k][t].begin(), t_d1[k][t].end());
        ix.ra_.insert(ix.ra_.end(), t_ra[k][t].begin(), t_ra[k][t].end());
        ix.dec_.insert(ix.dec_.end(), t_dec[k][t].begin(), t_dec[k][t].end());
      }
      if (ix.d1_.empty()) {
        out.clear();
        return false;
      }
      ix.finalise();
    }

    if (progress) progress(1.0);
    return true;
  }

  namespace {
    // "ASTAP quad index". The byte order mark and the version both have to
    // match: a file read with the wrong endianness or an older layout would not
    // fail, it would match against nonsense.
    const char kIndexMagic[8] = {'A', 'S', 'T', 'A', 'P', 'Q', 'I', 'X'};
    // 4: the bin grid is stored alongside the quads, and every array starts on a
    // kAlign boundary, so loading is a header parse and a mapping and the pages
    // a query never reaches are never read.
    //
    // Storing the grid is what makes the mapping worth having. Version 3 rebuilt
    // it on load, which meant scanning every quad's ratios — the one access
    // pattern that faults in the entire file, and measurably slower through a
    // mapping than the bulk read it replaced. It costs about 15% more on disk.
    //
    // Version 2 additionally packed the arrays end to end, which left d1, ra and
    // dec at a 4 byte offset whenever a tier held an odd number of quads:
    // readable through a copying read, not addressable as doubles.
    constexpr uint32_t kIndexVersion = 4;
    constexpr uint32_t kByteOrderMark = 0x01020304u;

    // 8 would be enough to make the doubles naturally aligned. A cache line
    // costs at most 63 bytes per array on a file of gigabytes and keeps the
    // start of one array off the last line of the one before it.
    constexpr uint64_t kAlign = 64;

    uint64_t align_up(uint64_t x) { return (x + kAlign - 1) / kAlign * kAlign; }

    template <typename T>
    void put(std::ostream &f, const T &v) {
      f.write(reinterpret_cast<const char *>(&v), sizeof(T));
    }

    template <typename T>
    void get(std::istream &f, T &v) {
      f.read(reinterpret_cast<char *>(&v), sizeof(T));
    }

    // Written from the index's views rather than its vectors, so that a tier
    // read from one cache file can be written back out to another.
    template <typename T>
    void put_view(std::ostream &f, const T *p, size_t n) {
      if (n) f.write(reinterpret_cast<const char *>(p), static_cast<std::streamsize>(n * sizeof(T)));
    }

    // Fixed size header, so tier blob offsets can be computed before writing.
    struct TierEntry {
      double density;
      uint64_t nquads;
      uint64_t offset;
    };

    // magic, version, byte order, database type, tier count, bin count, a spare
    // word, then the four doubles and the tier table. 64 bytes before the table,
    // which is the alignment, so a file with no padding at all is still legal.
    size_t header_bytes(size_t ntiers) {
      return sizeof(kIndexMagic) + 5 * sizeof(uint32_t) + sizeof(uint32_t) /* reserved */ +
             4 * sizeof(double) + ntiers * sizeof(TierEntry);
    }

    // Where one tier's six arrays sit, relative to the start of the file. The
    // writer and the reader both derive them from the tier's own offset, quad
    // count and the file's cell count, so the two cannot come to disagree about
    // the layout.
    struct TierLayout {
      uint64_t ratio, d1, ra, dec, cell_start, items, end;
    };

    TierLayout tier_layout(uint64_t base, uint64_t nquads, uint64_t ncells) {
      TierLayout l;
      l.ratio = base;
      l.d1 = align_up(l.ratio + nquads * 5 * sizeof(float));
      l.ra = align_up(l.d1 + nquads * sizeof(double));
      l.dec = align_up(l.ra + nquads * sizeof(double));
      l.cell_start = align_up(l.dec + nquads * sizeof(double));
      // One past the last cell, so a query can read cell_start[c + 1].
      l.items = align_up(l.cell_start + (ncells + 1) * sizeof(uint32_t));
      l.end = align_up(l.items + nquads * sizeof(uint32_t));
      return l;
    }

    uint64_t cells_for(uint32_t nbins) {
      return static_cast<uint64_t>(nbins) * nbins * nbins;
    }

    // Writes zeros up to `target`, keeping `pos` in step with the stream.
    void pad_to(std::ostream &f, uint64_t &pos, uint64_t target) {
      static const char zeros[kAlign] = {};
      while (pos < target) {
        const uint64_t n = std::min<uint64_t>(target - pos, kAlign);
        f.write(zeros, static_cast<std::streamsize>(n));
        pos += n;
      }
    }

    bool fail(std::string *error, const std::string &why) {
      if (error) *error = why;
      return false;
    }
  } // namespace

  bool save_index_file(const std::string &path, const std::vector<QuadIndex> &tiers) {
    if (tiers.empty()) return false;
    std::ofstream f(path, std::ios::binary | std::ios::trunc);
    if (!f) return false;

    const QuadIndexSettings &s0 = tiers.front().settings();
    // The bin count follows from the tolerance, which the whole file shares, so
    // every tier's grid is the same shape.
    const uint32_t nbins = static_cast<uint32_t>(tiers.front().nbins_);
    const uint64_t ncells = cells_for(nbins);
    if (nbins == 0) return false; // never finalised, so there is no grid to write

    std::vector<TierEntry> entries(tiers.size());
    uint64_t offset = align_up(header_bytes(tiers.size()));
    for (size_t k = 0; k < tiers.size(); k++) {
      const uint64_t n = tiers[k].size();
      entries[k] = {tiers[k].settings().star_density, n, offset};
      offset = tier_layout(offset, n, ncells).end;
    }

    f.write(kIndexMagic, sizeof(kIndexMagic));
    put(f, kIndexVersion);
    put(f, kByteOrderMark);
    put(f, static_cast<uint32_t>(tiers.front().settings().radius_deg >= 180 ? 0 : 1));
    put(f, static_cast<uint32_t>(tiers.size()));
    put(f, nbins);
    put(f, static_cast<uint32_t>(0)); // reserved, and keeps the doubles aligned
    put(f, s0.quad_tolerance);
    put(f, s0.centre_ra);
    put(f, s0.centre_dec);
    put(f, s0.radius_deg);
    for (const TierEntry &e : entries) put(f, e);

    uint64_t pos = header_bytes(tiers.size());
    for (size_t k = 0; k < tiers.size(); k++) {
      const QuadIndex &ix = tiers[k];
      const uint64_t n = entries[k].nquads;
      if (static_cast<uint32_t>(ix.nbins_) != nbins) return false; // mismatched grids
      const TierLayout l = tier_layout(entries[k].offset, n, ncells);
      pad_to(f, pos, l.ratio);
      put_view(f, ix.v_ratio_, n * 5);
      pos += n * 5 * sizeof(float);
      pad_to(f, pos, l.d1);
      put_view(f, ix.v_d1_, n);
      pos += n * sizeof(double);
      pad_to(f, pos, l.ra);
      put_view(f, ix.v_ra_, n);
      pos += n * sizeof(double);
      pad_to(f, pos, l.dec);
      put_view(f, ix.v_dec_, n);
      pos += n * sizeof(double);
      pad_to(f, pos, l.cell_start);
      put_view(f, ix.v_cell_start_, ncells + 1);
      pos += (ncells + 1) * sizeof(uint32_t);
      pad_to(f, pos, l.items);
      put_view(f, ix.v_items_, n);
      pos += n * sizeof(uint32_t);
      pad_to(f, pos, l.end); // so the next tier, and the file, end aligned too
    }
    f.flush();
    return static_cast<bool>(f);
  }

  namespace {
    // Shared by the header reader and the loader.
    bool open_and_read_header(const std::string &path, std::ifstream &f, QuadIndexFile &info,
                              std::vector<TierEntry> &entries, std::string *error) {
      f.open(path, std::ios::binary);
      if (!f) return fail(error, "cannot open " + path);

      char magic[8];
      f.read(magic, sizeof(magic));
      if (!f || std::memcmp(magic, kIndexMagic, sizeof(magic)) != 0)
        return fail(error, path + " is not a quad index file");

      uint32_t version = 0, bom = 0, reserved = 0, ntiers = 0, nbins = 0, spare = 0;
      get(f, version);
      get(f, bom);
      get(f, reserved);
      get(f, ntiers);
      if (version != kIndexVersion)
        return fail(error, path + " was written by a different version, rebuild it");
      if (bom != kByteOrderMark) return fail(error, path + " has the opposite byte order");
      if (ntiers == 0 || ntiers > 4096) return fail(error, path + " has an implausible tier count");
      get(f, nbins);
      get(f, spare);
      // The grid is indexed by cell number, so an implausible bin count would
      // send a query reading somewhere the tier does not reach. 4096 bins is a
      // tolerance of 0.00024, far finer than anything that matches.
      if (nbins < 2 || nbins > 4096) return fail(error, path + " has an implausible bin count");

      info.version = version;
      info.nbins = nbins;
      info.database_type = static_cast<int>(reserved);
      get(f, info.quad_tolerance);
      get(f, info.centre_ra);
      get(f, info.centre_dec);
      get(f, info.radius_deg);

      entries.resize(ntiers);
      info.densities.clear();
      info.quads.clear();
      info.bytes = 0;
      for (uint32_t k = 0; k < ntiers; k++) {
        get(f, entries[k]);
        info.densities.push_back(entries[k].density);
        info.quads.push_back(entries[k].nquads);
        info.bytes += entries[k].nquads * (5 * sizeof(float) + 3 * sizeof(double) + sizeof(uint32_t)) +
                      (cells_for(nbins) + 1) * sizeof(uint32_t);
      }
      if (!f) return fail(error, path + " is truncated");
      return true;
    }
  } // namespace

  bool read_index_file_header(const std::string &path, QuadIndexFile &info, std::string *error) {
    std::ifstream f;
    std::vector<TierEntry> entries;
    return open_and_read_header(path, f, info, entries, error);
  }

  bool load_index_file(const std::string &path, std::vector<QuadIndex> &out, double min_density,
                       double max_density, std::string *error) {
    out.clear();
    std::ifstream f;
    QuadIndexFile info;
    std::vector<TierEntry> entries;
    if (!open_and_read_header(path, f, info, entries, error)) return false;
    f.close();

    // One mapping for the whole file, shared by every tier taken out of it. The
    // quad arrays are then used where they lie: nothing is copied into this
    // process, and a page no query ever reaches is never read at all.
    //
    // Deliberately without advise_random(). MADV_RANDOM switches off readahead,
    // and readahead is what turns a cold cache into one request instead of thirty-two.
    // The image that finds no solution is the case that decides this: it sweeps the whole ladder
    const auto map = std::make_shared<MappedFile>();
    if (!map->open(path, error)) return false;
    const uint8_t *base = map->data();
    const uint64_t ncells = cells_for(info.nbins);

    for (const TierEntry &e : entries) {
      if (min_density > 0 && e.density < min_density) continue;
      if (max_density > 0 && e.density > max_density) continue;

      // The offsets come out of the file, so they are checked rather than
      // trusted: a wrong one would otherwise be a pointer into nothing, and a
      // misaligned one a double read the hardware may refuse.
      if (e.offset % kAlign != 0) return fail(error, path + " has a misaligned tier");
      const TierLayout l = tier_layout(e.offset, e.nquads, ncells);
      if (l.end > map->size()) return fail(error, path + " is truncated");

      QuadIndex ix;
      ix.settings_.star_density = e.density;
      ix.settings_.quad_tolerance = info.quad_tolerance;
      ix.settings_.centre_ra = info.centre_ra;
      ix.settings_.centre_dec = info.centre_dec;
      ix.settings_.radius_deg = info.radius_deg;
      ix.map_ = map;
      ix.nbins_ = static_cast<int>(info.nbins);
      ix.nquads_ = static_cast<size_t>(e.nquads);
      ix.nitems_ = static_cast<size_t>(e.nquads);
      ix.v_ratio_ = reinterpret_cast<const float *>(base + l.ratio);
      ix.v_d1_ = reinterpret_cast<const double *>(base + l.d1);
      ix.v_ra_ = reinterpret_cast<const double *>(base + l.ra);
      ix.v_dec_ = reinterpret_cast<const double *>(base + l.dec);
      // The grid comes out of the file too, which is the whole point: rebuilding
      // it would mean reading every quad, and then nothing would be lazy.
      ix.v_cell_start_ = reinterpret_cast<const uint32_t *>(base + l.cell_start);
      ix.v_items_ = reinterpret_cast<const uint32_t *>(base + l.items);
      out.push_back(std::move(ix));
    }
    if (out.empty()) return fail(error, path + " holds no tier in the requested density range");
    return true;
  }

  std::string default_index_cache_path(const std::string &db_name, int database_type,
                                       double quad_tolerance) {
    const char *xdg = std::getenv("XDG_CACHE_HOME");
#ifdef _WIN32
    // The index is a large file that is expensive to rebuild, so it belongs in
    // LocalAppData rather than in TEMP, which Windows is free to clear. This is
    // checked before HOME because a Git Bash or MSYS shell sets HOME as well,
    // and the cache should not move depending on which shell started the solve.
    const char *local = std::getenv("LOCALAPPDATA");
#else
    const char *local = nullptr;
#endif
    const char *home = std::getenv("HOME");
    std::filesystem::path dir;
    if (xdg && *xdg) {
      dir = std::filesystem::path(xdg) / "faster-astap";
    } else if (local && *local) {
      dir = std::filesystem::path(local) / "faster-astap" / "cache";
    } else if (home && *home) {
      dir = std::filesystem::path(home) / ".cache" / "faster-astap";
    } else {
      dir = std::filesystem::temp_directory_path() / "faster-astap";
    }

    // The tolerance is part of the name because it sets the bin width: an index
    // built at one tolerance cannot serve a solve at another.
    char name[160];
    std::snprintf(name, sizeof(name), "%s_%d_t%.4f.qix",
                  db_name.empty() ? "unknown" : db_name.c_str(), database_type, quad_tolerance);
    return (dir / name).string();
  }

  std::string index_tier_cache_path(const std::string &db_name, int database_type,
                                    double quad_tolerance, double density) {
    // Same directory and the same naming rules, with the rung appended. %g keeps
    // the sparse rungs readable (d0.5) and the deep ones exact (d3600).
    std::string path = default_index_cache_path(db_name, database_type, quad_tolerance);
    path.erase(path.size() - 4); // ".qix"
    char suffix[48];
    std::snprintf(suffix, sizeof(suffix), "_d%g.qix", density);
    return path + suffix;
  }

  bool ensure_parent_directory(const std::string &path) {
    const std::filesystem::path dir = std::filesystem::path(path).parent_path();
    if (dir.empty()) return true;
    std::error_code ec;
    if (std::filesystem::exists(dir, ec)) return std::filesystem::is_directory(dir, ec);
    std::filesystem::create_directories(dir, ec);
    return !ec;
  }

  void QuadIndex::query(const float *r, std::vector<uint32_t> &hits) const {
    if (nitems_ == 0) return;
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
          for (uint32_t p = v_cell_start_[c]; p < v_cell_start_[c + 1]; p++) {
            const uint32_t q = v_items_[p];
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
