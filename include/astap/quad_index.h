// A searchable index of database quads, built once and queried per image.
// See the index solver sections of README.md for the design and its risks.

#pragma once

#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "astap/mapped_file.h"
#include "astap/star_database.h"
#include "astap/types.h"

namespace astap {
  struct QuadIndexSettings {
    // Stars per square degree to take from the database. This is the depth
    // tier: it has to match the depth of the images being solved, because a
    // quad is only findable when all four of its stars are bright enough to be
    // detected. See the "depth tiers" section of the design note.
    double star_density = 300.0;

    // Must match the solver's matching tolerance; it also sets the hash bin
    // width, so an index is only valid for the tolerance it was built with.
    double quad_tolerance = 0.007;

    // Below this density, build every quad from each star's six nearest
    // neighbours (C(6,4) = 15 per star) instead of just its three nearest.
    //
    // A quad only matches when the image and the catalogue picked the same four
    // stars, and at a few stars per square degree that agreement is fragile: one
    // detection the catalogue subset does not contain changes which three
    // neighbours a star has, and so replaces its quad entirely. A ten degree
    // field yields only ~76 quads, with no redundancy to absorb that. The larger
    // group is a strict superset - C(6,4) contains the three-nearest quad - so it
    // can only add matches, and the sparse tiers are small enough to afford it.
    double many_quads_below_density = 5.0;

    // Restrict the build to a cap around this position. radius_deg >= 180
    // covers the sky.
    double centre_ra = 0; // radians
    double centre_dec = 0; // radians
    double radius_deg = 180;
  };

  // Quads are stored as parallel arrays: hash ratios in fp32 (the drift against
  // the tolerance is 3.2e-7, a 21800x margin), scale and position in fp64.
  class QuadIndex {
  public:
    QuadIndex() = default;

    // The quad arrays are addressed through pointers, which a copy or a move of
    // the owning vectors invalidates, so both have to re-point them. Tiers are
    // copied into and moved around vectors on every load, so this is not a
    // theoretical case.
    QuadIndex(const QuadIndex &o) { *this = o; }
    QuadIndex(QuadIndex &&o) noexcept { *this = std::move(o); }

    QuadIndex &operator=(const QuadIndex &o);

    QuadIndex &operator=(QuadIndex &&o) noexcept;

    // Builds the index by walking the database tiles. Returns false when the
    // database cannot be read. `progress` is called with a 0..1 fraction.
    bool build(StarDatabase &db, const QuadIndexSettings &s,
               const std::function<void(double)> &progress = nullptr);

    size_t size() const { return nquads_; }
    const QuadIndexSettings &settings() const { return settings_; }

    // Ratios 1..5 of one quad.
    const float *ratios(size_t i) const { return v_ratio_ + i * 5; }
    double d1(size_t i) const { return v_d1_[i]; } // longest side, arcsec
    double centre_ra(size_t i) const { return v_ra_[i]; } // quad centre, radians
    double centre_dec(size_t i) const { return v_dec_[i]; }

    // True when the quads are read straight out of a mapped cache file rather
    // than held in this process's own memory.
    bool mapped() const { return map_ != nullptr; }

    // Appends the indices of every quad whose five ratios lie within
    // `quad_tolerance` of the given ones.
    void query(const float *r, std::vector<uint32_t> &hits) const;

    // Rough memory footprint in bytes.
    size_t bytes() const;

  private:
    friend bool build_tiers(StarDatabase &, const QuadIndexSettings &, const std::vector<double> &,
                            std::vector<QuadIndex> &, const std::function<void(double)> &);

    friend bool save_index_file(const std::string &, const std::vector<QuadIndex> &);

    friend bool load_index_file(const std::string &, std::vector<QuadIndex> &, double, double,
                                std::string *);

    // Sorts the accumulated quads into the 3D bin grid. Called once the quad
    // arrays are complete, by either build path.
    void finalise();

    // Aims the views below at whatever storage this instance holds. Call after
    // the owning vectors are filled, or after they move.
    void repoint();

    // Bin on the first three ratios; probing the 3x3x3 neighbourhood covers the
    // tolerance ball in those dimensions, the other two are checked exactly.
    int bin_of(float v) const;

    uint32_t cell_of(int b0, int b1, int b2) const;

    QuadIndexSettings settings_;

    // The quads sit in one of two places: these vectors, when the index was
    // built in this process, or a mapped cache file shared with the other tiers
    // that came out of it. Everything that reads a quad goes through the views,
    // so the two cases are distinguished here and nowhere else.
    std::shared_ptr<const MappedFile> map_;
    std::vector<float> ratio_; // 5 per quad
    std::vector<double> d1_;
    std::vector<double> ra_, dec_;
    std::vector<uint32_t> cell_start_; // CSR over the 3D bin grid
    std::vector<uint32_t> items_;

    const float *v_ratio_ = nullptr;
    const double *v_d1_ = nullptr;
    const double *v_ra_ = nullptr;
    const double *v_dec_ = nullptr;
    const uint32_t *v_cell_start_ = nullptr;
    const uint32_t *v_items_ = nullptr;
    size_t nquads_ = 0;
    size_t nitems_ = 0;
    int nbins_ = 0;
  };

  // Builds one index per requested depth tier in a single pass over the database.
  //
  // A ladder of tiers is what makes the solver as capable as the spiral search,
  // which adapts its depth at every position. Measured over the corpus, every
  // image the index solver missed was missed for this reason alone: a field at
  // 1000 stars/deg^2 or at 5 stars/deg^2 has no findable quads in an index built
  // at 300, because a quad is only findable when all four of its stars were
  // bright enough to be detected in the image.
  //
  // Building them together rather than one at a time matters: each tile's stars
  // are read and projected once, at the deepest density requested, and every
  // tier takes a prefix of that list. The read is the expensive part.
  //
  // `base` supplies the tolerance and the sky cap; `base.star_density` is
  // ignored. `out` comes back in the same order as `densities`.
  bool build_tiers(StarDatabase &db, const QuadIndexSettings &base,
                   const std::vector<double> &densities, std::vector<QuadIndex> &out,
                   const std::function<void(double)> &progress = nullptr);

  // --- on-disk index ---------------------------------------------------------
  //
  // A ladder takes seconds to build and is identical for every image and every
  // run, so a released solver builds it once and reads it back. Only the quad
  // arrays are stored; the bin grid is rebuilt on load, which is faster than
  // reading it.
  //
  // The file records the tolerance it was binned with and the database it came
  // from. Loading refuses a file written by a different build, so a stale cache
  // fails loudly instead of producing quietly wrong matches.
  //
  // Every array starts on a 64 byte boundary, padded apart rather than packed
  // end to end. That costs at most 63 bytes per array on a file of gigabytes and
  // is what lets the file be memory mapped and the arrays used where they lie:
  // packed, a tier holding an odd number of quads puts its d1/ra/dec arrays at a
  // 4 byte offset, which a copying read does not care about and a `const double
  // *` into a mapping does.

  struct QuadIndexFile {
    uint32_t version = 0;
    uint32_t nbins = 0; // bins per ratio axis; the grid holds nbins^3 cells
    int database_type = 0;
    double quad_tolerance = 0;
    double centre_ra = 0, centre_dec = 0, radius_deg = 180;
    std::vector<double> densities; // one per tier, in file order
    std::vector<uint64_t> quads; // quads per tier
    uint64_t bytes = 0;
  };

  bool save_index_file(const std::string &path, const std::vector<QuadIndex> &tiers);

  // Reads the tiers whose density lies in [min_density, max_density]; 0 and 0
  // load every tier. Loading a subset is what lets a caller who knows the field
  // size pay for one tier instead of the ladder. `error` receives a reason on
  // failure.
  bool load_index_file(const std::string &path, std::vector<QuadIndex> &out,
                       double min_density = 0, double max_density = 0,
                       std::string *error = nullptr);

  // Reads the header only: which tiers a file holds, without loading them.
  bool read_index_file_header(const std::string &path, QuadIndexFile &info,
                              std::string *error = nullptr);

  // Where a built ladder is cached by default: `$XDG_CACHE_HOME/faster-astap`,
  // or `$HOME/.cache/faster-astap`. One file per star database, named for the
  // database and the tolerance, because an index is only reusable by a solve
  // that would have built the same one. Several databases can therefore be
  // cached side by side.
  std::string default_index_cache_path(const std::string &db_name, int database_type,
                                       double quad_tolerance);

  // Where one rung is cached, in the same directory and with the density in the
  // name. Tiers are cached one per file so that ladders compose: a run that
  // wants a deeper ceiling builds and stores only the rungs it adds, and the
  // ones it shares with the default ladder are read back rather than rebuilt.
  // A single file for the whole ladder cannot do that - it belongs to exactly
  // one ladder, so asking for one more rung discards all the others.
  std::string index_tier_cache_path(const std::string &db_name, int database_type,
                                    double quad_tolerance, double density);

  // Creates the directory holding `path` when it does not exist.
  bool ensure_parent_directory(const std::string &path);
} // namespace astap
