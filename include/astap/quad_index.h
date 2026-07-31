// A searchable index of database quads, built once and queried per image.
// See docs/index_solver.md for the design and its risks.

#pragma once

#include <cstdint>
#include <functional>
#include <string>
#include <vector>

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

    // Restrict the build to a cap around this position. radius_deg >= 180
    // covers the sky.
    double centre_ra = 0;   // radians
    double centre_dec = 0;  // radians
    double radius_deg = 180;
  };

  // Quads are stored as parallel arrays: hash ratios in fp32 (the drift against
  // the tolerance is 3.2e-7, a 21800x margin), scale and position in fp64.
  class QuadIndex {
  public:
    // Builds the index by walking the database tiles. Returns false when the
    // database cannot be read. `progress` is called with a 0..1 fraction.
    bool build(StarDatabase &db, const QuadIndexSettings &s,
               const std::function<void(double)> &progress = nullptr);

    size_t size() const { return d1_.size(); }
    const QuadIndexSettings &settings() const { return settings_; }

    // Ratios 1..5 of one quad.
    const float *ratios(size_t i) const { return &ratio_[i * 5]; }
    double d1(size_t i) const { return d1_[i]; }        // longest side, arcsec
    double centre_ra(size_t i) const { return ra_[i]; }  // quad centre, radians
    double centre_dec(size_t i) const { return dec_[i]; }

    // Appends the indices of every quad whose five ratios lie within
    // `quad_tolerance` of the given ones.
    void query(const float *r, std::vector<uint32_t> &hits) const;

    // Rough memory footprint in bytes.
    size_t bytes() const;

  private:
    // Bin on the first three ratios; probing the 3x3x3 neighbourhood covers the
    // tolerance ball in those dimensions, the other two are checked exactly.
    int bin_of(float v) const;
    uint32_t cell_of(int b0, int b1, int b2) const;

    QuadIndexSettings settings_;
    std::vector<float> ratio_;   // 5 per quad
    std::vector<double> d1_;
    std::vector<double> ra_, dec_;

    int nbins_ = 0;
    std::vector<uint32_t> cell_start_;  // CSR over the 3D bin grid
    std::vector<uint32_t> items_;
  };
} // namespace astap
