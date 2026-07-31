// Solver built on a pre-built quad index instead of a spiral search.
// See the index solver sections of README.md.
//
// The matching front end is new — one pass over the image quads against a
// whole-sky index, a joint vote on scale and position, then a RANSAC consensus
// — but the back end is the port's: the same overdetermined system solved by
// `lsq_fit`, the same acceptance checks, so a solution has to clear the same bar
// the reference implementation sets.

#pragma once

#include <string>
#include <vector>

#include "astap/quad_index.h"
#include "astap/types.h"

namespace astap {
  struct IndexSolveSettings {
    // The matching tolerance is deliberately absent: it belongs to the
    // QuadIndex, which is binned by it and is only valid for the value it was
    // built with.

    // Minimum surviving quads. The port uses 3 + nrstars_image/140; 3 gives the
    // three quad centres that are the bare minimum for the fit.
    int minimum_quads = 3;
    // Votes are cast into (log scale, sky cell) buckets. This is the log-space
    // bin width, so 0.01 is one percent of plate scale.
    double scale_bin = 0.01;
    // Sky cell size for the position vote, in degrees. 0 derives it from the
    // image: a cell has to be at least as large as the field, or a field spread
    // across several cells splits its own vote.
    double sky_cell_deg = 0;
    // How many vote peaks to verify. A peak is cheap to reject, so this is
    // generous; the joint vote makes the true peak rank highly but not always
    // first.
    int top_k = 40;

    // --- RANSAC verification -------------------------------------------------
    // A pair is an inlier when the model puts its image quad centre within this
    // many pixels of its database quad centre.
    double inlier_tolerance_px = 4.0;
    // Consensus needed before a candidate position is fitted at all. Two pairs
    // define the model, so this is the number of independent confirmations plus
    // two.
    int min_inliers = 4;
    // Minimal samples drawn per candidate. Small clusters are enumerated
    // exhaustively and never reach this.
    int max_samples = 3000;
  };

  struct IndexSolveResult {
    bool solved = false;
    std::string reason;  // why it failed, when it did

    double ra0 = 0, dec0 = 0;      // image centre, radians
    double crpix1 = 0, crpix2 = 0;
    double cdelt1 = 0, cdelt2 = 0; // deg/px
    double crota1 = 0, crota2 = 0; // deg
    double cd1_1 = 0, cd1_2 = 0, cd2_1 = 0, cd2_2 = 0;

    int nr_matches = 0;     // candidate pairs from the index
    int nr_references = 0;  // pairs used in the final fit
    int nr_inliers = 0;     // consensus size the fit was built from
    int peaks_tried = 0;    // vote peaks verified before one passed
    double scale_arcsec_px = 0;
    double tier_density = 0;  // depth tier that produced the solution
    int tiers_tried = 0;
    // Set when the solution needed the larger quad groups of the second pass.
    bool many_quads_pass = false;
    // The tolerance the matching index was built with, carried so a later stage
    // can match on the same terms.
    double quad_tolerance_used = 0.007;
    // Set once refine_with_database has improved this solution.
    bool refined = false;
  };

  // `image_quads` comes from find_quads on the detected stars, in image pixel
  // coordinates. `width`/`height` are the dimensions those coordinates refer to.
  IndexSolveResult solve_with_index(const QuadIndex &index, const RowList &image_quads, int width,
                                    int height, const IndexSolveSettings &s = {});

  // Solves against a ladder of depth tiers (see build_tiers). One index commits
  // to one depth, and a quad is only findable when all four of its stars were
  // bright enough to be detected, so a single tier can only solve images within
  // about a factor of two of its own density. The spiral search adapts its depth
  // at every position; sweeping the ladder is how the index solver earns the
  // same range back.
  //
  // `density_hint` is the image's detected stars per square degree when the
  // caller knows the field size; it only reorders the sweep, so a wrong hint
  // costs time rather than a solution. A tier that reaches a strong consensus
  // ends the sweep; otherwise the best consensus across all tiers wins.
  IndexSolveResult solve_with_tiers(const std::vector<QuadIndex> &tiers, const RowList &image_quads,
                                    int width, int height, const IndexSolveSettings &s = {},
                                    double density_hint = 0);

  // Solves from a detected star list rather than from pre-built quads.
  //
  // Quads are first built one per star from its three nearest neighbours, which
  // solves most images. If that finds nothing, they are rebuilt from every
  // combination of each star's six nearest — fifteen per star, a strict superset
  // — which is what sparse and wide fields need: a quad matches only when the
  // image and the catalogue chose the same four stars, and at a few stars per
  // square degree a single differing detection replaces a star's quad outright.
  //
  // Prefer this over solve_with_tiers when the star list is available.
  IndexSolveResult solve_stars_with_tiers(const std::vector<QuadIndex> &tiers, const RowList &stars,
                                          int width, int height, const IndexSolveSettings &s = {},
                                          double density_hint = 0);

  struct IndexRefineResult {
    bool ok = false;
    int nr_quads = 0;      // quads matched in the second pass
    int nr_candidates = 0; // database quads it had to match against
    bool sip_valid = false;
    bool kept = false;  // false when the refit was discarded as less well supported
    SipCoefficients sip;
    std::string reason;
  };

  // Second pass, once the field is known.
  //
  // The index gives a position from quad centres alone, which is enough to
  // identify the field but leaves accuracy on the table and rarely produces the
  // twenty matched quads a cubic distortion fit needs. This reads the database
  // once at the solved position, at a depth matched to the image's own star
  // count, and redoes the match there — the same work the spiral does at each of
  // its many positions, done once at the right one. It costs 0.3 to 0.9 ms all
  // in, against a 5 ms solve.
  //
  // On success `io` is updated with the improved linear solution. SIP is fitted
  // from the matched quad centres, exactly as the port does, and needs at least
  // twenty of them.
  //
  // The refined solution is kept only when it rests on at least as many quads as
  // the consensus it would replace. Measured over the corpus, every case where
  // refining made the position worse was one where it matched fewer quads than
  // the index already had; a fit from more independent points is the better
  // constrained one, and a weakly supported fit should not displace a strong one.
  //
  // `stars` and `io` must be in the same pixel frame, and the returned SIP
  // coefficients are in that frame too; a caller working on a binned copy has to
  // rescale them (see scale_sip_to_original in the index solver front end).
  IndexRefineResult refine_with_database(StarDatabase &db, const RowList &stars, int width,
                                         int height, IndexSolveResult &io,
                                         const IndexSolveSettings &s = {}, bool want_sip = false);
} // namespace astap
