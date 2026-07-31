// Solver built on a pre-built quad index instead of a spiral search.
// See docs/index_solver.md.
//
// The matching front end is new — one pass over the image quads against a
// whole-sky index, then a vote on scale and position — but the back end is the
// port's: the same overdetermined system solved by `lsq_fit`, the same
// acceptance checks, so a solution has to clear the same bar the reference
// implementation sets.

#pragma once

#include <string>
#include <vector>

#include "astap/quad_index.h"
#include "astap/types.h"

namespace astap {
  struct IndexSolveSettings {
    double quad_tolerance = 0.007;
    // Minimum surviving quads. The port uses 3 + nrstars_image/140; 3 gives the
    // three quad centres that are the bare minimum for the fit.
    int minimum_quads = 3;
    // Scale votes are histogrammed in log space with this bin width, so 0.01 is
    // one percent of plate scale.
    double scale_bin = 0.01;
    // Sky cells for the position vote, in degrees.
    double sky_cell_deg = 1.0;
    // How many position peaks to try before giving up. The best cell is not
    // always the right one when the vote is thin.
    int top_k = 5;
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
    double scale_arcsec_px = 0;
  };

  // `image_quads` comes from find_quads on the detected stars, in image pixel
  // coordinates. `width`/`height` are the dimensions those coordinates refer to.
  IndexSolveResult solve_with_index(const QuadIndex &index, const RowList &image_quads, int width,
                                    int height, const IndexSolveSettings &s = {});
} // namespace astap
