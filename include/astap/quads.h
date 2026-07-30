// Step 3 to 5 of the ASTAP method: build the smallest irregular tetrahedrons
// (quads) out of four close stars, calculate the six distances between the four
// stars and the mean x,y position of the quad, sort the distances and scale them
// to the longest one. The result is the quad hash code.
//
// Ported from unit_command_line_solving.pas (find_quads, find_many_quads).

#pragma once

#include "astap/types.h"

namespace astap {
  // Sorts the star list in X only. The SNR row is *not* updated, exactly like the
  // original: after this call row 2 no longer belongs to the star in rows 0 and 1.
  void quicksort_starlist(RowList &a, long lo, long hi);

  // Builds quads from groups of the `mode` (5, 6 or 7) closest stars, giving
  // 5, 15 respectively 35 quads per star. Used when the image has few stars.
  void find_many_quads(const RowList &starlist, RowList &quads, int mode);

  // Builds quads from the three closest stars of every star.
  // `nrstars_image` is the star count of the *image*, also when building the
  // database quads: the group size has to be based on the image, not on the
  // possibly larger database field.
  // Output rows: 0 = longest distance d1, 1..5 = d2/d1 .. d6/d1,
  // 6 = mean x of the quad, 7 = mean y of the quad.
  void find_quads(int nrstars_image, RowList &starlist, RowList &quads);
} // namespace astap
