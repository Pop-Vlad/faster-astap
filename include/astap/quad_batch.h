// Batched quad construction.
//
// The spiral search evaluates a batch of sky positions at a time. Building the
// database quads for a whole batch through one call lets the thread pool spread
// the positions over the cores: per position the work is only a few hundred
// stars, far too little to parallelise on its own.
//
// Produces exactly what calling find_quads() on each list separately produces,
// including the in-place sort of the star list.

#pragma once

#include <vector>

#include "astap/types.h"

namespace astap {
  // For every i, equivalent to find_quads(nrstars_image, lists[i], out[i]).
  // `lists` is modified in place exactly as find_quads would modify it.
  void build_quads_batch(int nrstars_image, std::vector<RowList> &lists, std::vector<RowList> &out);
} // namespace astap
