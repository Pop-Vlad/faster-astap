#include "astap/quad_batch.h"

#include "astap/parallel.h"
#include "astap/quads.h"

namespace astap {
  void build_quads_batch(int nrstars_image, std::vector<RowList> &lists, std::vector<RowList> &out) {
    out.resize(lists.size());
    parallel_for(0, lists.size(), [&](size_t i, unsigned) {
      find_quads(nrstars_image, lists[i], out[i]);
    });
  }
} // namespace astap
