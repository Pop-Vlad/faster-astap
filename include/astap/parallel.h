// Minimal threading helpers. No dependencies beyond the standard library.
//
// Every parallel region in this port is written so that the result is bit
// identical to the sequential version: work is split into disjoint index
// ranges, reductions are either integer (order independent) or combined in a
// fixed thread order, and no floating point accumulation order is changed.

#pragma once

#include <algorithm>
#include <cstddef>
#include <functional>
#include <thread>
#include <vector>

namespace astap {
  // Number of worker threads to use. 0 (the default) means one per hardware
  // thread the process is allowed to run on: all logical processors of the
  // affinity mask, SMT siblings included. Set to 1 to force sequential execution.
  void set_thread_count(unsigned n);

  unsigned thread_count();

  // Splits [begin, end) into one contiguous chunk per thread and calls
  // body(chunk_begin, chunk_end, thread_index). Chunks are always contiguous and
  // assigned in increasing order, so a caller may rely on chunk k covering lower
  // indices than chunk k+1.
  void parallel_ranges(size_t begin, size_t end,
                       const std::function<void(size_t, size_t, unsigned)> &body);

  // Convenience wrapper: calls body(i, thread_index) for every i in [begin, end).
  inline void parallel_for(size_t begin, size_t end,
                           const std::function<void(size_t, unsigned)> &body) {
    parallel_ranges(begin, end, [&](size_t b, size_t e, unsigned t) {
      for (size_t i = b; i < e; i++) body(i, t);
    });
  }

  // Number of chunks parallel_ranges would create for a range of `count` items.
  unsigned range_chunks(size_t count);
} // namespace astap
