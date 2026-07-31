// Batched quad construction.
//
// The spiral search evaluates a batch of sky positions at a time. Building the
// database quads for a whole batch through one interface gives an accelerator
// something worth launching for: per position the work is only a few hundred
// stars, far too little to be worth a kernel launch on its own.
//
// Every implementation must produce exactly what calling find_quads() on each
// list separately produces, including the in-place sort of the star list.

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "astap/types.h"

namespace astap {
  class QuadBatchBuilder {
  public:
    virtual ~QuadBatchBuilder() = default;

    // For every i, equivalent to find_quads(nrstars_image, lists[i], out[i]).
    // `lists` is modified in place exactly as find_quads would modify it.
    virtual void build(int nrstars_image, std::vector<RowList> &lists,
                       std::vector<RowList> &out) = 0;

    // Short name for logging, e.g. "cpu" or "sycl:NVIDIA GeForce RTX 3090".
    virtual std::string backend() const = 0;

    // How many positions this builder would like per call. 0 means it has no
    // preference and the caller's default applies. An accelerator wants a large
    // batch so the kernels dominate the launch and transfer cost; the CPU does
    // not care, and a smaller batch keeps the search's early exit cheap.
    virtual size_t preferred_batch() const { return 0; }
  };

  // Always available. Runs find_quads over the batch on the thread pool.
  std::unique_ptr<QuadBatchBuilder> make_cpu_quad_builder();

  // Returns nullptr when the build has no SYCL support or no suitable device
  // was found, so callers can fall back to the CPU builder.
  std::unique_ptr<QuadBatchBuilder> make_sycl_quad_builder();

  // True when the binary was built with SYCL support at all.
  bool sycl_available();
} // namespace astap
