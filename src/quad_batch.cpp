#include "astap/quad_batch.h"

#include "astap/parallel.h"
#include "astap/quads.h"

namespace astap {
  namespace {
    class CpuQuadBuilder : public QuadBatchBuilder {
    public:
      void build(int nrstars_image, std::vector<RowList> &lists, std::vector<RowList> &out) override {
        out.resize(lists.size());
        parallel_for(0, lists.size(), [&](size_t i, unsigned) {
          find_quads(nrstars_image, lists[i], out[i]);
        });
      }

      std::string backend() const override { return "cpu"; }
    };
  } // namespace

  std::unique_ptr<QuadBatchBuilder> make_cpu_quad_builder() {
    return std::unique_ptr<QuadBatchBuilder>(new CpuQuadBuilder());
  }

#ifndef ASTAP_WITH_SYCL
  std::unique_ptr<QuadBatchBuilder> make_sycl_quad_builder() { return nullptr; }
  bool sycl_available() { return false; }
#endif
} // namespace astap
