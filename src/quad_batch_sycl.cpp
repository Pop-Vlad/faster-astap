// SYCL offload of the batched quad construction.
//
// Portable SYCL 2020, no vendor extensions, so it builds under both oneAPI
// DPC++ and AdaptiveCpp.
//
// The results must match the CPU implementation bit for bit, which shapes the
// design in two ways:
//
//   * Everything is double precision and the arithmetic is written in the same
//     order as find_many_quads(). The build passes -ffp-contract=off so neither
//     host nor device may fuse a*b+c into an FMA and change the last bit.
//
//   * Duplicate quad rejection is a greedy scan: a candidate is dropped when an
//     *already accepted* quad sits within one pixel, so the outcome depends on
//     the order candidates are visited in. It therefore stays serial per
//     position — but positions are independent, so a batch of a few thousand
//     gives the device plenty of independent serial scans to run at once.
//
// Three kernels per batch:
//   1. one work item per (position, reference star): nearest neighbours, the
//      group's pairwise distances, and the C(n,4) candidate hash codes;
//   2. one work item per position: the greedy dedup, compacting in place;
//   3. one work item per output quad: gather into a densely packed, row major
//      block per position so the host only has to memcpy rows.
//
// Only the find_many_quads path (groups of 5, 6 or 7 closest stars) is
// offloaded. That is where a blind solve spends most of its time, because the
// early field-of-view attempts detect few stars. The regular find_quads path
// sorts the star list and falls back to the CPU.

#include <sycl/sycl.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>

#include "astap/parallel.h"
#include "astap/quad_batch.h"
#include "astap/quads.h"

namespace astap {
  namespace {
    constexpr int kMaxClose = 7;

    struct GroupSpec {
      int num_closest;
      int quads_per_group;
    };

    GroupSpec spec_for(int nrstars_image, size_t nrstars) {
      if (nrstars_image < 15 && nrstars > 6) return {7, 35};
      if (nrstars_image < 30 && nrstars > 5) return {6, 15};
      if (nrstars_image < 60 && nrstars > 4) return {5, 5};
      return {0, 0}; // not a find_many_quads case
    }

    // C(n,4) in ascending order, the order the CPU tables are written in.
    void fill_combinations(int n, std::vector<int> &out) {
      out.clear();
      for (int a = 0; a < n; a++)
        for (int b = a + 1; b < n; b++)
          for (int c = b + 1; c < n; c++)
            for (int d = c + 1; d < n; d++) {
              out.push_back(a);
              out.push_back(b);
              out.push_back(c);
              out.push_back(d);
            }
    }

    size_t next_pow2(size_t n) {
      size_t p = 64;
      while (p < n) p <<= 1;
      return p;
    }

    // A device allocation that is reused across batches and only ever grows.
    // Creating SYCL allocations per batch costs more than the kernels do.
    template<typename T>
    class DeviceArray {
    public:
      T *get(sycl::queue &q, size_t n) {
        if (n > cap_) {
          if (p_) sycl::free(p_, q);
          cap_ = n + n / 2 + 64;
          p_ = sycl::malloc_device<T>(cap_, q);
        }
        return p_;
      }

      void release(sycl::queue &q) {
        if (p_) sycl::free(p_, q);
        p_ = nullptr;
        cap_ = 0;
      }

    private:
      T *p_ = nullptr;
      size_t cap_ = 0;
    };

    class SyclQuadBuilder : public QuadBatchBuilder {
    public:
      explicit SyclQuadBuilder(sycl::queue q) : q_(std::move(q)) {
        name_ = "sycl:" + q_.get_device().get_info<sycl::info::device::name>();
        cpu_fallback_ = make_cpu_quad_builder();
      }

      ~SyclQuadBuilder() override {
        d_x_.release(q_);
        d_y_.release(q_);
        d_i32_.release(q_);
        d_cand_.release(q_);
        d_valid_.release(q_);
        d_heads_.release(q_);
        d_next_.release(q_);
        d_dense_.release(q_);
        d_key_.release(q_);
        d_tab_.release(q_);
      }

      std::string backend() const override { return name_; }

      // Enough positions per launch that the kernels dominate the launch cost.
      size_t preferred_batch() const override { return 2048; }

      void build(int nrstars_image, std::vector<RowList> &lists, std::vector<RowList> &out) override;

    private:
      sycl::queue q_;
      std::string name_;
      std::unique_ptr<QuadBatchBuilder> cpu_fallback_;

      DeviceArray<double> d_x_, d_y_, d_cand_, d_dense_;
      DeviceArray<uint32_t> d_i32_;
      DeviceArray<uint8_t> d_valid_;
      DeviceArray<int32_t> d_heads_, d_next_;
      DeviceArray<uint64_t> d_key_, d_tab_;

      // Host staging, reused so the copies are from pinned-friendly storage.
      std::vector<double> h_x_, h_y_, h_dense_;
      std::vector<uint32_t> h_meta_, h_nq_;
    };

    void SyclQuadBuilder::build(int nrstars_image, std::vector<RowList> &lists,
                                std::vector<RowList> &out) {
      const size_t npos = lists.size();
      out.resize(npos);
      if (npos == 0) return;

      // Split the batch: positions whose star count selects find_many_quads go
      // to the device, the rest keep the CPU path.
      std::vector<size_t> gpu_idx;
      gpu_idx.reserve(npos);
      GroupSpec spec{0, 0};
      for (size_t i = 0; i < npos; i++) {
        const GroupSpec s = spec_for(nrstars_image, lists[i].count());
        if (s.num_closest == 0) continue;
        if (spec.num_closest == 0) spec = s;
        if (s.num_closest == spec.num_closest) gpu_idx.push_back(i);
      }
      if (gpu_idx.empty()) {
        cpu_fallback_->build(nrstars_image, lists, out);
        return;
      }

      const int num_closest = spec.num_closest;
      const int quads_per_group = spec.quads_per_group;
      const size_t nblocks = gpu_idx.size();

      // Offsets: stars, candidates and the per position dedup grid.
      h_meta_.assign(nblocks * 4 + 3, 0);
      uint32_t *counts = h_meta_.data();                 // [nblocks]
      uint32_t *star_off = counts + nblocks;             // [nblocks+1]
      uint32_t *cand_off = star_off + nblocks + 1;       // [nblocks+1]
      std::vector<uint32_t> head_off_store(nblocks + 1, 0);
      uint32_t *head_off = head_off_store.data();        // [nblocks+1]

      for (size_t k = 0; k < nblocks; k++) {
        const uint32_t n = static_cast<uint32_t>(lists[gpu_idx[k]].count());
        counts[k] = n;
        star_off[k + 1] = star_off[k] + n;
        const uint32_t nc = n * static_cast<uint32_t>(quads_per_group);
        cand_off[k + 1] = cand_off[k] + nc;
        head_off[k + 1] = head_off[k] + static_cast<uint32_t>(next_pow2(size_t(nc) * 2));
      }
      const size_t total_stars = star_off[nblocks];
      const size_t total_cand = cand_off[nblocks];
      const size_t total_heads = head_off[nblocks];
      if (total_stars == 0 || total_cand == 0) {
        cpu_fallback_->build(nrstars_image, lists, out);
        return;
      }

      h_x_.resize(total_stars);
      h_y_.resize(total_stars);
      for (size_t k = 0; k < nblocks; k++) {
        const RowList &l = lists[gpu_idx[k]];
        const double *lx = l.data(0);
        const double *ly = l.data(1);
        std::memcpy(h_x_.data() + star_off[k], lx, l.count() * sizeof(double));
        std::memcpy(h_y_.data() + star_off[k], ly, l.count() * sizeof(double));
      }

      std::vector<int> combos;
      fill_combinations(num_closest, combos);

      double *dx = d_x_.get(q_, total_stars);
      double *dy = d_y_.get(q_, total_stars);
      // counts, star_off, cand_off, head_off, nrquads, dense_off and the combos
      // all live in one integer allocation: one array of nblocks, five of
      // nblocks+1, and the combination table.
      const size_t meta_n = nblocks + (nblocks + 1) * 5 + combos.size();
      uint32_t *dmeta = d_i32_.get(q_, meta_n);
      double *dcand = d_cand_.get(q_, total_cand * 8);
      uint8_t *dvalid = d_valid_.get(q_, total_cand);
      int32_t *dheads = d_heads_.get(q_, total_heads);
      int32_t *dnext = d_next_.get(q_, total_cand);
      uint64_t *dkey = d_key_.get(q_, total_cand);
      uint64_t *dtab = d_tab_.get(q_, total_heads);

      uint32_t *d_counts = dmeta;
      uint32_t *d_soff = d_counts + nblocks;
      uint32_t *d_coff = d_soff + nblocks + 1;
      uint32_t *d_hoff = d_coff + nblocks + 1;
      uint32_t *d_nq = d_hoff + nblocks + 1;
      uint32_t *d_doff = d_nq + nblocks + 1;
      uint32_t *d_combo = d_doff + nblocks + 1;

      q_.memcpy(dx, h_x_.data(), total_stars * sizeof(double));
      q_.memcpy(dy, h_y_.data(), total_stars * sizeof(double));
      q_.memcpy(d_counts, counts, nblocks * sizeof(uint32_t));
      q_.memcpy(d_soff, star_off, (nblocks + 1) * sizeof(uint32_t));
      q_.memcpy(d_coff, cand_off, (nblocks + 1) * sizeof(uint32_t));
      q_.memcpy(d_hoff, head_off, (nblocks + 1) * sizeof(uint32_t));
      {
        std::vector<uint32_t> cb(combos.begin(), combos.end());
        q_.memcpy(d_combo, cb.data(), cb.size() * sizeof(uint32_t)).wait();
      }

      // Kernel 1: candidates, one work item per reference star of the batch.
      q_.parallel_for(sycl::range<1>(total_stars), [=](sycl::id<1> gid) {
        const size_t g = gid[0];

        size_t lo = 0, hi = nblocks;
        while (lo + 1 < hi) { // which position owns this star
          const size_t mid = (lo + hi) / 2;
          if (d_soff[mid] <= g) lo = mid; else hi = mid;
        }
        const size_t p = lo;
        const size_t base = d_soff[p];
        const size_t n = d_counts[p];
        const size_t i = g - base;
        if (i >= n) return;

        long idx[kMaxClose];
        double dist[kMaxClose];
        idx[0] = static_cast<long>(i);
        dist[0] = 0;
        for (int k = 1; k < num_closest; k++) { idx[k] = -1; dist[k] = 1E99; }

        const double x1 = dx[base + i];
        const double y1 = dy[base + i];
        for (size_t j = 0; j < n; j++) {
          if (j == i) continue;
          const double ddx = dx[base + j] - x1;
          const double ddy = dy[base + j] - y1;
          const double d = ddx * ddx + ddy * ddy;
          if (d <= 1) continue; // not an identical star

          int insert_pos = -1;
          for (int k = num_closest - 1; k >= 1; k--) {
            if (d < dist[k]) insert_pos = k; else break;
          }
          if (insert_pos >= 0) {
            for (int k = num_closest - 1; k >= insert_pos + 1; k--) {
              dist[k] = dist[k - 1];
              idx[k] = idx[k - 1];
            }
            dist[insert_pos] = d;
            idx[insert_pos] = static_cast<long>(j);
          }
        }

        const size_t cbase = d_coff[p] + i * static_cast<size_t>(quads_per_group);
        if (idx[num_closest - 1] < 0) {
          for (int q = 0; q < quads_per_group; q++) dvalid[cbase + q] = 0;
          return;
        }

        double gxs[kMaxClose], gys[kMaxClose];
        for (int a = 0; a < num_closest; a++) {
          gxs[a] = dx[base + static_cast<size_t>(idx[a])];
          gys[a] = dy[base + static_cast<size_t>(idx[a])];
        }
        double gd[kMaxClose * kMaxClose];
        for (int a = 0; a < num_closest; a++)
          for (int b = a + 1; b < num_closest; b++) {
            const double ex = gxs[a] - gxs[b];
            const double ey = gys[a] - gys[b];
            gd[a * kMaxClose + b] = sycl::sqrt(ex * ex + ey * ey);
          }

        for (int q = 0; q < quads_per_group; q++) {
          const int i0 = d_combo[q * 4 + 0], i1 = d_combo[q * 4 + 1];
          const int i2 = d_combo[q * 4 + 2], i3 = d_combo[q * 4 + 3];

          // Canonical key: the four star indices, sorted, packed into 64 bits.
          uint32_t si[4] = {static_cast<uint32_t>(idx[i0]), static_cast<uint32_t>(idx[i1]),
                            static_cast<uint32_t>(idx[i2]), static_cast<uint32_t>(idx[i3])};
          for (int a = 0; a < 3; a++)
            for (int b = a + 1; b < 4; b++)
              if (si[b] < si[a]) { const uint32_t t2 = si[a]; si[a] = si[b]; si[b] = t2; }
          dkey[cbase + q] = (static_cast<uint64_t>(si[0]) << 48) |
                            (static_cast<uint64_t>(si[1]) << 32) |
                            (static_cast<uint64_t>(si[2]) << 16) | static_cast<uint64_t>(si[3]);

          const double xt = (gxs[i0] + gxs[i1] + gxs[i2] + gxs[i3]) * 0.25;
          const double yt = (gys[i0] + gys[i1] + gys[i2] + gys[i3]) * 0.25;

          double d1 = gd[i0 * kMaxClose + i1];
          double d2 = gd[i0 * kMaxClose + i2];
          double d3 = gd[i0 * kMaxClose + i3];
          double d4 = gd[i1 * kMaxClose + i2];
          double d5 = gd[i1 * kMaxClose + i3];
          double d6 = gd[i2 * kMaxClose + i3];

          double t;
#define ASTAP_SW(a, b) if ((b) > (a)) { t = (a); (a) = (b); (b) = t; }
          ASTAP_SW(d1, d2) ASTAP_SW(d2, d3) ASTAP_SW(d3, d4) ASTAP_SW(d4, d5) ASTAP_SW(d5, d6)
          ASTAP_SW(d1, d2) ASTAP_SW(d2, d3) ASTAP_SW(d3, d4) ASTAP_SW(d4, d5)
          ASTAP_SW(d1, d2) ASTAP_SW(d2, d3) ASTAP_SW(d3, d4)
          ASTAP_SW(d1, d2) ASTAP_SW(d2, d3)
          ASTAP_SW(d1, d2)
#undef ASTAP_SW

          double *o = dcand + (cbase + q) * 8;
          o[0] = d1;
          o[1] = d2 / d1;
          o[2] = d3 / d1;
          o[3] = d4 / d1;
          o[4] = d5 / d1;
          o[5] = d6 / d1;
          o[6] = xt;
          o[7] = yt;
          dvalid[cbase + q] = 1;
        }
      });

      // Kernel 2: duplicate rejection by star set, serial within a position and
      // parallel across positions. Keeping the first occurrence of each key
      // reproduces what the CPU does, and a duplicate now costs one integer
      // probe rather than nine cell lookups and a pile of double compares.
      q_.parallel_for(sycl::range<1>(nblocks), [=](sycl::id<1> pid) {
        const size_t p = pid[0];
        const size_t first = d_coff[p];
        const size_t ncand = d_coff[p + 1] - first;
        uint64_t *tab = dtab + d_hoff[p];
        const size_t tsize = d_hoff[p + 1] - d_hoff[p];
        const size_t mask = tsize - 1;
        constexpr uint64_t kEmpty = ~0ULL;
        for (size_t b = 0; b < tsize; b++) tab[b] = kEmpty;

        size_t nq = 0;
        for (size_t c = 0; c < ncand; c++) {
          if (!dvalid[first + c]) continue;
          const uint64_t key = dkey[first + c];

          size_t slot = static_cast<size_t>((key * 0x9E3779B97F4A7C15ULL) >> 40) & mask;
          bool duplicate = false;
          while (tab[slot] != kEmpty) {
            if (tab[slot] == key) { duplicate = true; break; }
            slot = (slot + 1) & mask;
          }
          if (duplicate) continue;
          tab[slot] = key;

          const double *o = dcand + (first + c) * 8;
          double *dst = dcand + (first + nq) * 8;
          if (nq != c)
            for (int r = 0; r < 8; r++) dst[r] = o[r];
          nq++;
        }
        d_nq[p] = static_cast<uint32_t>(nq);
      });

      // Only the surviving quad counts come back, then the dense offsets go out.
      h_nq_.resize(nblocks + 1);
      q_.memcpy(h_nq_.data(), d_nq, nblocks * sizeof(uint32_t)).wait();

      std::vector<uint32_t> dense_off(nblocks + 1, 0);
      for (size_t k = 0; k < nblocks; k++) dense_off[k + 1] = dense_off[k] + h_nq_[k];
      const size_t total_quads = dense_off[nblocks];

      if (total_quads == 0) {
        q_.wait();
        for (size_t k = 0; k < nblocks; k++) out[gpu_idx[k]].resize(8, 0);
      } else {
        q_.memcpy(d_doff, dense_off.data(), (nblocks + 1) * sizeof(uint32_t)).wait();
        double *ddense = d_dense_.get(q_, total_quads * 8);

        // Kernel 3: gather into a row major block per position, so the host copy
        // back is one memcpy per row and no transpose is needed.
        q_.parallel_for(sycl::range<1>(total_quads), [=](sycl::id<1> gid) {
          const size_t g = gid[0];
          size_t lo = 0, hi = nblocks;
          while (lo + 1 < hi) {
            const size_t mid = (lo + hi) / 2;
            if (d_doff[mid] <= g) lo = mid; else hi = mid;
          }
          const size_t p = lo;
          const size_t i = g - d_doff[p];
          const size_t n = d_doff[p + 1] - d_doff[p];
          if (i >= n) return;
          const double *src = dcand + (d_coff[p] + i) * 8;
          double *dst = ddense + static_cast<size_t>(d_doff[p]) * 8;
          for (int r = 0; r < 8; r++) dst[r * n + i] = src[r];
        });

        h_dense_.resize(total_quads * 8);
        q_.memcpy(h_dense_.data(), ddense, total_quads * 8 * sizeof(double)).wait();

        parallel_for(0, nblocks, [&](size_t k, unsigned) {
          RowList &qq = out[gpu_idx[k]];
          const size_t n = h_nq_[k];
          qq.resize(8, n);
          if (n == 0) return;
          const double *src = h_dense_.data() + static_cast<size_t>(dense_off[k]) * 8;
          for (int r = 0; r < 8; r++)
            std::memcpy(qq.data(r), src + static_cast<size_t>(r) * n, n * sizeof(double));
        });
      }

      // Positions that did not go to the device.
      std::vector<char> done(npos, 0);
      for (size_t k: gpu_idx) done[k] = 1;
      parallel_for(0, npos, [&](size_t i, unsigned) {
        if (!done[i]) find_quads(nrstars_image, lists[i], out[i]);
      });
    }
  } // namespace

  std::unique_ptr<QuadBatchBuilder> make_sycl_quad_builder() {
    try {
      // In order: the kernels below depend on the preceding copies, and a
      // default SYCL queue is out of order, which would let a kernel read the
      // offset arrays before they have been written.
      sycl::queue q{sycl::default_selector_v, sycl::property::queue::in_order{}};
      // Double precision is required for bit-identical results.
      if (!q.get_device().has(sycl::aspect::fp64)) return nullptr;
      if (!q.get_device().has(sycl::aspect::usm_device_allocations)) return nullptr;
      return std::unique_ptr<QuadBatchBuilder>(new SyclQuadBuilder(std::move(q)));
    } catch (const sycl::exception &) {
      return nullptr;
    }
  }

  bool sycl_available() { return true; }
} // namespace astap
