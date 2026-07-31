// SYCL offload of the batched quad construction.
//
// Portable SYCL 2020, no vendor extensions, so it builds under both oneAPI
// DPC++ and AdaptiveCpp.
//
// The results must match the CPU implementation bit for bit, which shapes the
// design in two ways:
//
//   * Everything is double precision and the arithmetic is written in the same
//     order as find_many_quads(). The build passes -ffp-contract=off so the
//     device cannot fuse a*b+c into an FMA and change the last bit.
//
//   * Duplicate quad rejection is a greedy scan: a candidate is dropped when an
//     *already accepted* quad sits within one pixel, so the outcome depends on
//     the order candidates are visited in. That part therefore stays serial per
//     position — but positions are independent, so a batch still gives the
//     device thousands of independent serial scans to run at once.
//
// Only the find_many_quads path (groups of 5, 6 or 7 closest stars) is
// offloaded. That is the path a blind solve spends most of its time in, because
// the early field-of-view attempts detect few stars. The regular find_quads
// path sorts the star list, which is a poor fit here, and falls back to the CPU.

#include <sycl/sycl.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

#include "astap/parallel.h"
#include "astap/quad_batch.h"
#include "astap/quads.h"

namespace astap {
  namespace {
    // The same C(n,4) index tables the CPU uses, flattened for the device.
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

    // Builds the combination table for `n` choose 4 in ascending order, which is
    // exactly the order the CPU tables are written in.
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

    class SyclQuadBuilder : public QuadBatchBuilder {
    public:
      explicit SyclQuadBuilder(sycl::queue q) : q_(std::move(q)) {
        name_ = "sycl:" + q_.get_device().get_info<sycl::info::device::name>();
        cpu_fallback_ = make_cpu_quad_builder();
      }

      std::string backend() const override { return name_; }

      void build(int nrstars_image, std::vector<RowList> &lists, std::vector<RowList> &out) override;

    private:
      sycl::queue q_;
      std::string name_;
      std::unique_ptr<QuadBatchBuilder> cpu_fallback_;
    };

    void SyclQuadBuilder::build(int nrstars_image, std::vector<RowList> &lists,
                                std::vector<RowList> &out) {
      const size_t npos = lists.size();
      out.resize(npos);
      if (npos == 0) return;

      // Split the batch: positions whose star count selects the find_many_quads
      // path go to the device, the rest keep the CPU path.
      std::vector<size_t> gpu_idx;
      gpu_idx.reserve(npos);
      GroupSpec spec{0, 0};
      for (size_t i = 0; i < npos; i++) {
        const GroupSpec s = spec_for(nrstars_image, lists[i].count());
        if (s.num_closest == 0) continue;
        // Every position of a batch sees the same nrstars_image, so the group
        // size only varies when the star counts straddle a threshold. Keep the
        // first one seen and leave the odd ones out to the CPU.
        if (spec.num_closest == 0) spec = s;
        if (s.num_closest == spec.num_closest) gpu_idx.push_back(i);
      }

      if (gpu_idx.empty()) {
        cpu_fallback_->build(nrstars_image, lists, out);
        return;
      }

      const int num_closest = spec.num_closest;
      const int quads_per_group = spec.quads_per_group;

      // Flatten the star lists of the offloaded positions.
      std::vector<uint32_t> counts(gpu_idx.size());
      std::vector<uint32_t> star_off(gpu_idx.size() + 1, 0);
      std::vector<uint32_t> cand_off(gpu_idx.size() + 1, 0);
      for (size_t k = 0; k < gpu_idx.size(); k++) {
        const uint32_t n = static_cast<uint32_t>(lists[gpu_idx[k]].count());
        counts[k] = n;
        star_off[k + 1] = star_off[k] + n;
        cand_off[k + 1] = cand_off[k] + n * static_cast<uint32_t>(quads_per_group);
      }
      const size_t total_stars = star_off.back();
      const size_t total_cand = cand_off.back();
      if (total_stars == 0 || total_cand == 0) {
        cpu_fallback_->build(nrstars_image, lists, out);
        return;
      }

      std::vector<double> hx(total_stars), hy(total_stars);
      for (size_t k = 0; k < gpu_idx.size(); k++) {
        const RowList &l = lists[gpu_idx[k]];
        for (size_t i = 0; i < l.count(); i++) {
          hx[star_off[k] + i] = l(0, i);
          hy[star_off[k] + i] = l(1, i);
        }
      }

      std::vector<int> combos;
      fill_combinations(num_closest, combos);

      // Device buffers. Candidates carry the six hash values plus the centre,
      // laid out as eight doubles per candidate exactly like a quad column.
      sycl::buffer<double> b_x(hx.data(), sycl::range<1>(total_stars));
      sycl::buffer<double> b_y(hy.data(), sycl::range<1>(total_stars));
      sycl::buffer<uint32_t> b_counts(counts.data(), sycl::range<1>(counts.size()));
      sycl::buffer<uint32_t> b_soff(star_off.data(), sycl::range<1>(star_off.size()));
      sycl::buffer<uint32_t> b_coff(cand_off.data(), sycl::range<1>(cand_off.size()));
      sycl::buffer<int> b_combo(combos.data(), sycl::range<1>(combos.size()));

      std::vector<double> cand(total_cand * 8);
      std::vector<uint8_t> valid(total_cand, 0);
      sycl::buffer<double> b_cand(cand.data(), sycl::range<1>(cand.size()));
      sycl::buffer<uint8_t> b_valid(valid.data(), sycl::range<1>(valid.size()));

      const size_t nblocks = gpu_idx.size();

      // Kernel: one work item per (position, reference star). Finds the closest
      // `num_closest` stars, computes the group's pairwise distances once, then
      // emits one candidate per combination.
      q_.submit([&](sycl::handler &h) {
        auto x = b_x.get_access<sycl::access_mode::read>(h);
        auto y = b_y.get_access<sycl::access_mode::read>(h);
        auto cnt = b_counts.get_access<sycl::access_mode::read>(h);
        auto soff = b_soff.get_access<sycl::access_mode::read>(h);
        auto coff = b_coff.get_access<sycl::access_mode::read>(h);
        auto combo = b_combo.get_access<sycl::access_mode::read>(h);
        auto outc = b_cand.get_access<sycl::access_mode::write>(h);
        auto outv = b_valid.get_access<sycl::access_mode::write>(h);

        // A flat range over all reference stars of the batch; each item looks up
        // which position it belongs to from the offsets.
        h.parallel_for(sycl::range<1>(total_stars), [=](sycl::id<1> gid) {
          const size_t g = gid[0];

          // Binary search for the position owning this star.
          size_t lo = 0, hi = nblocks;
          while (lo + 1 < hi) {
            const size_t mid = (lo + hi) / 2;
            if (soff[mid] <= g)
              lo = mid;
            else
              hi = mid;
          }
          const size_t p = lo;
          const size_t base = soff[p];
          const size_t n = cnt[p];
          const size_t i = g - base;
          if (i >= n) return;

          // Closest `num_closest` stars, insertion sorted, exactly as the CPU does.
          long idx[kMaxClose];
          double dist[kMaxClose];
          idx[0] = static_cast<long>(i);
          dist[0] = 0;
          for (int k = 1; k < num_closest; k++) {
            idx[k] = -1;
            dist[k] = 1E99;
          }

          const double x1 = x[base + i];
          const double y1 = y[base + i];
          for (size_t j = 0; j < n; j++) {
            if (j == i) continue;
            const double dx = x[base + j] - x1;
            const double dy = y[base + j] - y1;
            const double d = dx * dx + dy * dy;
            if (d <= 1) continue; // not an identical star

            int insert_pos = -1;
            for (int k = num_closest - 1; k >= 1; k--) {
              if (d < dist[k])
                insert_pos = k;
              else
                break;
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

          const size_t cbase = coff[p] + i * static_cast<size_t>(quads_per_group);
          if (idx[num_closest - 1] < 0) {
            for (int q = 0; q < quads_per_group; q++) outv[cbase + q] = 0;
            return;
          }

          // Group coordinates and their pairwise distances, computed once.
          double gxs[kMaxClose], gys[kMaxClose];
          for (int a = 0; a < num_closest; a++) {
            gxs[a] = x[base + static_cast<size_t>(idx[a])];
            gys[a] = y[base + static_cast<size_t>(idx[a])];
          }
          double gd[kMaxClose * kMaxClose];
          for (int a = 0; a < num_closest; a++)
            for (int b = a + 1; b < num_closest; b++) {
              const double dx = gxs[a] - gxs[b];
              const double dy = gys[a] - gys[b];
              gd[a * kMaxClose + b] = sycl::sqrt(dx * dx + dy * dy);
            }

          for (int q = 0; q < quads_per_group; q++) {
            const int i0 = combo[q * 4 + 0], i1 = combo[q * 4 + 1];
            const int i2 = combo[q * 4 + 2], i3 = combo[q * 4 + 3];

            const double xt = (gxs[i0] + gxs[i1] + gxs[i2] + gxs[i3]) * 0.25;
            const double yt = (gys[i0] + gys[i1] + gys[i2] + gys[i3]) * 0.25;

            double d1 = gd[i0 * kMaxClose + i1];
            double d2 = gd[i0 * kMaxClose + i2];
            double d3 = gd[i0 * kMaxClose + i3];
            double d4 = gd[i1 * kMaxClose + i2];
            double d5 = gd[i1 * kMaxClose + i3];
            double d6 = gd[i2 * kMaxClose + i3];

            // The same unrolled descending sort network as the CPU.
            double t;
#define ASTAP_SW(a, b) if ((b) > (a)) { t = (a); (a) = (b); (b) = t; }
            ASTAP_SW(d1, d2) ASTAP_SW(d2, d3) ASTAP_SW(d3, d4) ASTAP_SW(d4, d5) ASTAP_SW(d5, d6)
            ASTAP_SW(d1, d2) ASTAP_SW(d2, d3) ASTAP_SW(d3, d4) ASTAP_SW(d4, d5)
            ASTAP_SW(d1, d2) ASTAP_SW(d2, d3) ASTAP_SW(d3, d4)
            ASTAP_SW(d1, d2) ASTAP_SW(d2, d3)
            ASTAP_SW(d1, d2)
#undef ASTAP_SW

            double *o = &outc[(cbase + q) * 8];
            o[0] = d1;
            o[1] = d2 / d1;
            o[2] = d3 / d1;
            o[3] = d4 / d1;
            o[4] = d5 / d1;
            o[5] = d6 / d1;
            o[6] = xt;
            o[7] = yt;
            outv[cbase + q] = 1;
          }
        });
      });
      q_.wait();

      // Host side: the greedy duplicate rejection, in candidate order. This is
      // the order dependent part; it is cheap and stays exact here.
      {
        sycl::host_accessor hc(b_cand, sycl::read_only);
        sycl::host_accessor hv(b_valid, sycl::read_only);

        parallel_for(0, gpu_idx.size(), [&](size_t k, unsigned) {
          RowList &q = out[gpu_idx[k]];
          const size_t first = cand_off[k];
          const size_t ncand = cand_off[k + 1] - first;
          q.resize(8, ncand);

          // Same grid as find_many_quads: a cell size of one means two centres
          // within a pixel are at most one cell apart, so probing the 3x3
          // neighbourhood is exactly equivalent to scanning every accepted quad,
          // without the quadratic cost.
          static thread_local std::vector<int32_t> heads;
          static thread_local std::vector<int32_t> next_of;
          size_t bucket_count = 64;
          while (bucket_count < ncand) bucket_count <<= 1;
          const size_t bucket_mask = bucket_count - 1;
          heads.assign(bucket_count, -1);
          next_of.resize(ncand);
          auto cell_of = [bucket_mask](long gx, long gy) {
            return static_cast<size_t>(static_cast<unsigned long>(gx) * 73856093UL ^
                                       static_cast<unsigned long>(gy) * 19349663UL) &
                   bucket_mask;
          };

          size_t nrquads = 0;
          for (size_t c = 0; c < ncand; c++) {
            if (!hv[first + c]) continue;
            const double *o = &hc[(first + c) * 8];
            const double xt = o[6], yt = o[7];

            const long gx = static_cast<long>(std::floor(xt));
            const long gy = static_cast<long>(std::floor(yt));
            bool identical = false;
            const double *qx = q.data(6);
            const double *qy = q.data(7);
            for (long dy = -1; dy <= 1 && !identical; dy++)
              for (long dx = -1; dx <= 1 && !identical; dx++)
                for (int32_t e = heads[cell_of(gx + dx, gy + dy)]; e >= 0; e = next_of[e]) {
                  if (std::fabs(xt - qx[e]) < 1 && std::fabs(yt - qy[e]) < 1) {
                    identical = true;
                    break;
                  }
                }
            if (identical) continue;

            for (int r = 0; r < 8; r++) q(r, nrquads) = o[r];
            const size_t b = cell_of(gx, gy);
            next_of[nrquads] = heads[b];
            heads[b] = static_cast<int32_t>(nrquads);
            nrquads++;
          }
          q.resize(8, nrquads);
        });
      }

      // Positions that did not go to the device.
      std::vector<char> done(npos, 0);
      for (size_t k : gpu_idx) done[k] = 1;
      parallel_for(0, npos, [&](size_t i, unsigned) {
        if (!done[i]) find_quads(nrstars_image, lists[i], out[i]);
      });
    }
  } // namespace

  std::unique_ptr<QuadBatchBuilder> make_sycl_quad_builder() {
    try {
      sycl::queue q{sycl::default_selector_v};
      // Double precision is required for bit-identical results.
      if (!q.get_device().has(sycl::aspect::fp64)) return nullptr;
      return std::unique_ptr<QuadBatchBuilder>(new SyclQuadBuilder(std::move(q)));
    } catch (const sycl::exception &) {
      return nullptr;
    }
  }

  bool sycl_available() { return true; }
} // namespace astap
