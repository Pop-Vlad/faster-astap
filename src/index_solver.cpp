#include "astap/index_solver.h"

#include <algorithm>
#include <cmath>
#include <map>

#include "astap/astro_math.h"
#include "astap/matching.h"

namespace astap {
  namespace {
    struct Pair {
      uint32_t iq;    // image quad
      uint32_t dq;    // database quad
      double scale;   // d1_database / d1_image, arcsec per pixel
    };

    // Attempts a fit for one candidate sky position: the surviving pairs are
    // turned into the same overdetermined system the port builds, and solved
    // with the same routine.
    bool fit_at(const QuadIndex &index, const RowList &image_quads, const std::vector<Pair> &pairs,
                double ref_ra, double ref_dec, int width, int height,
                const IndexSolveSettings &s, IndexSolveResult &out) {
      if (static_cast<int>(pairs.size()) < s.minimum_quads) return false;

      // A holds the image quad centres in pixels, b the database quad centres in
      // standard coordinates (arcsec) about the reference position. This is
      // exactly the system solve_image() fills from the spiral's matches.
      RowList a(3, pairs.size());
      std::vector<double> bx(pairs.size()), by(pairs.size());
      for (size_t k = 0; k < pairs.size(); k++) {
        a(0, k) = image_quads(6, pairs[k].iq);
        a(1, k) = image_quads(7, pairs[k].iq);
        a(2, k) = 1;
        equatorial_standard(ref_ra, ref_dec, index.centre_ra(pairs[k].dq),
                            index.centre_dec(pairs[k].dq), 1, bx[k], by[k]);
      }

      SolutionVector sx, sy;
      if (!lsq_fit(a, bx, sx)) return false;
      if (!lsq_fit(a, by, sy)) return false;

      // The port's final sanity check: x and y must carry the same scale, or the
      // fit is not a rigid transformation and the match was coincidence.
      const double xy_ratio = (sqr(sx[0]) + sqr(sx[1])) / (0.00000001 + sqr(sy[0]) + sqr(sy[1]));
      if (xy_ratio < 0.9 || xy_ratio > 1.1) return false;

      // WCS, derived exactly as solve_image() does.
      const double centerX = (width - 1) / 2.0;
      const double centerY = (height - 1) / 2.0;

      double ra_solved, dec_solved;
      standard_equatorial(ref_ra, ref_dec, sx[0] * centerX + sx[1] * centerY + sx[2],
                          sy[0] * centerX + sy[1] * centerY + sy[2], 1, ra_solved, dec_solved);

      const double flipped = (sx[0] * sy[1] - sx[1] * sy[0] > 0) ? -1.0 : 1.0;

      double ra7, dec7;
      standard_equatorial(ref_ra, ref_dec, sx[0] * centerX + sx[1] * (centerY + 1) + sx[2],
                          sy[0] * centerX + sy[1] * (centerY + 1) + sy[2], 1, ra7, dec7);
      const double crota2_rad = -position_angle(ra7, dec7, ra_solved, dec_solved);
      const double cdelt1_arcsec = flipped * std::sqrt(sqr(sx[0]) + sqr(sx[1]));
      const double cdelt2_arcsec = std::sqrt(sqr(sy[0]) + sqr(sy[1]));

      standard_equatorial(ref_ra, ref_dec, sx[0] * (centerX + flipped) + sx[1] * centerY + sx[2],
                          sy[0] * (centerX + flipped) + sy[1] * centerY + sy[2], 1, ra7, dec7);
      double crota1_rad = kPi / 2 - position_angle(ra7, dec7, ra_solved, dec_solved);
      if (crota1_rad > kPi) crota1_rad -= 2 * kPi;

      out.solved = true;
      out.ra0 = ra_solved;
      out.dec0 = dec_solved;
      out.crpix1 = centerX + 1;
      out.crpix2 = centerY + 1;
      out.cdelt1 = cdelt1_arcsec / 3600;
      out.cdelt2 = cdelt2_arcsec / 3600;
      out.cd1_1 = +out.cdelt1 * std::cos(crota1_rad);
      out.cd1_2 = -out.cdelt1 * std::sin(crota1_rad) * flipped;
      out.cd2_1 = +out.cdelt2 * std::sin(crota2_rad) * flipped;
      out.cd2_2 = +out.cdelt2 * std::cos(crota2_rad);
      out.crota1 = crota1_rad * 180 / kPi;
      out.crota2 = crota2_rad * 180 / kPi;
      out.nr_references = static_cast<int>(pairs.size());
      out.scale_arcsec_px = cdelt2_arcsec;
      return true;
    }
  } // namespace

  IndexSolveResult solve_with_index(const QuadIndex &index, const RowList &image_quads, int width,
                                    int height, const IndexSolveSettings &s) {
    IndexSolveResult out;
    if (image_quads.count() == 0 || index.size() == 0) {
      out.reason = "no quads";
      return out;
    }

    // 1. One pass over the image quads against the whole sky.
    std::vector<Pair> pairs;
    std::vector<uint32_t> hits;
    for (size_t q = 0; q < image_quads.count(); q++) {
      float r[5];
      for (int k = 0; k < 5; k++) r[k] = static_cast<float>(image_quads(k + 1, q));
      hits.clear();
      index.query(r, hits);
      const double d1 = image_quads(0, q);
      if (d1 <= 0) continue;
      for (uint32_t h : hits)
        pairs.push_back({static_cast<uint32_t>(q), h, index.d1(h) / d1});
    }
    out.nr_matches = static_cast<int>(pairs.size());
    if (static_cast<int>(pairs.size()) < s.minimum_quads) {
      out.reason = "too few candidate pairs";
      return out;
    }

    // 2. Vote on plate scale. A true solution puts all of its pairs at one
    // scale; coincidental matches scatter. This is the port's median-ratio
    // filter, applied over the whole sky at once instead of per position.
    std::map<long, int> scale_hist;
    for (const Pair &p : pairs)
      scale_hist[static_cast<long>(std::floor(std::log(p.scale) / s.scale_bin))]++;
    long best_bin = 0;
    int best_n = 0;
    for (const auto &kv : scale_hist)
      if (kv.second > best_n) {
        best_n = kv.second;
        best_bin = kv.first;
      }
    const double peak_scale = std::exp((best_bin + 0.5) * s.scale_bin);

    std::vector<Pair> at_scale;
    for (const Pair &p : pairs)
      if (std::fabs(p.scale - peak_scale) < 2 * s.scale_bin * peak_scale) at_scale.push_back(p);
    if (static_cast<int>(at_scale.size()) < s.minimum_quads) {
      out.reason = "no scale consensus";
      return out;
    }

    // 3. Vote on position: every surviving pair says where on the sky its image
    // quad sits, so the true ones pile into one cell.
    const double cell = s.sky_cell_deg * kPi / 180;
    auto cell_of = [&](size_t k) {
      const double dec = index.centre_dec(at_scale[k].dq);
      const long cy = static_cast<long>(std::floor(dec / cell));
      const long cx = static_cast<long>(
        std::floor(index.centre_ra(at_scale[k].dq) * std::cos(dec) / cell));
      return std::make_pair(cx, cy);
    };
    std::map<std::pair<long, long>, std::vector<size_t>> sky;
    for (size_t k = 0; k < at_scale.size(); k++) sky[cell_of(k)].push_back(k);

    std::vector<std::pair<size_t, std::pair<long, long>>> peaks;
    for (const auto &kv : sky) peaks.push_back({kv.second.size(), kv.first});
    std::sort(peaks.rbegin(), peaks.rend());

    // 4. Try the strongest peaks in turn. The best cell is not always right when
    // the vote is thin, and a wrong one is cheap to reject: the fit's scale
    // consistency check throws it out.
    const int tries = std::min<int>(s.top_k, static_cast<int>(peaks.size()));
    for (int t = 0; t < tries; t++) {
      const auto &members = sky[peaks[t].second];
      // Include the neighbouring cells: a field lying on a cell boundary would
      // otherwise have its votes split.
      std::vector<Pair> cluster;
      const long cx = peaks[t].second.first, cy = peaks[t].second.second;
      for (size_t k = 0; k < at_scale.size(); k++) {
        const auto c = cell_of(k);
        if (std::labs(c.first - cx) <= 1 && std::labs(c.second - cy) <= 1)
          cluster.push_back(at_scale[k]);
      }
      if (static_cast<int>(cluster.size()) < s.minimum_quads) continue;

      double ref_ra = 0, ref_dec = 0;
      for (size_t k : members) {
        ref_ra += index.centre_ra(at_scale[k].dq);
        ref_dec += index.centre_dec(at_scale[k].dq);
      }
      ref_ra /= members.size();
      ref_dec /= members.size();

      if (fit_at(index, image_quads, cluster, ref_ra, ref_dec, width, height, s, out)) return out;
    }

    out.reason = "no position passed the fit";
    return out;
  }
} // namespace astap
