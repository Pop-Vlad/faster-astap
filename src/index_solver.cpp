#include "astap/index_solver.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <unordered_map>
#include <vector>

#include "astap/astro_math.h"
#include "astap/matching.h"

namespace astap {
  namespace {
    struct Pair {
      uint32_t iq;   // image quad
      uint32_t dq;   // database quad
      double scale;  // d1_database / d1_image, arcsec per pixel
      long sbin;     // log-space scale bin
    };

    // A similarity between image pixels and standard coordinates about some
    // reference point, written as a complex map so two correspondences solve it
    // in closed form:
    //
    //   parity +1 (no mirror)   z' = alpha * z      + beta
    //   parity -1 (mirrored)    z' = alpha * conj(z) + beta
    //
    // with z = x + iy in pixels and z' = xi + i*eta in arcseconds. Sky images
    // come in both parities, so both are tried.
    struct Model {
      double p = 0, q = 0;    // alpha = p + iq
      double tx = 0, ty = 0;  // beta
      int parity = 1;
      double scale() const { return std::sqrt(p * p + q * q); }

      void apply(double x, double y, double &xi, double &eta) const {
        if (parity > 0) {
          xi = p * x - q * y + tx;
          eta = q * x + p * y + ty;
        } else {
          xi = p * x + q * y + tx;
          eta = q * x - p * y + ty;
        }
      }
    };

    // Model through two correspondences. Returns false when the two image points
    // are too close for the division to be meaningful.
    bool model_from(double x1, double y1, double xi1, double eta1, double x2, double y2, double xi2,
                    double eta2, int parity, double min_baseline, Model &m) {
      const double dx = parity > 0 ? x2 - x1 : x2 - x1;
      const double dy = parity > 0 ? y2 - y1 : -(y2 - y1);  // conj for the mirrored case
      const double den = dx * dx + dy * dy;
      if (den < min_baseline * min_baseline) return false;

      // alpha = (z2' - z1') / w, with w = (z2 - z1) or conj(z2 - z1).
      const double nx = xi2 - xi1, ny = eta2 - eta1;
      m.p = (nx * dx + ny * dy) / den;
      m.q = (ny * dx - nx * dy) / den;
      m.parity = parity;

      double px, py;
      m.tx = 0;
      m.ty = 0;
      m.apply(x1, y1, px, py);
      m.tx = xi1 - px;
      m.ty = eta1 - py;
      return true;
    }

    // Unit vector mean of a set of sky positions, converted back to RA/Dec. Used
    // for cluster centroids: averaging RA directly breaks at the 0h wraparound
    // and near the poles.
    void mean_direction(const std::vector<double> &ra, const std::vector<double> &dec, double &mra,
                        double &mdec) {
      double vx = 0, vy = 0, vz = 0;
      for (size_t i = 0; i < ra.size(); i++) {
        const double cd = std::cos(dec[i]);
        vx += cd * std::cos(ra[i]);
        vy += cd * std::sin(ra[i]);
        vz += std::sin(dec[i]);
      }
      const double n = std::sqrt(vx * vx + vy * vy + vz * vz);
      if (n <= 0) {
        mra = ra.empty() ? 0 : ra[0];
        mdec = dec.empty() ? 0 : dec[0];
        return;
      }
      mra = std::atan2(vy / n, vx / n);
      if (mra < 0) mra += 2 * kPi;
      mdec = std::asin(std::max(-1.0, std::min(1.0, vz / n)));
    }

    // Solves the same overdetermined system the port builds, with the same
    // routine: A holds the image quad centres in pixels, b the database quad
    // centres in standard coordinates (arcsec) about the reference position.
    bool linear_fit(const QuadIndex &index, const RowList &image_quads,
                    const std::vector<Pair> &pairs, double ref_ra, double ref_dec,
                    SolutionVector &sx, SolutionVector &sy) {
      RowList a(3, pairs.size());
      std::vector<double> bx(pairs.size()), by(pairs.size());
      for (size_t k = 0; k < pairs.size(); k++) {
        a(0, k) = image_quads(6, pairs[k].iq);
        a(1, k) = image_quads(7, pairs[k].iq);
        a(2, k) = 1;
        equatorial_standard(ref_ra, ref_dec, index.centre_ra(pairs[k].dq),
                            index.centre_dec(pairs[k].dq), 1, bx[k], by[k]);
      }
      if (!lsq_fit(a, bx, sx)) return false;
      if (!lsq_fit(a, by, sy)) return false;

      // The port's final sanity check: x and y must carry the same scale, or the
      // fit is not a rigid transformation and the match was coincidence.
      const double xy_ratio = (sqr(sx[0]) + sqr(sx[1])) / (0.00000001 + sqr(sy[0]) + sqr(sy[1]));
      return xy_ratio >= 0.9 && xy_ratio <= 1.1;
    }

    // Turns a fitted solution into a WCS, derived exactly as solve_image() does.
    void wcs_from(const SolutionVector &sx, const SolutionVector &sy, double ref_ra, double ref_dec,
                  int width, int height, IndexSolveResult &out) {
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
      out.scale_arcsec_px = cdelt2_arcsec;
    }

    // Fits one candidate sky position from its consensus set.
    //
    // A refinement pass was tried here and removed: re-scanning every candidate
    // pair at the fitted solution and fitting again changed the corpus mean
    // error by nothing at all (0.248 px either way, 1 image better and 1 worse).
    // The RANSAC consensus already contains every pair the solution explains, so
    // the residual 0.2 px is the precision of a quad centre, not a shortage of
    // points. Do not retry this; a real improvement has to come from the star
    // positions, not the quad centres.
    bool fit_at(const QuadIndex &index, const RowList &image_quads,
                const std::vector<Pair> &consensus, double ref_ra, double ref_dec, int width,
                int height, const IndexSolveSettings &s, IndexSolveResult &out) {
      if (static_cast<int>(consensus.size()) < s.minimum_quads) return false;

      SolutionVector sx, sy;
      if (!linear_fit(index, image_quads, consensus, ref_ra, ref_dec, sx, sy)) return false;

      wcs_from(sx, sy, ref_ra, ref_dec, width, height, out);
      out.nr_references = static_cast<int>(consensus.size());
      out.solved = true;
      return true;
    }

    // A candidate that survived verification, waiting to be fitted.
    struct Candidate {
      std::vector<Pair> inliers;
      double ref_ra = 0, ref_dec = 0;
    };

    // RANSAC over one cluster. Every image point is paired with its database quad
    // centre in standard coordinates about `ref`; a model through two of them is
    // scored by how many of the rest it explains.
    //
    // This is the step the vote alone cannot supply. A cell can hold thirty
    // coincidental pairs that share a plate scale purely by accident, but they
    // cannot also share one rotation and one translation.
    size_t ransac(const QuadIndex &index, const RowList &image_quads,
                  const std::vector<Pair> &cluster, double ref_ra, double ref_dec,
                  double voted_scale, double tol_arcsec, double min_baseline_px,
                  const IndexSolveSettings &s, std::vector<Pair> &best_inliers) {
      const size_t m = cluster.size();
      std::vector<double> px(m), py(m), xi(m), eta(m);
      for (size_t k = 0; k < m; k++) {
        px[k] = image_quads(6, cluster[k].iq);
        py[k] = image_quads(7, cluster[k].iq);
        equatorial_standard(ref_ra, ref_dec, index.centre_ra(cluster[k].dq),
                            index.centre_dec(cluster[k].dq), 1, xi[k], eta[k]);
      }

      // Enumerate every sample for a small cluster; fall back to a deterministic
      // pseudo-random draw when that would be too many. Determinism matters: the
      // solver has to give the same answer for the same input every run.
      const size_t exhaustive = m * (m - 1) / 2;
      const bool enumerate = exhaustive <= static_cast<size_t>(s.max_samples);
      const size_t samples = enumerate ? exhaustive : static_cast<size_t>(s.max_samples);
      uint64_t rng = 0x9E3779B97F4A7C15ULL ^ m;

      size_t best = 0;
      std::vector<size_t> inl, best_idx;
      size_t ii = 0, jj = 1;
      for (size_t t = 0; t < samples; t++) {
        size_t i, j;
        if (enumerate) {
          i = ii;
          j = jj;
          if (++jj >= m) {
            ii++;
            jj = ii + 1;
          }
        } else {
          rng = rng * 6364136223846793005ULL + 1442695040888963407ULL;
          i = (rng >> 33) % m;
          rng = rng * 6364136223846793005ULL + 1442695040888963407ULL;
          j = (rng >> 33) % m;
          if (i == j) continue;
        }
        // Two pairs sharing an image quad describe one point, not two.
        if (cluster[i].iq == cluster[j].iq) continue;

        for (int parity = 1; parity >= -1; parity -= 2) {
          Model mo;
          if (!model_from(px[i], py[i], xi[i], eta[i], px[j], py[j], xi[j], eta[j], parity,
                          min_baseline_px, mo))
            continue;
          // The sample has to reproduce the scale the cluster voted for, which
          // rejects most bad draws before any inlier is counted.
          if (mo.scale() < voted_scale * 0.97 || mo.scale() > voted_scale * 1.03) continue;

          inl.clear();
          for (size_t k = 0; k < m; k++) {
            double ex, ey;
            mo.apply(px[k], py[k], ex, ey);
            if (std::fabs(ex - xi[k]) <= tol_arcsec && std::fabs(ey - eta[k]) <= tol_arcsec)
              inl.push_back(k);
          }
          if (inl.size() > best) {
            best = inl.size();
            best_idx = inl;
          }
        }
      }

      // One image quad may have matched several database quads. Only one of them
      // can be right, so keep the closest and let the fit see independent points.
      best_inliers.clear();
      if (best_idx.empty()) return 0;
      std::unordered_map<uint32_t, size_t> chosen;  // image quad -> position in best_inliers
      for (size_t k : best_idx) {
        auto it = chosen.find(cluster[k].iq);
        if (it == chosen.end()) {
          chosen.emplace(cluster[k].iq, best_inliers.size());
          best_inliers.push_back(cluster[k]);
        }
      }
      return best_inliers.size();
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
      const double d1 = image_quads(0, q);
      if (d1 <= 0) continue;
      float r[5];
      for (int k = 0; k < 5; k++) r[k] = static_cast<float>(image_quads(k + 1, q));
      hits.clear();
      index.query(r, hits);
      for (uint32_t h : hits) {
        const double sc = index.d1(h) / d1;
        pairs.push_back({static_cast<uint32_t>(q), h, sc,
                         static_cast<long>(std::floor(std::log(sc) / s.scale_bin))});
      }
    }
    out.nr_matches = static_cast<int>(pairs.size());
    if (static_cast<int>(pairs.size()) < s.minimum_quads) {
      out.reason = "too few candidate pairs";
      return out;
    }

    // 2. Joint vote on scale *and* position.
    //
    // Voting on scale alone and then filtering, as this used to, throws the
    // solution away whenever a noise bin outranks the true one — which happens,
    // because scale is one dimensional and coincidental matches pile up in it.
    // Measured on the galactic pole field: a pure-noise bin took 22 votes to the
    // true bin's 18, and every true pair was discarded before the position vote
    // ever ran. Voting on the two together spreads the noise across both axes
    // while the true pairs stay in one bucket.
    const double diag_px = std::sqrt(static_cast<double>(width) * width +
                                     static_cast<double>(height) * height);

    // A cell has to be at least as large as the field or the field splits its own
    // vote, and the field's angular size depends on the scale being voted for, so
    // the cell size is a function of the scale bin.
    auto cell_size = [&](long sbin) {
      const double sc = std::exp((sbin + 0.5) * s.scale_bin);
      const double field_rad = diag_px * sc / 3600 * kPi / 180;
      const double floor_rad = s.sky_cell_deg > 0 ? s.sky_cell_deg * kPi / 180 : 0.25 * kPi / 180;
      return std::max(floor_rad, field_rad);
    };

    // Cells are cut on the unit sphere rather than in RA/Dec: no wraparound at 0h
    // and no degeneracy at the poles.
    struct Key {
      long s, x, y, z;
      bool operator==(const Key &o) const { return s == o.s && x == o.x && y == o.y && z == o.z; }
    };
    struct KeyHash {
      size_t operator()(const Key &k) const {
        uint64_t h = 1469598103934665603ULL;
        for (long v : {k.s, k.x, k.y, k.z}) {
          h ^= static_cast<uint64_t>(v);
          h *= 1099511628211ULL;
        }
        return static_cast<size_t>(h);
      }
    };

    std::unordered_map<Key, std::vector<uint32_t>, KeyHash> vote;
    vote.reserve(pairs.size() * 2);
    for (uint32_t k = 0; k < pairs.size(); k++) {
      const double c = cell_size(pairs[k].sbin);
      const double dec = index.centre_dec(pairs[k].dq), ra = index.centre_ra(pairs[k].dq);
      const double cd = std::cos(dec);
      const Key key{pairs[k].sbin, static_cast<long>(std::floor(cd * std::cos(ra) / c)),
                    static_cast<long>(std::floor(cd * std::sin(ra) / c)),
                    static_cast<long>(std::floor(std::sin(dec) / c))};
      vote[key].push_back(k);
    }

    std::vector<std::pair<size_t, const std::vector<uint32_t> *>> peaks;
    peaks.reserve(vote.size());
    for (const auto &kv : vote) peaks.push_back({kv.second.size(), &kv.second});
    // Sort by vote count, breaking ties on the first member so the order does not
    // depend on how the hash table happened to lay itself out.
    std::sort(peaks.begin(), peaks.end(), [](const auto &a, const auto &b) {
      if (a.first != b.first) return a.first > b.first;
      return (*a.second)[0] < (*b.second)[0];
    });

    // 3. Verify the strongest peaks. Verification is the point: a bucket is a
    // hypothesis, not an answer, and one that cannot support a consistent
    // transformation is rejected before least squares ever sees it.
    std::vector<Candidate> candidates;
    const int tries = std::min<int>(s.top_k, static_cast<int>(peaks.size()));
    std::vector<double> cra, cdec;
    std::vector<Pair> cluster, inliers;

    for (int t = 0; t < tries; t++) {
      const std::vector<uint32_t> &members = *peaks[t].second;
      if (static_cast<int>(members.size()) < 2) break;  // nothing below this can help either
      out.peaks_tried = t + 1;

      // Centroid of the bucket, then everything at a compatible scale within one
      // field of it — which recovers pairs the cell boundary happened to cut off.
      cra.clear();
      cdec.clear();
      for (uint32_t k : members) {
        cra.push_back(index.centre_ra(pairs[k].dq));
        cdec.push_back(index.centre_dec(pairs[k].dq));
      }
      double ref_ra, ref_dec;
      mean_direction(cra, cdec, ref_ra, ref_dec);

      const long sbin = pairs[members[0]].sbin;
      const double voted_scale = std::exp((sbin + 0.5) * s.scale_bin);
      const double gather = cell_size(sbin);

      cluster.clear();
      for (const Pair &p : pairs) {
        if (std::labs(p.sbin - sbin) > 1) continue;
        double sep;
        ang_sep(index.centre_ra(p.dq), index.centre_dec(p.dq), ref_ra, ref_dec, sep);
        if (sep <= gather) cluster.push_back(p);
      }
      if (static_cast<int>(cluster.size()) < s.min_inliers) continue;

      const double tol_arcsec = s.inlier_tolerance_px * voted_scale;
      const size_t n = ransac(index, image_quads, cluster, ref_ra, ref_dec, voted_scale, tol_arcsec,
                              0.05 * diag_px, s, inliers);
      if (static_cast<int>(n) < s.min_inliers) continue;

      candidates.push_back({inliers, ref_ra, ref_dec});
      // A consensus this size is not going to be beaten by a later peak, and the
      // whole point of the index solver is that it finishes in milliseconds.
      if (static_cast<int>(n) >= 3 * s.min_inliers) break;
    }

    if (candidates.empty()) {
      out.reason = "no candidate position reached consensus";
      return out;
    }

    // 4. Fit the strongest consensus first. A candidate can still fail the port's
    // scale-consistency check, in which case the next one gets its turn.
    std::stable_sort(candidates.begin(), candidates.end(),
                     [](const Candidate &a, const Candidate &b) {
                       return a.inliers.size() > b.inliers.size();
                     });
    for (const Candidate &c : candidates) {
      out.nr_inliers = static_cast<int>(c.inliers.size());
      if (fit_at(index, image_quads, c.inliers, c.ref_ra, c.ref_dec, width, height, s, out))
        return out;
    }

    out.solved = false;
    out.nr_inliers = 0;
    out.reason = "no consensus survived the fit";
    return out;
  }

  IndexSolveResult solve_with_tiers(const std::vector<QuadIndex> &tiers, const RowList &image_quads,
                                    int width, int height, const IndexSolveSettings &s,
                                    double density_hint) {
    IndexSolveResult best;
    best.reason = "no depth tier solved this image";
    if (tiers.empty()) return best;

    std::vector<size_t> order(tiers.size());
    for (size_t i = 0; i < order.size(); i++) order[i] = i;
    if (density_hint > 0) {
      // Closest tier in log density first: depth has to match to within about a
      // factor of two, so distance in ratio is the meaningful ordering.
      std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) {
        const double da = std::fabs(std::log(tiers[a].settings().star_density / density_hint));
        const double db = std::fabs(std::log(tiers[b].settings().star_density / density_hint));
        return da < db;
      });
    }

    int tried = 0;
    for (size_t i : order) {
      if (tiers[i].size() == 0) continue;
      tried++;
      IndexSolveResult r = solve_with_index(tiers[i], image_quads, width, height, s);
      r.tier_density = tiers[i].settings().star_density;
      if (r.solved && r.nr_inliers > best.nr_inliers) {
        best = r;
        // A consensus this strong is not a coincidence, and the remaining tiers
        // cost milliseconds each that the caller asked the index solver to save.
        if (r.nr_inliers >= 3 * s.min_inliers) break;
      }
    }
    best.tiers_tried = tried;
    return best;
  }
} // namespace astap
