#include "astap/matching.h"

#include <algorithm>
#include <cmath>

#include "astap/astro_math.h"

namespace astap {
  namespace {
    void say(const std::function < void(const std::string &) > &log, const std::string & s) {
      if (log) log(s);
    }

    // Shared tail of find_fit and find_fit_using_hash: throw out the quads whose
    // longest side deviates from the median ratio, then fill the equations.
    bool finish_fit(MatchState &st, const std::vector<size_t> &match1, const std::vector<size_t> &match2,
                    double quad_tolerance, const std::function<void(const std::string &)> &log) {
      const int n = st.nr_references2;

      // Calculate the median of the longest length ratio for the matching quads.
      std::vector<double> ratios(static_cast<size_t>(n));
      for (int k = 0; k < n; k++)
        ratios[static_cast<size_t>(k)] = st.quad_star_distances1(0, match1[static_cast<size_t>(k)]) /
                                         st.quad_star_distances2(0, match2[static_cast<size_t>(k)]);

      // SMedian sorts in place, so keep `ratios` in match order and sort a copy.
      std::vector<double> ratios_sorted = ratios;
      double median_ratio = smedian(ratios_sorted, static_cast<size_t>(n));

      st.nr_references = 0;
      std::vector<size_t> keep1(static_cast<size_t>(n));
      std::vector<size_t> keep2(static_cast<size_t>(n));
      for (int k = 0; k < n; k++) {
        if (std::fabs(median_ratio - ratios[static_cast<size_t>(k)]) <= quad_tolerance * median_ratio) {
          keep1[static_cast<size_t>(st.nr_references)] = match1[static_cast<size_t>(k)];
          keep2[static_cast<size_t>(st.nr_references)] = match2[static_cast<size_t>(k)];
          st.nr_references++;
        } else if (log) {
          say(log, "Quad outlier removed due to abnormal size: " +
                   float_to_str(100 * ratios[static_cast<size_t>(k)] / median_ratio, 2) + "%");
        }
      }
      // Outliers in the longest length removed.

      if (st.nr_references < 3) return false; // 3 quad centre positions are the bare minimum

      st.a_xy_positions.resize(3, static_cast<size_t>(st.nr_references));
      st.b_xrefpositions.assign(static_cast<size_t>(st.nr_references), 0.0);
      st.b_yrefpositions.assign(static_cast<size_t>(st.nr_references), 0.0);

      for (int k = 0; k < st.nr_references; k++) {
        size_t kk = static_cast<size_t>(k);
        st.a_xy_positions(0, kk) = st.quad_star_distances2(6, keep2[kk]); // average x of the quad
        st.a_xy_positions(1, kk) = st.quad_star_distances2(7, keep2[kk]); // average y of the quad
        st.a_xy_positions(2, kk) = 1;

        st.b_xrefpositions[kk] = st.quad_star_distances1(6, keep1[kk]); // x of the database quad
        st.b_yrefpositions[kk] = st.quad_star_distances1(7, keep1[kk]); // y of the database quad
      }
      return true;
    }
  } // namespace

  bool lsq_fit(const RowList &a_matrix, std::vector<double> b_matrix, SolutionVector &x_matrix) {
    constexpr double tiny = 1E-10; // accuracy

    const size_t nr_equations = a_matrix.count();
    const int nr_columns = a_matrix.rows(); // should be 3 for this application

    // Duplicate a_matrix to protect the caller's data.
    RowList temp_matrix = a_matrix;

    for (int j = 0; j < nr_columns; j++) {
      // Eliminate the matrix elements A[i,j] with i>j from column j.
      for (size_t i = static_cast<size_t>(j) + 1; i < nr_equations; i++) {
        if (temp_matrix(j, i) == 0) continue;

        double p, q, h;
        if (std::fabs(temp_matrix(j, static_cast<size_t>(j))) < tiny * std::fabs(temp_matrix(j, i))) {
          p = 0;
          q = 1;
          temp_matrix(j, static_cast<size_t>(j)) = -temp_matrix(j, i);
          temp_matrix(j, i) = 0;
        } else {
          // Zero the left bottom corner of the matrix. Residuals are r1..rn and
          // the sum of sqr(residuals) has to be minimised. Take two numbers where
          // (p^2+q^2) = 1, then (r1^2+r2^2) = (p^2+q^2)*(r1^2+r2^2).
          // Choose p = +A11/h and q = -A21/h where h = +-sqrt(A11^2+A21^2), so
          // A21 = q*A11+p*A21 = (-A21*A11 + A21*A11)/h = 0.
          h = std::sqrt(temp_matrix(j, static_cast<size_t>(j)) * temp_matrix(j, static_cast<size_t>(j)) +
                        temp_matrix(j, i) * temp_matrix(j, i));
          if (temp_matrix(j, static_cast<size_t>(j)) < 0) h = -h;
          p = temp_matrix(j, static_cast<size_t>(j)) / h;
          q = -temp_matrix(j, i) / h;
          temp_matrix(j, static_cast<size_t>(j)) = h;
          temp_matrix(j, i) = 0;
        }

        // Calculate the rest of the line.
        for (int k = j + 1; k < nr_columns; k++) {
          h = p * temp_matrix(k, static_cast<size_t>(j)) - q * temp_matrix(k, i);
          temp_matrix(k, i) = q * temp_matrix(k, static_cast<size_t>(j)) + p * temp_matrix(k, i);
          temp_matrix(k, static_cast<size_t>(j)) = h;
        }
        h = p * b_matrix[static_cast<size_t>(j)] - q * b_matrix[i];
        b_matrix[i] = q * b_matrix[static_cast<size_t>(j)] + p * b_matrix[i];
        b_matrix[static_cast<size_t>(j)] = h;
      }
    }

    for (int i = 0; i < nr_columns; i++) x_matrix[i] = 0;

    for (int i = nr_columns - 1; i >= 0; i--) {
      // back substitution
      double h = b_matrix[static_cast<size_t>(i)];
      for (int k = i + 1; k < nr_columns; k++) h -= temp_matrix(k, static_cast<size_t>(i)) * x_matrix[k];
      if (std::fabs(temp_matrix(i, static_cast<size_t>(i))) > 1E-30)
        x_matrix[i] = h / temp_matrix(i, static_cast<size_t>(i));
      else
        return false; // prevent a division by zero, force a failure instead
      // solution vector x := x_matrix[0]*x + x_matrix[1]*y + x_matrix[2]
    }
    return true;
  }

  bool find_fit(MatchState &st, int minimum_count, double quad_tolerance,
                const std::function<void(const std::string &)> &log) {
    const size_t nrquads1 = st.quad_star_distances1.count();
    const size_t nrquads2 = st.quad_star_distances2.count();

    // minimum_count required, 6 for stacking, 3 for plate solving.
    if (nrquads1 < static_cast<size_t>(minimum_count) ||
        nrquads2 < static_cast<size_t>(minimum_count)) {
      st.nr_references = 0;
      return false;
    }

    std::vector<size_t> match1, match2;
    match1.reserve(1000);
    match2.reserve(1000);

    // Cache the image side row pointers, this turns a three level dereference
    // into a single one.
    const double *pImg1 = st.quad_star_distances2.data(1);
    const double *pImg2 = st.quad_star_distances2.data(2);
    const double *pImg3 = st.quad_star_distances2.data(3);
    const double *pImg4 = st.quad_star_distances2.data(4);
    const double *pImg5 = st.quad_star_distances2.data(5);

    st.nr_references2 = 0;
    for (size_t i = 0; i < nrquads1; i++) {
      // Load the database side ratios once per i, reused for all inner iterations.
      double db_r1 = st.quad_star_distances1(1, i);
      double db_r2 = st.quad_star_distances1(2, i);
      double db_r3 = st.quad_star_distances1(3, i);
      double db_r4 = st.quad_star_distances1(4, i);
      double db_r5 = st.quad_star_distances1(5, i);

      for (size_t j = 0; j < nrquads2; j++) {
        //       ==database==                ==image==
        // The short-circuit chain is deliberate: with a tolerance of ~0.007 and
        // ratios in [0,1] roughly 98.6% of the pairs already fail the first check.
        // Combining the five checks into one max() would force all five abs()
        // calls per pair, about 5x slower in the common case.
        if (std::fabs(db_r1 - pImg1[j]) <= quad_tolerance) // all lengths are scaled to the longest,
          if (std::fabs(db_r2 - pImg2[j]) <= quad_tolerance) // so this is scale independent
            if (std::fabs(db_r3 - pImg3[j]) <= quad_tolerance)
              if (std::fabs(db_r4 - pImg4[j]) <= quad_tolerance)
                if (std::fabs(db_r5 - pImg5[j]) <= quad_tolerance) {
                  match1.push_back(i); // store the match position
                  match2.push_back(j);
                  st.nr_references2++;
                }
      }
    }

    if (log) say(log, "Found " + std::to_string(st.nr_references2) + " references");

    if (st.nr_references2 < minimum_count) {
      st.nr_references = 0;
      return false;
    }
    return finish_fit(st, match1, match2, quad_tolerance, log);
  }

  bool find_fit_using_hash(MatchState &st, int minimum_count, double quad_tolerance,
                           const std::function<void(const std::string &)> &log) {
    constexpr int kNeighborBins = 1; // check +/-1 bin to cover the quad tolerance

    const size_t nrquads1 = st.quad_star_distances1.count();
    const size_t nrquads2 = st.quad_star_distances2.count();

    if (nrquads1 < static_cast<size_t>(minimum_count) ||
        nrquads2 < static_cast<size_t>(minimum_count)) {
      st.nr_references = 0;
      return false;
    }

    // The bin index is trunc(ratio[1] / quad_tolerance) and ratio[1] is in (0..1],
    // so the bin can only take about 1/quad_tolerance distinct values (~143 for
    // the default tolerance of 0.007). The +2 covers the slot reached by
    // ratio = 1.0 and gives one extra slot for the neighbour extension.
    int hash_bins = static_cast<int>(pround(1.0 / quad_tolerance)) + 2;
    if (hash_bins < 4) hash_bins = 4; // safety floor for large tolerances
    if (hash_bins > 10000) hash_bins = 10000; // safety cap for tiny tolerances

    std::vector<std::vector<size_t> > hash_table1(static_cast<size_t>(hash_bins));
    std::vector<std::vector<size_t> > hash_table2(static_cast<size_t>(hash_bins));

    size_t max_hash_count = 0;

    for (size_t i = 0; i < nrquads1; i++) {
      long bin = ptrunc(st.quad_star_distances1(1, i) / quad_tolerance);
      // Defensive clamp against numerical noise pushing a ratio above 1.
      if (bin < 0) bin = 0;
      if (bin >= hash_bins) bin = hash_bins - 1;
      hash_table1[static_cast<size_t>(bin)].push_back(i);
      max_hash_count = std::max(max_hash_count, hash_table1[static_cast<size_t>(bin)].size());
    }
    for (size_t j = 0; j < nrquads2; j++) {
      long bin = ptrunc(st.quad_star_distances2(1, j) / quad_tolerance);
      if (bin < 0) bin = 0;
      if (bin >= hash_bins) bin = hash_bins - 1;
      hash_table2[static_cast<size_t>(bin)].push_back(j);
      max_hash_count = std::max(max_hash_count, hash_table2[static_cast<size_t>(bin)].size());
    }

    std::vector<size_t> match1, match2;
    match1.reserve(nrquads1);
    match2.reserve(nrquads1);
    st.nr_references2 = 0;

    for (int bin = 0; bin < hash_bins; bin++) {
      if (hash_table1[static_cast<size_t>(bin)].empty()) continue;
      for (int delta_bin = -kNeighborBins; delta_bin <= kNeighborBins; delta_bin++) {
        // No wraparound: the ratio space [0,1] is not circular, so simply skip
        // out-of-range neighbours.
        int adjusted_bin = bin + delta_bin;
        if (adjusted_bin < 0 || adjusted_bin >= hash_bins) continue;
        if (hash_table2[static_cast<size_t>(adjusted_bin)].empty()) continue;
        for (size_t i: hash_table1[static_cast<size_t>(bin)]) {
          for (size_t j: hash_table2[static_cast<size_t>(adjusted_bin)]) {
            if (std::fabs(st.quad_star_distances1(1, i) - st.quad_star_distances2(1, j)) <= quad_tolerance)
              if (std::fabs(st.quad_star_distances1(2, i) - st.quad_star_distances2(2, j)) <= quad_tolerance)
                if (std::fabs(st.quad_star_distances1(3, i) - st.quad_star_distances2(3, j)) <= quad_tolerance)
                  if (std::fabs(st.quad_star_distances1(4, i) - st.quad_star_distances2(4, j)) <= quad_tolerance)
                    if (std::fabs(st.quad_star_distances1(5, i) - st.quad_star_distances2(5, j)) <= quad_tolerance) {
                      match1.push_back(i);
                      match2.push_back(j);
                      st.nr_references2++;
                    }
          }
        }
      }
    }

    if (log)
      say(log, "Found " + std::to_string(st.nr_references2) +
               " references, max hash bin size: " + std::to_string(max_hash_count));

    if (st.nr_references2 < minimum_count) {
      st.nr_references = 0;
      return false;
    }
    return finish_fit(st, match1, match2, quad_tolerance, log);
  }

  void reset_solution_vectors(MatchState &st, double factor) {
    st.solution_vector_x[0] = factor; // should be one; x := s[1]x + s[2]y + s[3]
    st.solution_vector_x[1] = 0;
    st.solution_vector_x[2] = 0;

    st.solution_vector_y[0] = 0; // y := s[1]x + s[2]y + s[3]
    st.solution_vector_y[1] = factor;
    st.solution_vector_y[2] = 0;
  }

  bool find_offset_and_rotation(MatchState &st, int minimum_quads, double tolerance,
                                const std::function<void(const std::string &)> &log) {
    const size_t nrquads = st.quad_star_distances1.count();
    // 3 quads are required, giving 3 quad centre references.
    if (nrquads < 120) {
      // use the brute force method
      if (!find_fit(st, minimum_quads, tolerance, log)) {
        reset_solution_vectors(st, 0.001); // nullify
        return false;
      }
    } else {
      // use the hash based routine
      if (!find_fit_using_hash(st, minimum_quads, tolerance, log)) {
        reset_solution_vectors(st, 0.001);
        return false;
      }
    }

    // In matrix calculations,
    //   b_refpositionX[0..2, 0..nr_equations-1] := solution_vectorX[0..2] * A_XYpositions[...]
    //   b_refpositionY[0..2, 0..nr_equations-1] := solution_vectorY[0..2] * A_XYpositions[...]

    // Solution vector for X := ax + by + c
    if (!lsq_fit(st.a_xy_positions, st.b_xrefpositions, st.solution_vector_x)) {
      reset_solution_vectors(st, 0.001);
      return false;
    }
    // Solution vector for Y := ax + by + c
    if (!lsq_fit(st.a_xy_positions, st.b_yrefpositions, st.solution_vector_y)) {
      reset_solution_vectors(st, 0.001);
      return false;
    }

    double xy_sqr_ratio = (sqr(st.solution_vector_x[0]) + sqr(st.solution_vector_x[1])) /
                          (0.00000001 + sqr(st.solution_vector_y[0]) + sqr(st.solution_vector_y[1]));
    if (xy_sqr_ratio < 0.9 || xy_sqr_ratio > 1.1) {
      // The dimensions x and y are not the same, something is wrong.
      reset_solution_vectors(st, 0.001);
      if (log) say(log, "Solution skipped on XY ratio: " + float_to_str(xy_sqr_ratio, 6));
      return false;
    }
    return true;
  }
} // namespace astap
