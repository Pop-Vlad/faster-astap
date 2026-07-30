#include "astap/calc_trans_cubic.h"

#include <algorithm>
#include <cmath>

namespace astap {
  namespace {
    constexpr double kMatrixTol = 1E-100;
    constexpr int kNum = 10; // ten coefficients per axis

    using Matrix = double[kNum][kNum];
    using Vector10 = double[kNum];

    // Finds the largest value in the matrix at or below `row` and switches rows so
    // that it ends up in `row`. Called by gauss_matrix.
    void gauss_pivot(Matrix matrix, int num, Vector10 vec, double *biggest_val, int row) {
      int pivot_row = row;
      double big = std::fabs(matrix[row][row] / biggest_val[row]);

      for (int i = row + 1; i < num; i++) {
        double other_big = std::fabs(matrix[i][row] / biggest_val[i]);
        if (other_big > big) {
          big = other_big;
          pivot_row = i;
        }
      }

      if (pivot_row != row) {
        for (int col = row; col < num; col++) std::swap(matrix[pivot_row][col], matrix[row][col]);
        std::swap(vec[pivot_row], vec[row]);
        std::swap(biggest_val[pivot_row], biggest_val[row]);
      }
    }

    // Solves matrix * solution = vector with Gaussian elimination (partial
    // pivoting) and back substitution. The solution replaces the contents of
    // `vec`.
    bool gauss_matrix(Matrix matrix, int num, Vector10 vec, std::string &err_mess) {
      err_mess.clear();
      double biggest_val[kNum];
      double solution_vector[kNum];

      // Step 1: find the largest value in each row, used to pivot the matrix.
      for (int i = 0; i < num; i++) {
        biggest_val[i] = std::fabs(matrix[i][0]);
        for (int j = 1; j < num; j++)
          if (std::fabs(matrix[i][j]) > biggest_val[i]) biggest_val[i] = std::fabs(matrix[i][j]);

        if (biggest_val[i] == 0.0) {
          err_mess = "Gauss_matrix: biggest val in row is zero";
          return false;
        }
      }

      // Step 2: convert the matrix into a triangular matrix.
      for (int i = 0; i < num - 1; i++) {
        gauss_pivot(matrix, num, vec, biggest_val, i);

        if (std::fabs(matrix[i][i] / biggest_val[i]) < kMatrixTol) {
          err_mess = "Gauss_matrix error: row has a too tiny value";
          return false;
        }

        for (int j = i + 1; j < num; j++) {
          double factor = matrix[j][i] / matrix[i][i];
          for (int k = i + 1; k < num; k++) matrix[j][k] -= factor * matrix[i][k];
          vec[j] -= factor * vec[i];
        }
      }

      if (std::fabs(matrix[num - 1][num - 1] / biggest_val[num - 1]) < kMatrixTol) {
        err_mess = "Gauss_matrix error: last row has a too tiny value";
        return false;
      }

      // Step 3: back substitution.
      solution_vector[num - 1] = vec[num - 1] / matrix[num - 1][num - 1];
      for (int i = num - 2; i >= 0; i--) {
        double sum = 0.0;
        for (int j = i + 1; j < num; j++) sum += matrix[i][j] * solution_vector[j];
        solution_vector[i] = (vec[i] - sum) / matrix[i][i];
      }

      // Step 4: return the solution in `vec`.
      for (int i = 0; i < num; i++) vec[i] = solution_vector[i];
      return true;
    }
  } // namespace

  bool calc_trans_cubic(const std::vector<SStar> &stars_reference,
                        const std::vector<SStar> &stars_distorted, TTrans &trans,
                        std::string &err_mess) {
    err_mess.clear();
    if (stars_reference.size() < 10) {
      // AT_MATCH_REQUIRE_CUBIC
      err_mess = "Calc_Trans_Cubic: Not enough equations.";
      return false;
    }

    // In the names below a '1' refers to a coordinate of the reference star,
    // which appears on both sides of the matrix equation, and a '2' to the
    // distorted star, which appears only on the left hand side.
    double sum = 0, sumx1 = 0, sumy1 = 0, sumx1sq = 0, sumx1y1 = 0, sumy1sq = 0;
    double sumx1cu = 0, sumx1sqy1 = 0, sumx1y1sq = 0, sumy1cu = 0;
    double sumx1qu = 0, sumx1cuy1 = 0, sumx1sqy1sq = 0, sumx1y1cu = 0, sumy1qu = 0;
    double sumx1pe = 0, sumx1quy1 = 0, sumx1cuy1sq = 0, sumx1sqy1cu = 0, sumx1y1qu = 0, sumy1pe = 0;
    double sumx1he = 0, sumx1pey1 = 0, sumx1quy1sq = 0, sumx1cuy1cu = 0, sumx1sqy1qu = 0,
        sumx1y1pe = 0, sumy1he = 0;
    double sumx2 = 0, sumx2x1 = 0, sumx2y1 = 0, sumx2x1sq = 0, sumx2x1y1 = 0, sumx2y1sq = 0,
        sumx2x1cu = 0, sumx2x1sqy1 = 0, sumx2x1y1sq = 0, sumx2y1cu = 0;
    double sumy2 = 0, sumy2x1 = 0, sumy2y1 = 0, sumy2x1sq = 0, sumy2x1y1 = 0, sumy2y1sq = 0,
        sumy2x1cu = 0, sumy2x1sqy1 = 0, sumy2x1y1sq = 0, sumy2y1cu = 0;

    // Take the minimum of the two arrays in case one list is longer.
    const size_t n = std::min(stars_reference.size(), stars_distorted.size());
    for (size_t i = 0; i < n; i++) {
      const double x1 = stars_reference[i].x;
      const double y1 = stars_reference[i].y;
      const double x2 = stars_distorted[i].x;
      const double y2 = stars_distorted[i].y;

      sumx2 += x2;
      sumx2x1 += x2 * x1;
      sumx2y1 += x2 * y1;
      sumx2x1sq += x2 * x1 * x1;
      sumx2x1y1 += x2 * x1 * y1;
      sumx2y1sq += x2 * y1 * y1;
      sumx2x1cu += x2 * x1 * x1 * x1;
      sumx2x1sqy1 += x2 * x1 * x1 * y1;
      sumx2x1y1sq += x2 * x1 * y1 * y1;
      sumx2y1cu += x2 * y1 * y1 * y1;

      sumy2 += y2;
      sumy2x1 += y2 * x1;
      sumy2y1 += y2 * y1;
      sumy2x1sq += y2 * x1 * x1;
      sumy2x1y1 += y2 * x1 * y1;
      sumy2y1sq += y2 * y1 * y1;
      sumy2x1cu += y2 * x1 * x1 * x1;
      sumy2x1sqy1 += y2 * x1 * x1 * y1;
      sumy2x1y1sq += y2 * x1 * y1 * y1;
      sumy2y1cu += y2 * y1 * y1 * y1;

      // Elements of the matrix.
      sum += 1.0;
      sumx1 += x1;
      sumy1 += y1;
      sumx1sq += x1 * x1;
      sumx1y1 += x1 * y1;
      sumy1sq += y1 * y1;
      sumx1cu += x1 * x1 * x1;
      sumx1sqy1 += x1 * x1 * y1;
      sumx1y1sq += x1 * y1 * y1;
      sumy1cu += y1 * y1 * y1;
      sumx1qu += x1 * x1 * x1 * x1;
      sumx1cuy1 += x1 * x1 * x1 * y1;
      sumx1sqy1sq += x1 * x1 * y1 * y1;
      sumx1y1cu += x1 * y1 * y1 * y1;
      sumy1qu += y1 * y1 * y1 * y1;
      sumx1pe += x1 * x1 * x1 * x1 * x1;
      sumx1quy1 += x1 * x1 * x1 * x1 * y1;
      sumx1cuy1sq += x1 * x1 * x1 * y1 * y1;
      sumx1sqy1cu += x1 * x1 * y1 * y1 * y1;
      sumx1y1qu += x1 * y1 * y1 * y1 * y1;
      sumy1pe += y1 * y1 * y1 * y1 * y1;
      sumx1he += x1 * x1 * x1 * x1 * x1 * x1;
      sumx1pey1 += x1 * x1 * x1 * x1 * x1 * y1;
      sumx1quy1sq += x1 * x1 * x1 * x1 * y1 * y1;
      sumx1cuy1cu += x1 * x1 * x1 * y1 * y1 * y1;
      sumx1sqy1qu += x1 * x1 * y1 * y1 * y1 * y1;
      sumx1y1pe += x1 * y1 * y1 * y1 * y1 * y1;
      sumy1he += y1 * y1 * y1 * y1 * y1 * y1;
    }

    Matrix matrix;
    Vector10 vec;

    // Fill the lower triangle and then transpose for the upper one. The same
    // normal-equations matrix is used for both axes.
    auto fill_matrix = [&]() {
      matrix[0][0] = sum;
      matrix[1][0] = sumx1;
      matrix[2][0] = sumy1;
      matrix[3][0] = sumx1sq;
      matrix[4][0] = sumx1y1;
      matrix[5][0] = sumy1sq;
      matrix[6][0] = sumx1cu;
      matrix[7][0] = sumx1sqy1;
      matrix[8][0] = sumx1y1sq;
      matrix[9][0] = sumy1cu;

      matrix[1][1] = sumx1sq;
      matrix[2][1] = sumx1y1;
      matrix[3][1] = sumx1cu;
      matrix[4][1] = sumx1sqy1;
      matrix[5][1] = sumx1y1sq;
      matrix[6][1] = sumx1qu;
      matrix[7][1] = sumx1cuy1;
      matrix[8][1] = sumx1sqy1sq;
      matrix[9][1] = sumx1y1cu;

      matrix[2][2] = sumy1sq;
      matrix[3][2] = sumx1sqy1;
      matrix[4][2] = sumx1y1sq;
      matrix[5][2] = sumy1cu;
      matrix[6][2] = sumx1cuy1;
      matrix[7][2] = sumx1sqy1sq;
      matrix[8][2] = sumx1y1cu;
      matrix[9][2] = sumy1qu;

      matrix[3][3] = sumx1qu;
      matrix[4][3] = sumx1cuy1;
      matrix[5][3] = sumx1sqy1sq;
      matrix[6][3] = sumx1pe;
      matrix[7][3] = sumx1quy1;
      matrix[8][3] = sumx1cuy1sq;
      matrix[9][3] = sumx1sqy1cu;

      matrix[4][4] = sumx1sqy1sq;
      matrix[5][4] = sumx1y1cu;
      matrix[6][4] = sumx1quy1;
      matrix[7][4] = sumx1cuy1sq;
      matrix[8][4] = sumx1sqy1cu;
      matrix[9][4] = sumx1y1qu;

      matrix[5][5] = sumy1qu;
      matrix[6][5] = sumx1cuy1sq;
      matrix[7][5] = sumx1sqy1cu;
      matrix[8][5] = sumx1y1qu;
      matrix[9][5] = sumy1pe;

      matrix[6][6] = sumx1he;
      matrix[7][6] = sumx1pey1;
      matrix[8][6] = sumx1quy1sq;
      matrix[9][6] = sumx1cuy1cu;

      matrix[7][7] = sumx1quy1sq;
      matrix[8][7] = sumx1cuy1cu;
      matrix[9][7] = sumx1sqy1qu;

      matrix[8][8] = sumx1sqy1qu;
      matrix[9][8] = sumx1y1pe;

      matrix[9][9] = sumy1he;

      for (int r = 0; r <= 8; r++)
        for (int c = r + 1; c <= 9; c++) matrix[r][c] = matrix[c][r];
    };

    // Coefficients x00, x10, x01, x20, x11, x02, x30, x21, x12, x03.
    fill_matrix();
    vec[0] = sumx2;
    vec[1] = sumx2x1;
    vec[2] = sumx2y1;
    vec[3] = sumx2x1sq;
    vec[4] = sumx2x1y1;
    vec[5] = sumx2y1sq;
    vec[6] = sumx2x1cu;
    vec[7] = sumx2x1sqy1;
    vec[8] = sumx2x1y1sq;
    vec[9] = sumx2y1cu;

    if (!gauss_matrix(matrix, kNum, vec, err_mess)) {
      err_mess += ", Calc_trans_cubic: can not solve for the x coefficients";
      return false;
    }

    trans.x00 = vec[0];
    trans.x10 = vec[1];
    trans.x01 = vec[2];
    trans.x20 = vec[3];
    trans.x11 = vec[4];
    trans.x02 = vec[5];
    trans.x30 = vec[6];
    trans.x21 = vec[7];
    trans.x12 = vec[8];
    trans.x03 = vec[9];

    // Coefficients y00, y10, y01, y20, y11, y02, y30, y21, y12, y03.
    fill_matrix();
    vec[0] = sumy2;
    vec[1] = sumy2x1;
    vec[2] = sumy2y1;
    vec[3] = sumy2x1sq;
    vec[4] = sumy2x1y1;
    vec[5] = sumy2y1sq;
    vec[6] = sumy2x1cu;
    vec[7] = sumy2x1sqy1;
    vec[8] = sumy2x1y1sq;
    vec[9] = sumy2y1cu;

    if (!gauss_matrix(matrix, kNum, vec, err_mess)) {
      err_mess += ", Calc_trans_cubic: can not solve for the y coefficients";
      return false;
    }

    trans.y00 = vec[0];
    trans.y10 = vec[1];
    trans.y01 = vec[2];
    trans.y20 = vec[3];
    trans.y11 = vec[4];
    trans.y02 = vec[5];
    trans.y30 = vec[6];
    trans.y21 = vec[7];
    trans.y12 = vec[8];
    trans.y03 = vec[9];
    return true;
  }
} // namespace astap
