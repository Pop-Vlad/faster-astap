// C++ port of ASTAP's plate solving algorithm.
//
// Original Pascal source (C) 2017-2026 by Han Kleijn, www.hnsky.org
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
// Ported from unit_command_line_solving.pas, unit_command_line_general.pas,
// unit_command_line_star_database.pas and unit_command_line_calc_trans_cubic.pas.

#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace astap {
  // Pascal: Timage_array = array of array of array of Single, indexed
  // [colour, y, x]. Stored here as one flat buffer; the fastest access order is
  // still row wise, exactly like the original.
  class ImageArray {
  public:
    ImageArray() = default;

    ImageArray(int colours, int height, int width) { resize(colours, height, width); }

    // Equivalent of SetLength(): new elements are zeroed.
    void resize(int colours, int height, int width) {
      colours_ = colours;
      height_ = height;
      width_ = width;
      data_.assign(static_cast<size_t>(colours) * height * width, 0.0f);
    }

    bool empty() const { return data_.empty(); }
    int colours() const { return colours_; }
    int height() const { return height_; }
    int width() const { return width_; }

    float *row(int c, int y) {
      return data_.data() + (static_cast<size_t>(c) * height_ + y) * width_;
    }

    const float *row(int c, int y) const {
      return data_.data() + (static_cast<size_t>(c) * height_ + y) * width_;
    }

    float &at(int c, int y, int x) { return row(c, y)[x]; }
    float at(int c, int y, int x) const { return row(c, y)[x]; }

  private:
    int colours_ = 0;
    int height_ = 0;
    int width_ = 0;
    std::vector<float> data_;
  };

  // Pascal: Tstar_list = array of array of double.
  // Row 0 = x, row 1 = y, row 2 = SNR for star lists.
  // Rows 0..5 = quad hash code, rows 6..7 = quad centre for quad lists.
  class RowList {
  public:
    RowList() = default;

    RowList(int rows, size_t count) { resize(rows, count); }

    // Equivalent of SetLength(list, rows, count): existing data is preserved,
    // new elements are zeroed.
    void resize(int rows, size_t count) {
      rows_.resize(static_cast<size_t>(rows));
      for (auto &r: rows_) r.resize(count, 0.0);
      count_ = count;
    }

    int rows() const { return static_cast<int>(rows_.size()); }
    size_t count() const { return count_; }
    bool empty() const { return count_ == 0; }

    double &operator()(int r, size_t i) { return rows_[static_cast<size_t>(r)][i]; }
    double operator()(int r, size_t i) const { return rows_[static_cast<size_t>(r)][i]; }

    std::vector<double> &operator[](int r) { return rows_[static_cast<size_t>(r)]; }
    const std::vector<double> &operator[](int r) const { return rows_[static_cast<size_t>(r)]; }

    void swap(RowList &o) {
      rows_.swap(o.rows_);
      std::swap(count_, o.count_);
    }

    double *data(int r) { return rows_[static_cast<size_t>(r)].data(); }
    const double *data(int r) const { return rows_[static_cast<size_t>(r)].data(); }

  private:
    std::vector<std::vector<double> > rows_;
    size_t count_ = 0;
  };

  // Pascal: solution_vector = array[0..2] of double.
  // Describes X := v[0]*x + v[1]*y + v[2]
  struct SolutionVector {
    double v[3] = {0.0, 0.0, 0.0};
    double &operator[](int i) { return v[i]; }
    double operator[](int i) const { return v[i]; }
  };

  // Pascal: Theader. Only the fields the solver actually needs are kept.
  struct Header {
    int bitpix = 16;
    int width = 0;
    int height = 0;
    int naxis = 0;
    int naxis3 = 1;

    double crpix1 = 0; // reference point X
    double crpix2 = 0; // reference point Y
    double cdelt1 = 0; // X pixel size (deg)
    double cdelt2 = 0; // Y pixel size (deg)
    double ra0 = 0; // centre / mount position, radians
    double dec0 = 0; // centre / mount position, radians

    double backgr = 0; // background value
    double star_level = 0; // star level for HFD ~2.25, above background
    double star_level2 = 0; // star level for HFD ~4.5, above background
    double noise_level = 0; // background noise level

    double crota1 = 0; // image rotation at centre in degrees
    double crota2 = 0;
    double cd1_1 = 0; // solution matrix
    double cd1_2 = 0;
    double cd2_1 = 0;
    double cd2_2 = 0;

    double xbinning = 1;
    double ybinning = 1;
    double xpixsz = 0; // pixel width in microns (after binning)
    double ypixsz = 0;

    double datamax_org = 0;

    // Raw 80 character header cards of the loaded file, used to write the .wcs
    // output. Kept in the original order.
    std::vector<std::string> cards;
  };

  // SIP (Simple Imaging Polynomial) distortion coefficients, third order.
  struct SipCoefficients {
    bool valid = false;
    int a_order = 0;
    int ap_order = 0;
    double a[4][4] = {}; // a[i][j] = A_i_j, pixel to sky
    double b[4][4] = {}; // b[i][j] = B_i_j, pixel to sky
    double ap[4][4] = {}; // ap[i][j] = AP_i_j, sky to pixel
    double bp[4][4] = {}; // bp[i][j] = BP_i_j, sky to pixel
  };

  // Exit statuses, identical to the ones of astap_cli.
  enum ErrorLevel : int {
    kErrNone = 0,
    kErrNoSolution = 1,
    kErrNotEnoughStars = 2,
    kErrImageRead = 16,
    kErrNoStarDatabase = 32,
    kErrStarDatabaseRead = 33,
  };
} // namespace astap
