// Step 1 and 2 of the ASTAP method: find background, noise and star level, then
// find the stars and their CCD x, y positions.
// Ported from unit_command_line_general.pas (get_hist, get_background, HFD) and
// unit_command_line_solving.pas (find_stars and friends).

#pragma once

#include <functional>
#include <vector>

#include "astap/types.h"

namespace astap {
  using LogFn = std::function<void(const std::string &)>;

  // Histogram of one colour plane, 0..65535. Mirrors the global `histogram` and
  // `his_mean` of the Pascal code.
  struct Histogram {
    std::vector<int> counts; // 65536 bins
    int mean = 0;

    Histogram() : counts(65536, 0) {
    }
  };

  // Pascal: get_hist. Builds the histogram of `colour` and its mean value.
  void get_hist(int colour, const ImageArray &img, Histogram &hist);

  // Pascal: get_background. Finds the background level, the noise level and the
  // two star levels used as detection thresholds.
  void get_background(int colour, const ImageArray &img, Header &head, bool calc_hist,
                      bool calc_noise_level, int max_stars, Histogram &hist,
                      const LogFn &log = nullptr);

  // Pascal: SigmaClippedMeanFromHistogram. Sigma clipped mean and standard
  // deviation of a rectangular sub section, used for the faint star retry.
  void sigma_clipped_mean_from_histogram(const ImageArray &img, int startx, int stopx, int starty,
                                         int stopy, int upperlimit, int max_iterations,
                                         double convergence_threshold, double &meanv, double &stdev);

  // Pascal: HFD. Calculates the half flux diameter, FWHM, SNR, flux and the
  // centre of gravity of the star around x1,y1. `rs` is the annulus radius.
  // All coordinates are zero based array positions.
  void hfd(const ImageArray &img, int x1, int y1, int rs, double &hfd1, double &star_fwhm,
           double &snr, double &flux, double &xc, double &yc);

  // Pascal: get_brightest_stars. Keeps only the brightest `nr_stars_required`
  // stars of the list, using a histogram of the SNR values.
  void get_brightest_stars(int nr_stars_required, double highest_snr, RowList &starlist);

  // Pascal: find_stars. Detects stars and returns them as [x, y, snr] rows.
  // The detection level is lowered in up to four retries until `max_stars` stars
  // are found.
  void find_stars(const ImageArray &img, const Header &head, double hfd_min, int max_stars,
                  RowList &starlist, double &mean_hfd, const LogFn &log = nullptr);
} // namespace astap
