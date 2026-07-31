#include "astap/star_detection.h"

#include <algorithm>
#include <cmath>
#include <string>

#include "astap/astro_math.h"
#include "astap/parallel.h"

namespace astap {
  namespace {
    void say(const LogFn &log, const std::string &s) {
      if (log) log(s);
    }
  } // namespace

  void get_hist(int colour, const ImageArray &img, Histogram &hist) {
    if (colour + 1 > img.colours()) colour = 0; // binned images are mono, use red

    std::fill(hist.counts.begin(), hist.counts.end(), 0);

    const int width5 = img.width();
    const int height5 = img.height();

    // If LibRaw was used, ignore the unused sensor areas up to 4.2% / 1.5%.
    const int offsetW = static_cast<int>(ptrunc(width5 * 0.042));
    const int offsetH = static_cast<int>(ptrunc(height5 * 0.015));

    const int first_row = offsetH;
    const int last_row = height5 - 1 - offsetH;
    if (last_row < first_row) {
      hist.mean = 0;
      return;
    }

    // Per-thread partial histograms, combined afterwards. The counts are
    // integers, and the pixel sum is accumulated as an integer as well (values
    // are below 65000 and there are at most a few billion of them, so it stays
    // exact), which makes the reduction independent of the thread count.
    const unsigned chunks = std::max(1u, range_chunks(static_cast<size_t>(last_row - first_row + 1)));
    std::vector<std::vector<int> > partial(chunks);
    std::vector<int64_t> sums(chunks, 0);
    std::vector<int64_t> counts(chunks, 0);

    parallel_ranges(static_cast<size_t>(first_row), static_cast<size_t>(last_row) + 1,
                    [&](size_t r0, size_t r1, unsigned t) {
                      std::vector<int> &h = partial[t];
                      h.assign(65536, 0);
                      int64_t sum = 0, cnt = 0;
                      for (size_t i = r0; i < r1; i++) {
                        const float *row = img.row(colour, static_cast<int>(i));
                        for (int j = offsetW; j <= width5 - 1 - offsetW; j++) {
                          int col = static_cast<int>(pround(row[j]));
                          // Ignore black overlap areas and bright stars.
                          if (col >= 1 && col < 65000) {
                            h[static_cast<size_t>(col)]++;
                            sum += col;
                            cnt++;
                          }
                        }
                      }
                      sums[t] = sum;
                      counts[t] = cnt;
                    });

    int64_t total_value = 0;
    int64_t count = 1; // prevent divide by zero
    for (unsigned t = 0; t < chunks; t++) {
      total_value += sums[t];
      count += counts[t];
      if (partial[t].empty()) continue;
      for (size_t v = 1; v < 65000; v++) hist.counts[v] += partial[t][v];
    }
    hist.mean = static_cast<int>(pround(static_cast<double>(total_value) / count));
  }

  void get_background(int colour, const ImageArray &img, Header &head, bool calc_hist,
                      bool calc_noise_level, int max_stars, Histogram &hist, const LogFn &log) {
    if (calc_hist) get_hist(colour, img, hist);

    head.backgr = img.at(0, 0, 0); // define something for images containing 0 or 65535 only

    // Find the peak in the histogram, which should be the average background.
    int pixels = 0;
    int max_range = hist.mean;
    if (max_range == 0) {
      head.backgr = 0; // empty colour
    } else {
      for (int i = 1; i <= max_range; i++) // ignore value 0 from oversize
        if (hist.counts[static_cast<size_t>(i)] > pixels) {
          pixels = hist.counts[static_cast<size_t>(i)];
          head.backgr = i;
        }
    }

    // Check the alternative mean value.
    if (hist.mean > 1.5 * head.backgr) {
      say(log, "Will use mean value " + std::to_string(hist.mean) +
               " as background rather than most common value " +
               std::to_string(static_cast<long>(pround(head.backgr))));
      head.backgr = hist.mean; // strange peak at a low value, use the mean instead
    }

    if (!calc_noise_level) return;

    // Calculate the noise level.
    const int width5 = img.width();
    const int height5 = img.height();
    int stepsize = static_cast<int>(pround(height5 / 71.0)); // about 71x71 = 5000 samples
    if (stepsize % 2 == 0) stepsize = stepsize + 1; // prevent problems with even raw OSC images
    if (stepsize <= 0) stepsize = 1;

    double sd = 99999;
    double sd_old;
    int iterations = 0;
    do {
      int counter = 0;
      sd_old = sd;
      sd = 0;
      for (int fitsX = 15; fitsX <= width5 - 1 - 15; fitsX += stepsize) {
        for (int fitsY = 15; fitsY <= height5 - 1 - 15; fitsY += stepsize) {
          double value = img.at(colour, fitsY, fitsX);
          // Not an outlier: noise should be symmetrical so it should be less than
          // twice the background.
          if (value < head.backgr * 2 && value != 0) {
            if (iterations == 0 || std::fabs(value - head.backgr) <= 3 * sd_old) {
              sd = sd + sqr(value - head.backgr);
              counter++;
            }
          }
        }
      }
      sd = counter != 0 ? std::sqrt(sd / counter) : 0.0;
      iterations++;
    } while (!((sd_old - sd) < 0.05 * sd || iterations >= 7));
    head.noise_level = sd;

    // Calculate the star levels.
    int i = (head.bitpix == 8 || head.bitpix == 24) ? 255 : 65001;
    head.star_level = 0;
    head.star_level2 = 0;
    // Empirical: number of pixels to test, produces about 700 stars at hfd=2.25
    // respectively hfd=4.5.
    double factor = 6.0 * max_stars;
    double factor2 = 24.0 * max_stars;
    double above = 0;
    while (head.star_level == 0 && i > head.backgr + 1 && i > 0) {
      i--;
      above = above + hist.counts[static_cast<size_t>(i)];
      if (above >= factor) head.star_level = i;
    }
    while (head.star_level2 == 0 && i > head.backgr + 1 && i > 0) {
      i--;
      above = above + hist.counts[static_cast<size_t>(i)];
      if (above >= factor2) head.star_level2 = i;
    }

    // Clip the calculated star levels:
    // 1) above 3.5*noise minimum, but also above the background value when there
    //    is no noise, so the minimum is 1
    // 2) below the saturated level, so subtract 1 for saturated images
    head.star_level = std::max(std::max(3.5 * sd, 1.0), head.star_level - head.backgr - 1);
    head.star_level2 = std::max(std::max(3.5 * sd, 1.0), head.star_level2 - head.backgr - 1);
  }

  // Pascal: get_hist2, the sub section histogram used by the sigma clipped mean.
  static void get_hist2(const ImageArray &img, int startx, int stopx, int starty, int stopy,
                        int upperlimit, std::vector<int> &histogram) {
    // Reused across calls: this runs for every tile of the faint star pass, and
    // the histogram spans 65500 bins.
    if (histogram.size() < static_cast<size_t>(upperlimit) + 1)
      histogram.resize(static_cast<size_t>(upperlimit) + 1);
    std::fill(histogram.begin(), histogram.begin() + upperlimit + 1, 0);
    for (int i = starty; i <= stopy; i++) {
      const float *row = img.row(0, i);
      for (int j = startx; j <= stopx; j++) {
        int col = static_cast<int>(pround(row[j]));
        if (col >= 1 && col < upperlimit) // ignore black overlap areas and bright stars
          histogram[static_cast<size_t>(col)]++;
      }
    }
  }

  void sigma_clipped_mean_from_histogram(const ImageArray &img, int startx, int stopx, int starty,
                                         int stopy, int upperlimit, int max_iterations,
                                         double convergence_threshold, double &meanv,
                                         double &stdev) {
    const double sigmaLow = 3.0; // standard values for astronomy
    const double sigmaHigh = 2.0;
    int iteration = 0;
    bool converged = false;
    meanv = 0;
    stdev = 0;

    int currentLowerLimit = 0;
    int currentUpperLimit = upperlimit;

    static thread_local std::vector<int> histogram;
    get_hist2(img, startx, stopx, starty, stopy, upperlimit, histogram);

    while (!converged && iteration < max_iterations) {
      double previousMean = meanv;
      double previousStdDev = stdev;

      double sum = 0.0;
      double sumSquares = 0.0;
      long totalCount = 0;

      for (int i = currentLowerLimit; i <= currentUpperLimit; i++) {
        if (i >= 0 && i < static_cast<int>(histogram.size())) {
          int val = histogram[static_cast<size_t>(i)];
          if (val > 0) {
            double binValue = i; // the bin index represents the pixel value
            sum += binValue * val;
            sumSquares += binValue * binValue * val;
            totalCount += val;
          }
        }
      }

      if (totalCount > 0) {
        meanv = sum / totalCount;
        if (totalCount > 1) {
          double variance = (sumSquares - (sum * sum) / totalCount) / (totalCount - 1);
          stdev = variance > 0 ? std::sqrt(variance) : 0.0;
        } else {
          stdev = 0.0;
        }
      } else {
        meanv = 0.0;
        stdev = 0.0;
        break; // no data left
      }

      iteration++;

      if (stdev > 0) {
        currentLowerLimit = std::max(0L, pround(meanv - sigmaLow * stdev));
        currentUpperLimit = std::min(static_cast<long>(upperlimit), pround(meanv + sigmaHigh * stdev));
      }

      if (iteration > 1 && std::fabs(meanv - previousMean) < convergence_threshold &&
          std::fabs(stdev - previousStdDev) < convergence_threshold)
        converged = true;
    }
  }

  void hfd(const ImageArray &img, int x1, int y1, int rs, double &hfd1, double &star_fwhm,
           double &snr, double &flux, double &xc, double &yc) {
    constexpr int kMaxRi = 74; // >= sqrt(sqr(rs+rs)+sqr(rs+rs))+1+2 assuming rs<=50
    int distance_histogram[kMaxRi + 1];
    double background[1001]; // 3*(2*PI*(50+3)) assuming rs<=50

    const int width2 = img.width();
    const int height2 = img.height();

    hfd1 = 999;
    snr = 0;
    star_fwhm = 0;
    flux = 0;
    xc = x1;
    yc = y1;

    // Calculates an image pixel value on sub pixel level.
    auto value_subpixel = [&](double sx, double sy) -> double {
      int x_trunc = static_cast<int>(ptrunc(sx));
      int y_trunc = static_cast<int>(ptrunc(sy));
      if (x_trunc <= 0 || x_trunc >= (width2 - 2) || y_trunc <= 0 || y_trunc >= (height2 - 2))
        return 0;
      double x_frac = pfrac(sx);
      double y_frac = pfrac(sy);
      double r = img.at(0, y_trunc, x_trunc) * (1 - x_frac) * (1 - y_frac); // left top
      r += img.at(0, y_trunc, x_trunc + 1) * (x_frac) * (1 - y_frac); // right top
      r += img.at(0, y_trunc + 1, x_trunc) * (1 - x_frac) * (y_frac); // left bottom
      r += img.at(0, y_trunc + 1, x_trunc + 1) * (x_frac) * (y_frac); // right bottom
      return r;
    };

    // rs should be <= 50 to prevent runtime errors.
    const int r1_square = rs * rs;
    const int r2 = rs + 1; // the annulus width is 1
    const int r2_square = r2 * r2;

    if (x1 - r2 <= 0 || x1 + r2 >= width2 - 1 || y1 - r2 <= 0 || y1 + r2 >= height2 - 1) return;

    double valmax = 0;

    // Median of the annulus is the local background.
    int counter = 0;
    for (int i = -r2; i <= r2; i++)
      for (int j = -r2; j <= r2; j++) {
        int distance = i * i + j * j; // sqr(distance) is faster than applying sqrt
        if (distance > r1_square && distance <= r2_square) {
          background[counter] = img.at(0, y1 + j, x1 + i);
          counter++;
        }
      }

    // Sorted and rewritten in place, as the original does: this runs once per
    // star candidate, so a heap allocation here is expensive.
    double star_bg = smedian(background, static_cast<size_t>(counter));
    for (int i = 0; i < counter; i++) background[i] = std::fabs(background[i] - star_bg);
    double mad_bg = smedian(background, static_cast<size_t>(counter)); // median absolute deviation
    // Conversion from MAD to SD for a normal distribution.
    double sd_bg = mad_bg * 1.4826;
    // Add some value for images with a zero noise background, otherwise the
    // background could be seen as a star.
    sd_bg = std::max(sd_bg, 1.0);

    bool boxed = false;
    int signal_counter = 0;
    do {
      // reduce the square annulus radius until symmetry, this removes nearby stars
      double SumVal = 0, SumValX = 0, SumValY = 0;
      signal_counter = 0;

      for (int i = -rs; i <= rs; i++)
        for (int j = -rs; j <= rs; j++) {
          double val = img.at(0, y1 + j, x1 + i) - star_bg;
          if (val > 3.0 * sd_bg) {
            SumVal += val;
            SumValX += val * i;
            SumValY += val * j;
            signal_counter++; // how many pixels are illuminated
          }
        }
      if (SumVal <= 12 * sd_bg) return; // no star found, too noisy, exit with hfd=999

      double Xg = SumValX / SumVal;
      double Yg = SumValY / SumVal;
      xc = x1 + Xg;
      yc = y1 + Yg;
      // centre of gravity found

      if (xc - rs < 0 || xc + rs > width2 - 1 || yc - rs < 0 || yc + rs > height2 - 1) return;

      // Are 2 out of 9 of the pixels inside the box illuminated? In general this
      // works better for solving than the ovality measurement used in the past.
      boxed = signal_counter >= (2.0 / 9.0) * sqr(rs + rs + 1);

      if (!boxed) {
        if (rs > 4)
          rs -= 2; // try a smaller window to exclude nearby stars
        else
          rs -= 1;
      }

      if (signal_counter <= 1) return; // one hot pixel
    } while (!(boxed || rs <= 1));

    rs += 2; // add some space

    // Build the signal histogram from the centre of gravity.
    for (int i = 0; i <= rs; i++) distance_histogram[i] = 0;
    for (int i = -rs; i <= rs; i++) {
      for (int j = -rs; j <= rs; j++) {
        int distance = static_cast<int>(pround(std::sqrt(static_cast<double>(i * i + j * j))));
        if (distance <= rs) {
          double val = value_subpixel(xc + i, yc + j) - star_bg;
          if (val > 3.0 * sd_bg) {
            // 3 * sd should be signal
            distance_histogram[distance]++;
            if (val > valmax) valmax = val; // record the peak value of the star
          }
        }
      }
    }

    int r_aperture = -1;
    int distance_top_value = 0;
    bool HistStart = false;
    int illuminated_pixels = 0;
    do {
      r_aperture++;
      illuminated_pixels += distance_histogram[r_aperture];
      // Continue until we find a value > 0: the centre of a defocused star image
      // can be black when the telescope has a central obstruction.
      if (distance_histogram[r_aperture] > 0) HistStart = true;
      if (distance_top_value < distance_histogram[r_aperture])
        distance_top_value = distance_histogram[r_aperture];
    } while (!(r_aperture >= rs ||
               (HistStart && distance_histogram[r_aperture] <= 0.1 * distance_top_value)));
    if (r_aperture >= rs) return; // star is equal to or larger than the box, abort

    if (r_aperture > 2 && illuminated_pixels < 0.35 * sqr(r_aperture + r_aperture - 2))
      return; // not a star disk but stars, abort with hfd 999

    // Get the HFD using the approximation that the HFD line divides the star in
    // equal portions of gravity.
    double SumVal = 0, SumValR = 0, pixel_counter = 0;
    for (int i = -r_aperture; i <= r_aperture; i++)
      for (int j = -r_aperture; j <= r_aperture; j++) {
        double val = value_subpixel(xc + i, yc + j) - star_bg;
        double r = std::sqrt(static_cast<double>(i * i + j * j)); // distance from the centre
        SumVal += val;
        // Method Kazuhisa Miyashita. Note the HFD is calculated over a square
        // area, which works more accurately than over a round area.
        SumValR += val * r;
        if (val >= valmax * 0.5) pixel_counter += 1; // pixels above half maximum
      }
    flux = std::max(SumVal, 0.00001); // prevent dividing by zero or negative values
    hfd1 = 2 * SumValR / flux;
    hfd1 = std::max(0.7, hfd1);

    // From the surface (counting pixels above half max) the diameter equals FWHM.
    star_fwhm = 2 * std::sqrt(pixel_counter / kPi);

    // snr := flux/sqrt(flux + r*r*pi*sd^2), valid for both shot noise limited
    // (bright stars) and sky background limited situations.
    snr = flux / std::sqrt(flux + sqr(r_aperture) * kPi * sqr(sd_bg));
  }

  void get_brightest_stars(int nr_stars_required, double highest_snr, RowList &starlist) {
    constexpr int kRange = 199;
    int snr_histogram[kRange + 1] = {0};

    double sqrtRange = std::sqrt(highest_snr);
    if (sqrtRange <= 0) return;

    for (size_t i = 0; i < starlist.count(); i++) {
      // Stretch the lower part with many similar stars by applying sqrt. This is
      // much faster than using ln().
      int snr_scaled = static_cast<int>(ptrunc(std::sqrt(starlist(2, i)) * kRange / sqrtRange));
      snr_scaled = std::max(0, std::min(kRange, snr_scaled));
      snr_histogram[snr_scaled]++;
    }

    long count = 0;
    int i = kRange;
    do {
      i--;
      count += snr_histogram[i];
    } while (!(i <= 0 || count >= nr_stars_required));

    double overshoot_correction =
        snr_histogram[i] != 0
          ? (count - nr_stars_required) / static_cast<double>(snr_histogram[i])
          : 0.0; // linear overshoot correction
    double snr_required = sqr(sqrtRange * (i + overshoot_correction) / kRange); // back from sqrt

    size_t kept = 0;
    size_t nrstars = starlist.count();
    for (size_t k = 0; k < nrstars; k++)
      if (starlist(2, k) >= snr_required) {
        // preserve the brightest stars
        starlist(0, kept) = starlist(0, k); // overwrite in the same array
        starlist(1, kept) = starlist(1, k);
        starlist(2, kept) = starlist(2, k);
        kept++;
      }
    starlist.resize(3, kept);
  }

  void find_stars(const ImageArray &img, const Header &head, double hfd_min, int max_stars,
                  RowList &starlist, double &mean_hfd, const LogFn &log) {
    constexpr size_t kBufferSize = 5000;
    constexpr int kRasterSteps = 12;

    const int width2 = img.width();
    const int height2 = img.height();

    starlist.resize(3, kBufferSize);

    // img_sa marks the areas already occupied by a detected star. The `retries`
    // value is stored so that the array never has to be cleared.
    ImageArray img_sa(1, height2, width2);

    double backgr = head.backgr;
    double noise_lev = head.noise_level;
    int retries = 4; // try up to four times to get enough stars from the image

    size_t nrstars = 0;
    double highest_snr = 0;
    double detection_level = 0;

    auto find_stars_routine = [&](int startX, int endX, int startY, int endY) {
      // Hoisted absolute thresholds, saves an operation per pixel.
      const float level_star = static_cast<float>(backgr + detection_level);
      const float level_cross = static_cast<float>(backgr + 4 * noise_lev);
      for (int fitsY = startY; fitsY <= endY; fitsY++) {
        const float *rowC = img.row(0, fitsY); // current row
        const float *rowM = img.row(0, fitsY - 1); // row above
        const float *rowP = img.row(0, fitsY + 1); // row below
        float *rowSA = img_sa.row(0, fitsY); // star area marking row
        for (int fitsX = startX; fitsX <= endX; fitsX++) {
          if (rowSA[fitsX] != retries /*star free area for this retry*/ &&
              rowC[fitsX] > level_star /*star*/) {
            int starpixels = 0;
            // Inspect in a cross around it, should be above 4*noise level.
            if (rowC[fitsX - 1] > level_cross) starpixels++;
            if (rowC[fitsX + 1] > level_cross) starpixels++;
            if (rowM[fitsX] > level_cross) starpixels++;
            if (rowP[fitsX] > level_cross) starpixels++;

            if (starpixels >= 2) {
              // at least 3 illuminated pixels, not a hot pixel
              double hfd1, star_fwhm, snr, flux, xc, yc;
              hfd(img, fitsX, fitsY, 14 /*annulus radius*/, hfd1, star_fwhm, snr, flux, xc, yc);
              if (hfd1 <= 30 && snr > 10 && hfd1 > hfd_min /*0.8 is two pixels minimum*/ &&
                  img_sa.at(0, static_cast<int>(pround(yc)), static_cast<int>(pround(xc))) !=
                  retries) {
                mean_hfd += hfd1; // sum up to calculate the mean HFD later

                // Mark the whole circular star area as occupied to prevent double
                // detections. A value between 2.5*hfd and 3.5*hfd gives the same
                // performance; in practice a star PSF has larger wings than a
                // Gaussian predicts.
                int radius = static_cast<int>(pround(3.0 * hfd1));
                int sqr_radius = radius * radius;
                int xci = static_cast<int>(pround(xc));
                int yci = static_cast<int>(pround(yc));
                for (int n = -radius; n <= radius; n++)
                  for (int m = -radius; m <= radius; m++) {
                    int j = n + yci;
                    int i = m + xci;
                    if (j >= 0 && i >= 0 && j < height2 && i < width2 &&
                        (m * m + n * n) <= sqr_radius)
                      img_sa.at(0, j, i) = static_cast<float>(retries);
                  }

                nrstars++;
                if (nrstars >= starlist.count()) starlist.resize(3, nrstars + kBufferSize);
                starlist(0, nrstars - 1) = xc; // store the star position
                starlist(1, nrstars - 1) = yc;
                starlist(2, nrstars - 1) = snr; // store the SNR

                if (snr > highest_snr) highest_snr = snr;
              }
            }
          }
        }
      }
    };

    do {
      mean_hfd = 0;
      highest_snr = 0;
      nrstars = 0;

      if (retries == 4) {
        if (head.star_level > 30 * noise_lev) {
          // stars are dominant
          detection_level = head.star_level;
          find_stars_routine(1, width2 - 1 - 1, 1, height2 - 1 - 1);
        } else {
          retries = 3; // skip
        }
      }
      if (retries == 3) {
        if (head.star_level2 > 30 * noise_lev) {
          // stars are dominant
          detection_level = head.star_level2;
          find_stars_routine(1, width2 - 1 - 1, 1, height2 - 1 - 1);
        } else {
          retries = 2; // skip
        }
      }
      if (retries == 2) {
        detection_level = 30 * noise_lev;
        find_stars_routine(1, width2 - 1 - 1, 1, height2 - 1 - 1);
      }
      if (retries == 1) {
        // Last try to find faint stars: divide the image in sections and find the
        // background and noise level for each section.
        int stepsX, stepsY;
        if (height2 < width2) {
          stepsX = kRasterSteps;
          stepsY = static_cast<int>(pround(kRasterSteps * height2 / static_cast<double>(width2)));
        } else {
          stepsY = kRasterSteps;
          stepsX = static_cast<int>(pround(kRasterSteps * width2 / static_cast<double>(height2)));
        }

        for (int yy = 0; yy <= stepsY; yy++)
          for (int xx = 0; xx <= stepsX; xx++) {
            int startX = 1 + static_cast<int>(pround(width2 * xx / static_cast<double>(stepsX + 1)));
            int endX = std::min<long>(width2 - 1 - 1,
                                      pround(width2 * (xx + 1) / static_cast<double>(stepsX + 1)));
            int startY = 1 + static_cast<int>(pround(height2 * yy / static_cast<double>(stepsY + 1)));
            int endY = std::min<long>(height2 - 1 - 1,
                                      pround(height2 * (yy + 1) / static_cast<double>(stepsY + 1)));

            sigma_clipped_mean_from_histogram(
              img, startX, endX, startY, endY,
              std::max<long>(65500, ptrunc(head.backgr * 2)), 6, 0.1, backgr, noise_lev);
            detection_level = 7 * noise_lev;
            find_stars_routine(startX, endX, startY, endY);
          }
      }

      if (log)
        say(log, std::to_string(nrstars) + " stars found of the requested " +
                 std::to_string(max_stars) + ". Background value is " +
                 std::to_string(static_cast<long>(pround(head.backgr))) +
                 ". Detection level used " +
                 std::to_string(static_cast<long>(pround(detection_level))) +
                 " above background. Star level is " +
                 std::to_string(static_cast<long>(pround(head.star_level))) +
                 " above background. Noise level is " + float_to_str(head.noise_level, 0));
      retries--; // try again with a lower detection level
    } while (!(nrstars >= static_cast<size_t>(max_stars) || retries <= 0));

    starlist.resize(3, nrstars);

    if (nrstars > static_cast<size_t>(max_stars)) {
      // reduce the number of stars if too high
      say(log, "Selecting the " + std::to_string(max_stars) + " brightest stars only.");
      get_brightest_stars(max_stars, highest_snr, starlist);
    }

    if (nrstars > 0) mean_hfd = mean_hfd / nrstars;
  }
} // namespace astap
