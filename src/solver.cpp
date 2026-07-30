#include "astap/solver.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>

#include "astap/astro_math.h"
#include "astap/calc_trans_cubic.h"
#include "astap/quads.h"

namespace astap {
  void Solver::bin_mono_and_crop(int &binning, double crop, const ImageArray &img, ImageArray &img2,
                                 const LogFn &log) {
    const int nrcolors = img.colours();
    const int width5 = img.width();
    const int height5 = img.height();

    const int w = static_cast<int>(ptrunc(crop * width5 / binning)); // dimensions after bin and crop
    const int h = static_cast<int>(ptrunc(crop * height5 / binning));

    img2.resize(1, h, w);

    // crop 0.9 means shifting by 0.05*width.
    const int shiftX = static_cast<int>(pround(width5 * (1 - crop) / 2));
    const int shiftY = static_cast<int>(pround(height5 * (1 - crop) / 2));

    const double norm = static_cast<double>(nrcolors) * binning * binning;
    for (int fitsY = 0; fitsY < h; fitsY++) {
      float *out = img2.row(0, fitsY);
      const int y = shiftY + fitsY * binning;
      for (int fitsX = 0; fitsX < w; fitsX++) {
        const int x = shiftX + fitsX * binning;
        double val = 0;
        for (int k = 0; k < nrcolors; k++) // all colours, this makes the result mono
          for (int i = 0; i < binning; i++) {
            const float *row = img.row(k, y + i);
            for (int j = 0; j < binning; j++) val += row[x + j];
          }
        out[fitsX] = static_cast<float>(val / norm);
      }
    }
    (void) log;
  }

  void Solver::convert_mono(ImageArray &img) {
    if (img.colours() < 3) return;
    ImageArray tmp(1, img.height(), img.width());
    for (int y = 0; y < img.height(); y++)
      for (int x = 0; x < img.width(); x++)
        tmp.at(0, y, x) = (img.at(0, y, x) + img.at(1, y, x) + img.at(2, y, x)) / 3;
    img = std::move(tmp);
  }

  void Solver::apply_check_pattern_filter(ImageArray &img, const LogFn &log) {
    if (img.colours() > 1) {
      if (log) log("Skipping check pattern filter. This filter works only for raw OSC images!");
      return;
    }
    if (log) log("Applying check pattern filter.");

    const int h = img.height();
    const int w = img.width();

    double value[4] = {0, 0, 0, 0};
    long counter[4] = {0, 0, 0, 0};

    // Use one quarter of the image to find the factors. This also works a little
    // better when no dark-flat is subtracted and when the border is black.
    for (int fitsY = h / 4; fitsY <= (h * 3) / 4; fitsY++)
      for (int fitsX = w / 4; fitsX <= (w * 3) / 4; fitsX++) {
        const int idx = (fitsX % 2) + 2 * (fitsY % 2);
        value[idx] += img.at(0, fitsY, fitsX);
        counter[idx]++; // separate counters in case of odd dimensions
      }

    double maxval = 0;
    for (int i = 0; i < 4; i++) {
      if (counter[i] != 0) value[i] /= counter[i];
      maxval = std::max(maxval, value[i]);
    }

    // Normalise the bayer pattern pixels.
    for (int fitsY = 0; fitsY < h; fitsY++)
      for (int fitsX = 0; fitsX < w; fitsX++) {
        const int idx = (fitsX % 2) + 2 * (fitsY % 2);
        if (value[idx] != 0)
          img.at(0, fitsY, fitsX) = static_cast<float>(img.at(0, fitsY, fitsX) * maxval / value[idx]);
      }
  }

  void Solver::bin_and_find_stars(const ImageArray &img, Header &head, int binfactor, double cropping,
                                  double hfd_min, int max_stars, RowList &starlist,
                                  std::string &short_warning) {
    short_warning.clear();

    const int width2 = img.width();
    const int height2 = img.height();
    double mean_hfd = 0;

    if (binfactor > 1 || cropping < 1) {
      if (binfactor > 1)
        say("Creating grayscale x " + std::to_string(binfactor) +
            " binning image for solving/star alignment.");
      if (cropping != 1) say("Cropping image x " + float_to_str(cropping, 2));

      ImageArray img_binned;
      bin_mono_and_crop(binfactor, cropping, img, img_binned, log_);

      get_background(0, img_binned, head, true, true, settings_.max_stars, histogram_, log_);
      find_stars(img_binned, head, hfd_min, max_stars, starlist, mean_hfd,
                 settings_.show_log ? log_ : nullptr);

      if (height2 < 960) {
        // Dimensions should be about 1280x960 or better.
        short_warning = "Warning, remaining image dimensions too low! ";
        say("Warning, remaining image dimensions too low! Try to REDUCE OR REMOVE DOWNSAMPLING.");
      }

      // Correct the star positions for binning and cropping. For zero based
      // indexing a star of 2x2 pixels at [2.5,2.5] is after 2x2 binning at [1,1];
      // doubling to [2,2] shifts it by 0.5 pixel. So the correction is
      // x := (binfactor-1)*0.5 + binfactor*x.
      for (size_t i = 0; i < starlist.count(); i++) {
        starlist(0, i) = (binfactor - 1) * 0.5 + starlist(0, i) * binfactor +
                         (width2 * (1 - cropping) / 2);
        starlist(1, i) = (binfactor - 1) * 0.5 + starlist(1, i) * binfactor +
                         (height2 * (1 - cropping) / 2);
      }
    } else {
      if (height2 > 2500) {
        short_warning = "Warning, increase downsampling!! ";
        say("Info: DOWNSAMPLING IS RECOMMENDED FOR LARGE IMAGES.");
      } else if (height2 < 960) {
        short_warning = "Warning, small image dimensions!! ";
        say("Warning, small image dimensions!!");
      }

      if (img.colours() >= 3) {
        // colour image
        ImageArray mono = img;
        say("Converting to mono.");
        convert_mono(mono);
        get_background(0, mono, head, true, true, settings_.max_stars, histogram_, log_);
        find_stars(mono, head, hfd_min, max_stars, starlist, mean_hfd,
                   settings_.show_log ? log_ : nullptr);
      } else {
        get_background(0, img, head, true, true, settings_.max_stars, histogram_, log_);
        find_stars(img, head, hfd_min, max_stars, starlist, mean_hfd,
                   settings_.show_log ? log_ : nullptr);
      }
    }
  }

  int Solver::report_binning_astrometric(double height, double arcsec_per_px) const {
    int result = settings_.downsample;
    if (result <= 0) {
      // auto
      result = height > 2500 ? 2 : 1;
      // The pixel scale should be larger than 1"/px.
      result = std::max<long>(result, pround(1.5 / arcsec_per_px));
    }
    return std::min(16, result); // 16 max, too much anyhow
  }

  bool Solver::read_stars(double telescope_ra, double telescope_dec, double search_field,
                          int nrstars_required, RowList &starlist) {
    int nrstars = 0;
    double ra2 = 0; // define a value, the first record read could be a header record
    double dec2 = 0, b_v = 0;

    starlist.resize(2, static_cast<size_t>(nrstars_required));

    if (database_.database_type() != kDatabaseWideField) {
      // .1476 or .290 files
      // Assume the search field is at a crossing of four tiles. The field area is
      // split over the tiles, so with 500 stars required and a 8/15/20/57% split
      // it retrieves 8% x 500 stars from the first tile and so on. This works as
      // long as the star density is reasonably homogeneous.
      int area[4];
      double frac[4];
      database_.find_areas(telescope_ra, telescope_dec, search_field, area[0], area[1], area[2],
                           area[3], frac[0], frac[1], frac[2], frac[3]);

      double frac_sum = 0;
      for (int a = 0; a < 4; a++) {
        frac_sum += frac[a];
        if (area[a] == 0) continue;
        if (!database_.open_area(telescope_dec, area[a])) return false;
        // Prevent round up errors resulting in an out of range star list.
        const int nrstars_required2 =
            std::min<long>(nrstars_required, ptrunc(nrstars_required * frac_sum));
        while (nrstars < nrstars_required2 &&
               database_.read_star(telescope_ra, telescope_dec, search_field, ra2, dec2, mag2_, b_v)) {
          // Store the star CCD x,y position.
          equatorial_standard(telescope_ra, telescope_dec, ra2, dec2, 1,
                              starlist(0, static_cast<size_t>(nrstars)),
                              starlist(1, static_cast<size_t>(nrstars)));
          nrstars++;
        }
      }
    } else {
      // wide field database
      if (!database_.read_stars_wide_field(settings_.database_path)) return false;
      const std::vector<float> &stars = database_.wide_field_stars();
      size_t count = 0;
      while (nrstars < nrstars_required && count < stars.size() / 3) {
        ra2 = stars[count * 3 + 1]; // contains mag1, ra1, dec1, mag2, ra2, dec2 ...
        dec2 = stars[count * 3 + 2];
        double sep;
        // The angular separation is required for a large field of view around the
        // pole, the simple formulas cannot be used any more.
        ang_sep(ra2, dec2, telescope_ra, telescope_dec, sep);
        // The factor 2/sqrt(pi) adapts the circular search field to a square
        // surface, 0.9 is a fiddle factor for trees, houses and dark corners and
        // pi/2 is the limit of equatorial_standard.
        if (sep < search_field * 0.5 * 0.9 * (2 / std::sqrt(kPi)) && sep < kPi / 2) {
          equatorial_standard(telescope_ra, telescope_dec, ra2, dec2, 1,
                              starlist(0, static_cast<size_t>(nrstars)),
                              starlist(1, static_cast<size_t>(nrstars)));
          nrstars++;
        }
        count++;
      }
      if (count > 0) mag2_ = stars[(count - 1) * 3]; // faintest magnitude used
    }

    // Fix the array length in case fewer stars were found.
    if (nrstars < nrstars_required) starlist.resize(2, static_cast<size_t>(nrstars));
    return true;
  }

  bool Solver::add_sip_coefficients(const Header &head, double ra_database, double dec_database) {
    // 1) Solve the image with the first order solver.
    // 2) Get the x,y coordinates of the detected stars = "stars_measured".
    // 3) Get the x,y coordinates of the reference stars = "stars_reference".
    // 4) Shift stars_measured to the centre of the image, so [0,0] is at
    //    CRPIX1, CRPIX2.
    // 5) Convert the reference star coordinates to the same coordinate system:
    //    quad x,y -> ra,dec -> image position using the first order solution.
    // 6) Both lists now match except for distortion.
    // 7) Calc_Trans_Cubic(stars_measured, stars_reference) works for pixel to sky.
    // 8) Calc_Trans_Cubic(stars_reference, stars_measured) works for sky to pixel.

    const size_t len = match_.b_xrefpositions.size();
    if (len < 20) {
      say("Not enough quads for calculating SIP.");
      return false;
    }

    std::vector<SStar> stars_measured(len);
    std::vector<SStar> stars_reference(len);

    const double SIN_dec_ref = std::sin(head.dec0);
    const double COS_dec_ref = std::cos(head.dec0);

    for (size_t i = 0; i < len; i++) {
      // Position as seen from the centre at crpix1, crpix2, in the FITS range 1..width.
      stars_measured[i].x = 1 + match_.a_xy_positions(0, i) - head.crpix1;
      stars_measured[i].y = 1 + match_.a_xy_positions(1, i) - head.crpix2;

      double ra_t, dec_t;
      standard_equatorial(ra_database, dec_database, match_.b_xrefpositions[i],
                          match_.b_yrefpositions[i], 1, ra_t, dec_t);

      // Conversion (RA,DEC) -> x,y image in the FITS range 1..max.
      const double SIN_dec_t = std::sin(dec_t);
      const double COS_dec_t = std::cos(dec_t);
      const double delta_ra = ra_t - head.ra0;
      const double SIN_delta_ra = std::sin(delta_ra);
      const double COS_delta_ra = std::cos(delta_ra);

      const double H = SIN_dec_t * SIN_dec_ref + COS_dec_t * COS_dec_ref * COS_delta_ra;
      const double dRA = (COS_dec_t * SIN_delta_ra / H) * 180 / kPi;
      const double dDEC =
          ((SIN_dec_t * COS_dec_ref - COS_dec_t * SIN_dec_ref * COS_delta_ra) / H) * 180 / kPi;

      const double det = head.cd2_2 * head.cd1_1 - head.cd1_2 * head.cd2_1;
      stars_reference[i].x = -(head.cd1_2 * dDEC - head.cd2_2 * dRA) / det;
      stars_reference[i].y = +(head.cd1_1 * dDEC - head.cd2_1 * dRA) / det;
    }

    std::string err_mess;
    TTrans trans_sky_to_pixel;
    if (!calc_trans_cubic(stars_reference, stars_measured, trans_sky_to_pixel, err_mess)) {
      say(err_mess);
      return false;
    }

    // Sky to pixel coefficients.
    sip_.ap_order = 3;
    sip_.ap[0][0] = trans_sky_to_pixel.x00;
    sip_.ap[0][1] = trans_sky_to_pixel.x01;
    sip_.ap[0][2] = trans_sky_to_pixel.x02;
    sip_.ap[0][3] = trans_sky_to_pixel.x03;
    sip_.ap[1][0] = -1 + trans_sky_to_pixel.x10;
    sip_.ap[1][1] = trans_sky_to_pixel.x11;
    sip_.ap[1][2] = trans_sky_to_pixel.x12;
    sip_.ap[2][0] = trans_sky_to_pixel.x20;
    sip_.ap[2][1] = trans_sky_to_pixel.x21;
    sip_.ap[3][0] = trans_sky_to_pixel.x30;

    sip_.bp[0][0] = trans_sky_to_pixel.y00;
    sip_.bp[0][1] = -1 + trans_sky_to_pixel.y01;
    sip_.bp[0][2] = trans_sky_to_pixel.y02;
    sip_.bp[0][3] = trans_sky_to_pixel.y03;
    sip_.bp[1][0] = trans_sky_to_pixel.y10;
    sip_.bp[1][1] = trans_sky_to_pixel.y11;
    sip_.bp[1][2] = trans_sky_to_pixel.y12;
    sip_.bp[2][0] = trans_sky_to_pixel.y20;
    sip_.bp[2][1] = trans_sky_to_pixel.y21;
    sip_.bp[3][0] = trans_sky_to_pixel.y30;

    // Inverse transformation: swap the arrays. This works as long as the offset
    // is small, like in this situation.
    TTrans trans_pixel_to_sky;
    if (!calc_trans_cubic(stars_measured, stars_reference, trans_pixel_to_sky, err_mess)) {
      say(err_mess);
      return false;
    }

    // Pixel to sky coefficients. SIP definitions:
    // https://irsa.ipac.caltech.edu/data/SPITZER/docs/files/spitzer/shupeADASS.pdf
    sip_.a_order = 3;
    sip_.a[0][0] = trans_pixel_to_sky.x00;
    sip_.a[0][1] = trans_pixel_to_sky.x01;
    sip_.a[0][2] = trans_pixel_to_sky.x02;
    sip_.a[0][3] = trans_pixel_to_sky.x03;
    sip_.a[1][0] = -1 + trans_pixel_to_sky.x10;
    sip_.a[1][1] = trans_pixel_to_sky.x11;
    sip_.a[1][2] = trans_pixel_to_sky.x12;
    sip_.a[2][0] = trans_pixel_to_sky.x20;
    sip_.a[2][1] = trans_pixel_to_sky.x21;
    sip_.a[3][0] = trans_pixel_to_sky.x30;

    sip_.b[0][0] = trans_pixel_to_sky.y00;
    sip_.b[0][1] = -1 + trans_pixel_to_sky.y01;
    sip_.b[0][2] = trans_pixel_to_sky.y02;
    sip_.b[0][3] = trans_pixel_to_sky.y03;
    sip_.b[1][0] = trans_pixel_to_sky.y10;
    sip_.b[1][1] = trans_pixel_to_sky.y11;
    sip_.b[1][2] = trans_pixel_to_sky.y12;
    sip_.b[2][0] = trans_pixel_to_sky.y20;
    sip_.b[2][1] = trans_pixel_to_sky.y21;
    sip_.b[3][0] = trans_pixel_to_sky.y30;

    sip_.valid = true;
    return true;
  }

  bool Solver::solve(ImageArray img, Header &head) {
    const auto startTick = std::chrono::steady_clock::now();
    warning_str_.clear();
    errorlevel_ = kErrNone;
    sip_ = SipCoefficients();

    const double quad_tolerance = settings_.quad_tolerance;
    int max_stars = settings_.max_stars;

    const int width2 = img.width();
    const int height2 = img.height();
    head.width = width2;
    head.height = height2;

    double fov_org;
    if (!settings_.fov_specified && head.cdelt2 != 0)
      fov_org = apply_arctan(height2 * std::fabs(head.cdelt2)); // PI can give a negative cdelt2
    else
      fov_org = std::min(180.0, settings_.search_fov); // 180 max to prevent runtime errors later

    std::string db_warning;
    if (!database_.select(settings_.database_path, settings_.star_database, fov_org, &db_warning)) {
      say("Error, no star database found at " + settings_.database_path +
          " ! Download and install a star database.");
      errorlevel_ = kErrNoStarDatabase;
      return false;
    }
    warning_str_ += db_warning;
    say("Using star database " + database_.name());

    if (fov_org > 30 && database_.database_type() != kDatabaseWideField)
      warning_str_ += "Wide field image, use W08 database! ";
    else if (fov_org > 6 && database_.database_type() == kDatabase1476)
      warning_str_ += "Large FOV, use G05 database! ";
    if (!warning_str_.empty()) say(warning_str_);

    if (settings_.check_pattern_filter) apply_check_pattern_filter(img, log_);

    // The FOV should stay below the database tile dimensions, otherwise a tile
    // beyond the next one could be selected.
    double max_fov;
    if (database_.database_type() == kDatabase1476)
      max_fov = 5.142857143;
    else if (database_.database_type() == kDatabase290)
      max_fov = 9.53;
    else
      max_fov = 180;

    const double min_star_size_arcsec = settings_.min_star_size;
    const bool autoFOV = (fov_org == 0);

    // The database name encodes the star density, e.g. d50 is 5000 stars per
    // square degree. The old V17/G17/G18/H17/H18 databases have no such number.
    int database_density = 9999;
    {
      const std::string &n = database_.name();
      if (n.size() >= 3 && std::isdigit(static_cast<unsigned char>(n[1])) &&
          std::isdigit(static_cast<unsigned char>(n[2]))) {
        int d = std::stoi(n.substr(1, 2));
        database_density = (d == 17 || d == 18) ? 9999 : d * 100;
      }
    }

    bool solution = false;
    double fov2 = fov_org;
    double fov_min = 0;
    double cropping = 1;
    double centerX = 0, centerY = 0;
    double ra_database = 0, dec_database = 0;
    double crota2_rad = 0, cdelt1_arcsec = 0, cdelt2_arcsec = 0, flipped_image = 1;
    double ra_solved = 0, dec_solved = 0;
    double ra_seed = head.ra0, dec_seed = head.dec0;
    int nr_quads = 0;
    std::string warning_downsample;
    RowList starlist1, starlist2;

    do {
      // autoFOV loop
      ra_seed = head.ra0;
      dec_seed = head.dec0;

      if (autoFOV) {
        if (fov_org == 0) {
          if (database_.database_type() != kDatabaseWideField) {
            fov_org = 9.5;
            fov_min = 0.38;
          } else {
            fov_org = 90;
            fov_min = 12;
          }
        } else {
          fov_org = fov_org / 1.5;
        }
        say("Trying FOV: " + float_to_str(fov_org, 1));
      }

      if (fov_org > max_fov) {
        cropping = max_fov / fov_org;
        fov2 = max_fov; // temporary cropped image, adjust the FOV to adapt
      } else {
        cropping = 1;
        fov2 = fov_org;
      }

      // Limit in stars per square degree: limit = density * surface_full_image.
      const long limit = pround(database_density * sqr(fov2) * width2 / static_cast<double>(height2));
      if (limit < max_stars) {
        max_stars = static_cast<int>(limit); // reduce the number of stars to use
        say("Database limit for this FOV is " + std::to_string(max_stars) + " stars.");
      }

      const double arcsec_per_px = fov_org * 3600 / height2; // unbinned
      const int binning = report_binning_astrometric(height2 * cropping, arcsec_per_px);
      // Ignore hot pixels which are too small.
      const double hfd_min = std::max(0.8, min_star_size_arcsec / (binning * arcsec_per_px));

      // Do this on every repeat since hfd_min is adapted.
      bin_and_find_stars(img, head, binning, cropping, hfd_min, max_stars, starlist2,
                         warning_downsample);
      const int nrstars_image = static_cast<int>(starlist2.count());

      say("Search radius: " + float_to_str(settings_.radius_search, 1) + " degrees, " +
          "start position: " + prepare_ra(head.ra0, ": ") + ", " + prepare_dec(head.dec0, "d ") +
          ", image height: " + float_to_str(fov_org, 2) + " degrees, binning: " +
          std::to_string(binning) + "x" + std::to_string(binning) + ", image dimensions: " +
          std::to_string(width2) + "x" + std::to_string(height2) + ", quad tolerance: " +
          float_to_str(quad_tolerance, 4) + ", minimum star size: " +
          float_to_str(min_star_size_arcsec, 1) + "\", speed: " +
          (settings_.force_oversize ? "slow" : "normal"));

      // A little less, the square search field is based on the height only.
      const int nrstars_required =
          static_cast<int>(pround(nrstars_image * (height2 / static_cast<double>(width2))));

      solution = false;
      bool go_ahead = (nrstars_image >= 5); // bare minimum, should be more but let's try

      double oversize = 1;
      double radius = settings_.radius_search;
      int minimum_quads = 3;

      if (go_ahead) {
        // enough stars, let's find quads
        find_quads(nrstars_image, starlist2, match_.quad_star_distances2);
        nr_quads = static_cast<int>(match_.quad_star_distances2.count());
        go_ahead = nr_quads >= 3;

        // The step size is fixed. If few stars are detected the search window (so
        // the database read area) is increased up to 200%, guaranteeing that all
        // image quads are compared with the database quads while stepping through
        // the sky.
        if (nrstars_image < 35)
          oversize = 2; // square search window twice the image height
        else if (nrstars_image > 140) // at least 100 quads
          oversize = 1; // square search window equal to the image height
        else
          // Between 35 (=2) and 140 (=1); quads are area related so take the sqrt.
          oversize = 2 * std::sqrt(35.0 / nrstars_image);

        if (settings_.force_oversize) oversize = 2;
        // Limit the request to one tile, otherwise a tile beyond the next one
        // could be selected.
        oversize = std::min(oversize, max_fov / fov2);

        // Prevent false detections for star rich images. Three quads give the
        // three quad centre references and are the bare minimum.
        minimum_quads = 3 + nrstars_image / 140;
      } else {
        say("Only " + std::to_string(nrstars_image) + " stars found in image. Abort");
        errorlevel_ = kErrNotEnoughStars;
      }

      if (!go_ahead) continue;

      const double search_field = fov2 * (kPi / 180);
      double step_size = search_field; // fixed step size for the search spiral
      int max_distance;
      if (database_.database_type() == kDatabaseWideField) {
        // Make small steps for wide field images, this is much more reliable.
        step_size = step_size * 0.1;
        max_distance = static_cast<int>(pround(radius / (0.1 * fov2 + 0.00001))); // in steps
        say("Wide field, making small steps for reliable solving.");
      } else {
        max_distance = static_cast<int>(pround(radius / (fov2 + 0.00001))); // in steps
      }

      say(std::to_string(nrstars_image) + " stars, " + std::to_string(nr_quads) +
          " quads selected in the image. " +
          std::to_string(pround(nrstars_required * sqr(oversize))) + " database stars, " +
          std::to_string(pround(nr_quads * nrstars_required * sqr(oversize) / nrstars_image)) +
          " database quads required for the " + float_to_str(oversize * fov2, 2) +
          "d square search window. Step size " + float_to_str(fov2, 2) + "d. Oversize " +
          float_to_str(oversize, 2));

      int match_nr = 0;
      do {
        // maximum accuracy loop, a second solve after a match on a corner
        solution = false; // could be true from a single lock, super rare
        int count = 0; // search field counter
        double distance = 0;
        int spiral_x = 0, spiral_y = 0;
        int spiral_dx = 0; // first step size x
        int spiral_dy = -1; // first step size y

        do {
          // search in a squared spiral
          if (count != 0) {
            // Start with [0 0], then [1 0], [1 1], [0 1], [-1 1], [-1 0],
            // [-1 -1], [0 -1], [1 -1], [2 -1], [2 0] ...
            if ((spiral_x == spiral_y) || ((spiral_x < 0) && (spiral_x == -spiral_y)) ||
                ((spiral_x > 0) && (spiral_x == 1 - spiral_y))) {
              // Turning point: swap dx by negative dy and dy by dx.
              int spiral_t = spiral_dx;
              spiral_dx = -spiral_dy;
              spiral_dy = spiral_t;
            }
            spiral_x += spiral_dx; // walk through the square
            spiral_y += spiral_dy;
          }

          // Adapt the search field to the matrix position.
          dec_database = step_size * spiral_y + dec_seed;
          double flip = 0;
          if (dec_database > +kPi / 2) {
            // crossed the pole
            dec_database = kPi - dec_database;
            flip = kPi;
          } else if (dec_database < -kPi / 2) {
            dec_database = -kPi - dec_database;
            flip = kPi;
          }

          // Use the distance furthest away from the pole.
          const double extra = dec_database > 0 ? step_size / 2 : -step_size / 2;

          // The step is larger near the pole. This is an offset from zero.
          const double ra_database_offset = (step_size * spiral_x / std::cos(dec_database - extra));
          if (ra_database_offset <= +kPi / 2 + step_size / 2 &&
              ra_database_offset >= -kPi / 2) {
            // step_size for overlap
            // Add the offset to RA after the if statement, otherwise the search
            // would not be symmetrical.
            ra_database = fnmodulo(flip + ra_seed + ra_database_offset, 2 * kPi);

            double seperation;
            ang_sep(ra_database, dec_database, ra_seed, dec_seed, seperation);

            // Use only the circular area within the square area.
            if (seperation <= radius * kPi / 180 + step_size / 2) {
              // Report a new distance only once per square spiral, it costs CPU time.
              if (seperation * 180 / kPi > distance + fov_org) {
                distance = seperation * 180 / kPi;
                say("Search distance: " + std::to_string(pround(distance)) + "d");
              }

              // Read nrstars_required stars from the database. When the search
              // field is oversized the number of required stars increases with
              // the power of the oversize factor, so that the star density stays
              // the same as in the image to solve.
              double oversize2;
              if (match_nr == 0) {
                oversize2 = oversize;
              } else {
                // Use the full image for the second solve, but limit it to one
                // tile to prevent tile selection problems.
                oversize2 = std::min(max_fov / fov2,
                                     std::max(oversize, std::sqrt(sqr(width2 / static_cast<double>(
                                                                        height2)) + 1)));
              }
              const int nrstars_required2 =
                  static_cast<int>(pround(nrstars_required * oversize2 * oversize2));

              if (!read_stars(ra_database, dec_database, search_field * oversize2, nrstars_required2,
                              starlist1)) {
                say("Error, no star database found at " + settings_.database_path +
                    " ! Download and install a star database.");
                errorlevel_ = kErrStarDatabaseRead;
                return false;
              }

              if (match_nr == 1) {
                // A first solution was found: keep only the stars visible in the
                // image, so that stars outside the image boundaries are not used
                // to create database quads.
                size_t nstars_visible = 0;
                for (size_t i = 0; i < starlist1.count(); i++) {
                  double xi, yi;
                  rotate(crota2_rad, starlist1(0, i) / cdelt1_arcsec,
                         starlist1(1, i) / cdelt2_arcsec, xi, yi); // rotate to screen orientation
                  xi = centerX - xi;
                  yi = centerY - yi;
                  if (xi > 0 && xi < width2 && yi > 0 && yi < height2) {
                    starlist1(0, nstars_visible) = starlist1(0, i);
                    starlist1(1, nstars_visible) = starlist1(1, i);
                    nstars_visible++;
                  }
                }
                starlist1.resize(2, nstars_visible);
              }

              find_quads(nrstars_image, starlist1, match_.quad_star_distances1);

              if (settings_.show_log)
                say("Search " + std::to_string(count) + ", [" + std::to_string(spiral_x) + "," +
                    std::to_string(spiral_y) + "], position: " + prepare_ra(ra_database, ": ") +
                    prepare_dec(dec_database, "d ") + "\t down to magn " +
                    float_to_str(mag2_ / 10, 1) + "\t " + std::to_string(starlist1.count()) +
                    " database stars\t " + std::to_string(match_.quad_star_distances1.count()) +
                    " database quads to compare.");

              solution = find_offset_and_rotation(match_, minimum_quads, quad_tolerance,
                                                  settings_.show_log ? log_ : nullptr);
            } // within the search circle, otherwise the search is within a kind of square
          } // RA in range

          count++; // step further in the spiral
        } while (!(solution || spiral_x > max_distance));

        if (solution) {
          centerX = (width2 - 1) / 2.0; // centre of the image in the 0..width-1 range
          centerY = (height2 - 1) / 2.0; // centre of the image in the 0..height-1 range

          standard_equatorial(
            ra_database, dec_database,
            (match_.solution_vector_x[0] * centerX + match_.solution_vector_x[1] * centerY +
             match_.solution_vector_x[2]),
            (match_.solution_vector_y[0] * centerX + match_.solution_vector_y[1] * centerY +
             match_.solution_vector_y[2]),
            1 /*CCD scale*/, ra_solved, dec_solved);
          ra_seed = ra_solved; // re-seed the maximum accuracy second pass
          dec_seed = dec_solved;

          // Flipped, either vertically or horizontally but not both. Flipped both
          // ways equals a 180 degree rotation and is not seen as flipped.
          if (match_.solution_vector_x[0] * match_.solution_vector_y[1] -
              match_.solution_vector_x[1] * match_.solution_vector_y[0] >
              0)
            flipped_image = -1;
          else
            flipped_image = +1;

          // The position one pixel away in the direction of crpix2.
          double ra7, dec7;
          standard_equatorial(
            ra_database, dec_database,
            (match_.solution_vector_x[0] * centerX + match_.solution_vector_x[1] * (centerY + 1) +
             match_.solution_vector_x[2]),
            (match_.solution_vector_y[0] * centerX + match_.solution_vector_y[1] * (centerY + 1) +
             match_.solution_vector_y[2]),
            1, ra7, dec7);

          // Position angle of the image +Y axis measured at the image centre from
          // north eastwards, negated for the CROTA2 convention.
          crota2_rad = -position_angle(ra7, dec7, ra_solved, dec_solved);
          cdelt1_arcsec = flipped_image * std::sqrt(sqr(match_.solution_vector_x[0]) +
                                                    sqr(match_.solution_vector_x[1]));
          cdelt2_arcsec =
              std::sqrt(sqr(match_.solution_vector_y[0]) + sqr(match_.solution_vector_y[1]));

          match_nr++;
        } else {
          match_nr = 0; // should not happen for the second solve, but just in case
        }
        // After a match possible on a corner, do a second solve using the found
        // position for maximum accuracy, using all stars.
      } while (!(!solution || match_nr >= 2));

      // Loop for autoFOV from 9.5 down to 0.37 degrees.
    } while (!(!autoFOV || solution || fov2 <= fov_min));

    const auto elapsed = std::chrono::steady_clock::now() - startTick;
    solved_seconds_ =
        std::round(std::chrono::duration<double>(elapsed).count() * 10) / 10;

    if (!solution) {
      say("No solution found!  :(");
      warning_str_ += warning_downsample;
      if (errorlevel_ == kErrNone) errorlevel_ = kErrNoSolution;
      return false;
    }

    // Calculate the search offset before updating head.ra0, head.dec0.
    ang_sep(ra_seed, dec_seed, head.ra0, head.dec0, sep_search_);
    head.ra0 = ra_solved;
    head.dec0 = dec_solved;
    head.crpix1 = centerX + 1; // centre of the image in the FITS range 1..width
    head.crpix2 = centerY + 1;

    say(std::to_string(match_.nr_references) + " of " + std::to_string(match_.nr_references2) +
        " quads selected matching within " + float_to_str(quad_tolerance, 4) + " tolerance.");
    say("Solution[\"] x:=" + float_to_str(match_.solution_vector_x[0], 6) + "*x+ " +
        float_to_str(match_.solution_vector_x[1], 6) + "*y+ " +
        float_to_str(match_.solution_vector_x[2], 6) +
        ",  y:=" + float_to_str(match_.solution_vector_y[0], 6) + "*x+ " +
        float_to_str(match_.solution_vector_y[1], 6) + "*y+ " +
        float_to_str(match_.solution_vector_y[2], 6));

    // The position 1*flipped_image pixels in the direction of crpix1.
    double ra7, dec7;
    standard_equatorial(
      ra_database, dec_database,
      (match_.solution_vector_x[0] * (centerX + flipped_image) +
       match_.solution_vector_x[1] * centerY + match_.solution_vector_x[2]),
      (match_.solution_vector_y[0] * (centerX + flipped_image) +
       match_.solution_vector_y[1] * centerY + match_.solution_vector_y[2]),
      1, ra7, dec7);

    // Position angle of the image +X axis, converted from north-referenced to the
    // CROTA1 convention.
    double crota1_rad = kPi / 2 - position_angle(ra7, dec7, ra_solved, dec_solved);
    if (crota1_rad > kPi) crota1_rad -= 2 * kPi; // keep within the range -pi..+pi

    head.cdelt1 = cdelt1_arcsec / 3600; // convert from arc seconds to degrees
    head.cdelt2 = cdelt2_arcsec / 3600;

    head.cd1_1 = +head.cdelt1 * std::cos(crota1_rad);
    head.cd1_2 = -head.cdelt1 * std::sin(crota1_rad) * flipped_image;
    head.cd2_1 = +head.cdelt2 * std::sin(crota2_rad) * flipped_image;
    head.cd2_2 = +head.cdelt2 * std::cos(crota2_rad);

    head.crota2 = crota2_rad * 180 / kPi;
    head.crota1 = crota1_rad * 180 / kPi;

    std::string mount_info;
    if (ra_mount_ < 99) {
      // mount position known and specified
      // Map to the range -pi..+pi.
      double delta_ra_mount = fnmodulo(ra_mount_ - head.ra0, 2 * kPi);
      if (delta_ra_mount > kPi) delta_ra_mount -= 2 * kPi;
      const std::string ra_offset = distance_to_string(
        sep_search_, delta_ra_mount * std::cos((head.dec0 + dec_mount_) * 0.5));
      const std::string dec_offset = distance_to_string(sep_search_, dec_mount_ - head.dec0);
      mount_info = " Mount offset RA=" + ra_offset + ", DEC=" + dec_offset + ".";
    }

    say("Solution found: " + prepare_ra(head.ra0, ": ") + " " + prepare_dec(head.dec0, "d ") +
        "\nSolved in " + float_to_str(solved_seconds_, 1) + " sec. Offset was " +
        distance_to_string(sep_search_, sep_search_) + "." + mount_info +
        " Used stars down to magnitude: " + float_to_str(mag2_ / 10, 1));

    if (settings_.add_sip) add_sip_coefficients(head, ra_database, dec_database);

    const double vfov = apply_arctan(height2 * head.cdelt2);
    if (fov_org > 1.05 * vfov || fov_org < 0.95 * vfov) {
      std::string suggest_str = "Warning scale was inaccurate! Set FOV=" +
                                float_to_str(height2 * head.cdelt2, 2) +
                                "d, scale=" + float_to_str(head.cdelt2 * 3600, 1) + "\"";
      if (head.xpixsz != 0)
        suggest_str += ", FL=" +
            std::to_string(pround(180 / (kPi * 1000) * head.xpixsz / head.cdelt2)) + "mm";
      say(suggest_str);
      warning_str_ = suggest_str + warning_str_;
    }

    warning_str_ += warning_downsample;
    return true;
  }

  // ---------------------------------------------------------------------------
  // Header card helpers and output files
  // ---------------------------------------------------------------------------

  namespace {
    // FITS style value field: right aligned in positions 11..30.
    std::string format_value(const std::string &v) {
      std::string s = v;
      if (s.size() < 20) s = std::string(20 - s.size(), ' ') + s;
      return s;
    }

    size_t find_key(const std::vector<std::string> &cards, const std::string &key) {
      for (size_t i = 0; i < cards.size(); i++)
        if (cards[i].compare(0, key.size(), key) == 0) return i;
      return cards.size();
    }

    // Inserts before the END card, or appends when there is none.
    void insert_card(std::vector<std::string> &cards, const std::string &card) {
      size_t end = find_key(cards, "END ");
      if (end < cards.size())
        cards.insert(cards.begin() + static_cast<long>(end), card);
      else
        cards.push_back(card);
    }
  } // namespace

  void update_text(std::vector<std::string> &cards, const std::string &key,
                   const std::string &value_and_comment) {
    std::string card = key;
    card.resize(9, ' ');
    card += value_and_comment;
    card.resize(80, ' ');
    size_t pos = find_key(cards, key);
    if (pos < cards.size())
      cards[pos] = card;
    else
      insert_card(cards, card);
  }

  void update_float(std::vector<std::string> &cards, const std::string &key,
                    const std::string &comment, double x) {
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%20.10G", x);
    update_text(cards, key, std::string(buf) + comment);
  }

  void update_integer(std::vector<std::string> &cards, const std::string &key,
                      const std::string &comment, long x) {
    update_text(cards, key, format_value(std::to_string(x)) + comment);
  }

  void remove_key(std::vector<std::string> &cards, const std::string &key) {
    size_t pos = find_key(cards, key);
    if (pos < cards.size()) cards.erase(cards.begin() + static_cast<long>(pos));
  }

  void update_solution_cards(std::vector<std::string> &cards, const Header &head,
                             const SipCoefficients &sip, bool solution, const std::string &comment) {
    if (!solution) {
      update_text(cards, "PLTSOLVD=", "                   F / No plate solution found.   ");
      remove_key(cards, "COMMENT 7");
      return;
    }

    if (sip.valid) {
      update_text(cards, "CTYPE1  =",
                  "'RA---TAN-SIP'       / TAN (gnomic) projection + SIP distortions      ");
      update_text(cards, "CTYPE2  =",
                  "'DEC--TAN-SIP'       / TAN (gnomic) projection + SIP distortions      ");
    } else {
      update_text(cards, "CTYPE1  =",
                  "'RA---TAN'           / first parameter RA,    projection TANgential   ");
      update_text(cards, "CTYPE2  =",
                  "'DEC--TAN'           / second parameter DEC,  projection TANgential   ");
    }
    update_text(cards, "CUNIT1  =",
                "'deg     '           / Unit of coordinates                            ");
    update_text(cards, "EQUINOX =",
                "              2000.0 / Equinox of coordinates                         ");

    update_float(cards, "CRPIX1  =", " / X of reference pixel                           ", head.crpix1);
    update_float(cards, "CRPIX2  =", " / Y of reference pixel                           ", head.crpix2);
    update_float(cards, "CRVAL1  =", " / RA of reference pixel (deg)                    ",
                 head.ra0 * 180 / kPi);
    update_float(cards, "CRVAL2  =", " / DEC of reference pixel (deg)                   ",
                 head.dec0 * 180 / kPi);
    update_float(cards, "CDELT1  =", " / X pixel size (deg)                             ", head.cdelt1);
    update_float(cards, "CDELT2  =", " / Y pixel size (deg)                             ", head.cdelt2);
    update_float(cards, "CROTA1  =", " / Image twist of X axis        (deg)             ", head.crota1);
    update_float(cards, "CROTA2  =", " / Image twist of Y axis        (deg)             ", head.crota2);
    update_float(cards, "CD1_1   =", " / CD matrix to convert (x,y) to (Ra, Dec)        ", head.cd1_1);
    update_float(cards, "CD1_2   =", " / CD matrix to convert (x,y) to (Ra, Dec)        ", head.cd1_2);
    update_float(cards, "CD2_1   =", " / CD matrix to convert (x,y) to (Ra, Dec)        ", head.cd2_1);
    update_float(cards, "CD2_2   =", " / CD matrix to convert (x,y) to (Ra, Dec)        ", head.cd2_2);
    update_text(cards, "PLTSOLVD=",
                "                   T / Astrometric solved by the ASTAP C++ port.      ");
    if (!comment.empty()) update_text(cards, "COMMENT 7", comment);

    if (!sip.valid) return;

    const char *names[10] = {"0_0", "1_0", "0_1", "2_0", "1_1", "0_2", "3_0", "2_1", "1_2", "0_3"};
    const int ij[10][2] = {
      {0, 0}, {1, 0}, {0, 1}, {2, 0}, {1, 1},
      {0, 2}, {3, 0}, {2, 1}, {1, 2}, {0, 3}
    };
    const char *comment_sip = " / SIP coefficient                                ";

    update_integer(cards, "A_ORDER =", " / Polynomial order, axis 1. Pixel to sky.        ", 3);
    for (int k = 0; k < 10; k++) {
      std::string key = std::string("A_") + names[k];
      key.resize(8, ' ');
      update_float(cards, key + "=", comment_sip, sip.a[ij[k][0]][ij[k][1]]);
    }
    update_integer(cards, "B_ORDER =", " / Polynomial order, axis 2. Pixel to sky.        ", 3);
    for (int k = 0; k < 10; k++) {
      std::string key = std::string("B_") + names[k];
      key.resize(8, ' ');
      update_float(cards, key + "=", comment_sip, sip.b[ij[k][0]][ij[k][1]]);
    }
    update_integer(cards, "AP_ORDER=", " / Inv polynomial order, axis 1. Sky to pixel.    ", 3);
    for (int k = 0; k < 10; k++) {
      std::string key = std::string("AP_") + names[k];
      key.resize(8, ' ');
      update_float(cards, key + "=", comment_sip, sip.ap[ij[k][0]][ij[k][1]]);
    }
    update_integer(cards, "BP_ORDER=", " / Inv polynomial order, axis 2. Sky to pixel.    ", 3);
    for (int k = 0; k < 10; k++) {
      std::string key = std::string("BP_") + names[k];
      key.resize(8, ' ');
      update_float(cards, key + "=", comment_sip, sip.bp[ij[k][0]][ij[k][1]]);
    }
  }

  bool write_ini(const std::string &filename, bool solution, const Header &head,
                 const std::string &cmdline, int errorlevel, const std::string &warning) {
    std::ofstream f(filename);
    if (!f.is_open()) return false;

    auto e = [](double x) {
      char buf[48];
      std::snprintf(buf, sizeof(buf), "%.10E", x);
      return std::string(buf);
    };

    if (solution) {
      f << "PLTSOLVD=T\n";
      f << "CRPIX1=" << e(head.crpix1) << "\n"; // X of reference pixel
      f << "CRPIX2=" << e(head.crpix2) << "\n"; // Y of reference pixel
      f << "CRVAL1=" << e(head.ra0 * 180 / kPi) << "\n"; // RA of reference pixel [deg]
      f << "CRVAL2=" << e(head.dec0 * 180 / kPi) << "\n"; // DEC of reference pixel [deg]
      f << "CDELT1=" << e(head.cdelt1) << "\n"; // X pixel size [deg]
      f << "CDELT2=" << e(head.cdelt2) << "\n"; // Y pixel size [deg]
      f << "CROTA1=" << e(head.crota1) << "\n"; // image twist of the X axis [deg]
      f << "CROTA2=" << e(head.crota2) << "\n"; // image twist of the Y axis [deg]
      f << "CD1_1=" << e(head.cd1_1) << "\n";
      f << "CD1_2=" << e(head.cd1_2) << "\n";
      f << "CD2_1=" << e(head.cd2_1) << "\n";
      f << "CD2_2=" << e(head.cd2_2) << "\n";
    } else {
      f << "\n";
      f << "PLTSOLVD=F\n";
    }
    f << "CMDLINE=" << cmdline << "\n";

    switch (errorlevel) {
      case kErrNotEnoughStars: f << "ERROR=Not enough stars.\n";
        break;
      case kErrImageRead: f << "ERROR=Error reading image file.\n";
        break;
      case kErrNoStarDatabase: f << "ERROR=No star database found.\n";
        break;
      case kErrStarDatabaseRead: f << "ERROR=Error reading star database.\n";
        break;
      default: break;
    }
    if (!warning.empty()) f << "WARNING=" << warning << "\n";
    return true;
  }
} // namespace astap
