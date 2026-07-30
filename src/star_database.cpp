#include "astap/star_database.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <sys/stat.h>

#include "astap/astro_math.h"

namespace astap {
  namespace {
    bool file_exists(const std::string &name) {
      struct stat st;
      return ::stat(name.c_str(), &st) == 0;
    }

    std::string to_lower(std::string s) {
      for (char &c: s) c = static_cast<char>(::tolower(static_cast<unsigned char>(c)));
      return s;
    }

    // The .290 layout: 18 rings of constant declination, the number of areas per
    // ring is given below. The declinations follow from
    // arcsin(1-1/289), arcsin(1-(1+8)/289), arcsin(1-(1+8+16)/289) ...
    const double kDecBoundaries290[19] = {
      -90 * kPi / 180, -85.23224404 * kPi / 180, -75.66348756 * kPi / 180,
      -65.99286637 * kPi / 180, -56.14497387 * kPi / 180, -46.03163067 * kPi / 180,
      -35.54307745 * kPi / 180, -24.53348115 * kPi / 180, -12.79440589 * kPi / 180,
      0, 12.79440589 * kPi / 180, 24.53348115 * kPi / 180,
      35.54307745 * kPi / 180, 46.03163067 * kPi / 180, 56.14497387 * kPi / 180,
      65.99286637 * kPi / 180, 75.66348756 * kPi / 180, 85.23224404 * kPi / 180,
      90 * kPi / 180
    };
    const int kCells290[18] = {1, 4, 8, 12, 16, 20, 24, 28, 32, 32, 28, 24, 20, 16, 12, 8, 4, 1};

    // The .1476 layout: 36 rings of about 5 degrees (90/17.5), each divided in
    // cells whose minimum width in RA is about 5 degrees.
    const double kDecBoundaries1476[37] = {
      -90.00000000 * kPi / 180, -87.42857143 * kPi / 180, -82.28571429 * kPi / 180,
      -77.14285714 * kPi / 180, -72.00000000 * kPi / 180, -66.85714286 * kPi / 180,
      -61.71428571 * kPi / 180, -56.57142857 * kPi / 180, -51.42857143 * kPi / 180,
      -46.28571429 * kPi / 180, -41.14285714 * kPi / 180, -36.00000000 * kPi / 180,
      -30.85714286 * kPi / 180, -25.71428571 * kPi / 180, -20.57142857 * kPi / 180,
      -15.42857143 * kPi / 180, -10.28571429 * kPi / 180, -5.142857143 * kPi / 180,
      0.0,
      5.142857143 * kPi / 180, 10.28571429 * kPi / 180, 15.42857143 * kPi / 180,
      20.57142857 * kPi / 180, 25.71428571 * kPi / 180, 30.85714286 * kPi / 180,
      36.00000000 * kPi / 180, 41.14285714 * kPi / 180, 46.28571429 * kPi / 180,
      51.42857143 * kPi / 180, 56.57142857 * kPi / 180, 61.71428571 * kPi / 180,
      66.85714286 * kPi / 180, 72.00000000 * kPi / 180, 77.14285714 * kPi / 180,
      82.28571429 * kPi / 180, 87.42857143 * kPi / 180, 90.00000000 * kPi / 180
    };
    const int kCells1476[36] = {
      1, 3, 9, 15, 21, 27, 33, 38, 43, 48, 52, 56,
      60, 63, 65, 67, 68, 69, 69, 68, 67, 65, 63, 60,
      56, 52, 48, 43, 38, 33, 27, 21, 15, 9, 3, 1
    };
  } // namespace

  // The file names are '<ring><cell>.<ext>' with two digits each, e.g.
  // 0101.290 or 3601.1476. Generated here instead of storing the 1476 entry
  // table of the original.
  std::string StarDatabase::area_filename(int database_type, int area) {
    const int *cells;
    int nrings;
    const char *ext;
    if (database_type == kDatabase290) {
      cells = kCells290;
      nrings = 18;
      ext = "290";
    } else {
      cells = kCells1476;
      nrings = 36;
      ext = "1476";
    }

    int remaining = area; // 1 based
    for (int ring = 0; ring < nrings; ring++) {
      if (remaining <= cells[ring]) {
        char buf[32];
        std::snprintf(buf, sizeof(buf), "%02d%02d.%s", ring + 1, remaining, ext);
        return std::string(buf);
      }
      remaining -= cells[ring];
    }
    return std::string(); // out of range, should never happen
  }

  // For a ra, dec position find the star database area number and the
  // corresponding boundary distances N, E, W, S. The long if-chain of the
  // original is expressed here as a loop over the ring table; the arithmetic is
  // identical.
  void StarDatabase::area_and_boundaries(double ra1, double dec1, int &area_nr, double &spaceE,
                                         double &spaceW, double &spaceN, double &spaceS) const {
    const double *bounds;
    const int *cells;
    int nrings;
    int total;
    if (database_type_ == kDatabase290) {
      bounds = kDecBoundaries290;
      cells = kCells290;
      nrings = 18;
      total = 290;
    } else {
      bounds = kDecBoundaries1476;
      cells = kCells1476;
      nrings = 36;
      total = 1476;
    }

    const double cos_dec1 = std::cos(dec1);

    if (dec1 > bounds[nrings - 1]) {
      // celestial north pole area
      area_nr = total;
      spaceS = dec1 - bounds[nrings - 1];
      // Minimum, could go beyond the celestial pole so above +90 degrees.
      spaceN = bounds[nrings] - bounds[nrings - 1];
      spaceW = kPi * 2;
      spaceE = kPi * 2;
      return;
    }

    // Ring r (0 based) covers the declinations (bounds[r], bounds[r+1]] and holds
    // cells[r] areas. Walk from the second highest ring downwards, as the
    // original does, keeping a running count of the areas below the current ring.
    int cumulative = total - cells[nrings - 1];
    for (int ring = nrings - 2; ring >= 1; ring--) {
      cumulative -= cells[ring];
      if (dec1 > bounds[ring]) {
        const int n = cells[ring];
        const double rot = ra1 * n / (2 * kPi);
        area_nr = cumulative + 1 + static_cast<int>(ptrunc(rot));
        spaceS = dec1 - bounds[ring];
        spaceN = bounds[ring + 1] - dec1;
        spaceW = (kPi * 2 / n) * pfrac(rot) * cos_dec1; // RA decreases in the direction west
        spaceE = (kPi * 2 / n) * (1 - pfrac(rot)) * cos_dec1;
        return;
      }
    }

    area_nr = 1; // celestial south pole area
    // Minimum, could go beyond the celestial pole so below -90 degrees.
    spaceS = bounds[1] - bounds[0];
    spaceN = bounds[1] - dec1;
    spaceW = kPi * 2;
    spaceE = kPi * 2;
  }

  void StarDatabase::find_areas(double ra1, double dec1, double fov, int &area1, int &area2,
                                int &area3, int &area4, double &frac1, double &frac2, double &frac3,
                                double &frac4) const {
    // The FOV must be smaller than the database tile dimensions, otherwise a tile
    // beyond the next one could be selected. This crop should never be needed.
    if (database_type_ == kDatabase290)
      fov = std::min(fov, 9.53 * kPi / 180);
    else
      fov = std::min(fov, 5.142857 * kPi / 180);

    const double fov_half = fov / 2;

    double dec_cornerN = dec1 + fov_half; // above +pi/2 does not matter, it is all the pole area
    double dec_cornerS = dec1 - fov_half;
    double ra_cornerWN = ra1 - fov_half / std::cos(dec_cornerN);
    if (ra_cornerWN < 0) ra_cornerWN += 2 * kPi; // for the direction west the RA decreases
    double ra_cornerEN = ra1 + fov_half / std::cos(dec_cornerN);
    if (ra_cornerEN >= 2 * kPi) ra_cornerEN -= 2 * kPi;
    double ra_cornerWS = ra1 - fov_half / std::cos(dec_cornerS);
    if (ra_cornerWS < 0) ra_cornerWS += 2 * kPi;
    double ra_cornerES = ra1 + fov_half / std::cos(dec_cornerS);
    if (ra_cornerES >= 2 * kPi) ra_cornerES -= 2 * kPi;

    double spaceE, spaceW, spaceN, spaceS;

    // The fraction is the part of the image covered by that database area.
    area_and_boundaries(ra_cornerEN, dec_cornerN, area1, spaceE, spaceW, spaceN, spaceS);
    frac1 = std::min(spaceW, fov) * std::min(spaceS, fov) / (fov * fov);

    area_and_boundaries(ra_cornerWN, dec_cornerN, area2, spaceE, spaceW, spaceN, spaceS);
    frac2 = std::min(spaceE, fov) * std::min(spaceS, fov) / (fov * fov);

    area_and_boundaries(ra_cornerES, dec_cornerS, area3, spaceE, spaceW, spaceN, spaceS);
    frac3 = std::min(spaceW, fov) * std::min(spaceN, fov) / (fov * fov);

    area_and_boundaries(ra_cornerWS, dec_cornerS, area4, spaceE, spaceW, spaceN, spaceS);
    frac4 = std::min(spaceE, fov) * std::min(spaceN, fov) / (fov * fov);

    if (area2 == area1) {
      area2 = 0;
      frac2 = 0;
    }
    if (area3 == area1) {
      area3 = 0;
      frac3 = 0;
    }
    if (area4 == area1) {
      area4 = 0;
      frac4 = 0;
    }
    if (area3 == area2) {
      area3 = 0;
      frac3 = 0;
    }
    if (area4 == area2) {
      area4 = 0;
      frac4 = 0;
    }
    if (area4 == area3) {
      area4 = 0;
      frac4 = 0;
    }

    if (frac1 < 0.01) {
      area1 = 0;
      frac1 = 0;
    } // too small, ignore
    if (frac2 < 0.01) {
      area2 = 0;
      frac2 = 0;
    }
    if (frac3 < 0.01) {
      area3 = 0;
      frac3 = 0;
    }
    if (frac4 < 0.01) {
      area4 = 0;
      frac4 = 0;
    }
  }

  bool StarDatabase::select(const std::string &database_path, const std::string &database, double fov,
                            std::string *warning) {
    database_path_ = database_path;
    database_type_ = kDatabase1476;
    bool old_database = false;

    const std::string db = to_lower(database);
    const char typ = db.empty() ? 'a' : db[0];

    auto exists = [&](const std::string &name) { return file_exists(database_path_ + name); };
    auto pick = [&](const std::string &name, int type) {
      name_database_ = name;
      database_type_ = type;
    };

    if (typ != 'a') {
      // manual setting
      if (typ == 'w') {
        if (exists(db + "_0101.001")) {
          pick(db, kDatabaseWideField);
          return true;
        }
      } else if (typ == 'd' || typ == 'v' || typ == 'i' || typ == 'h') {
        // d80, v50, h18
        if (exists(db + "_0101.1476")) {
          pick(db, kDatabase1476);
          return true;
        }
      }
      // No else: there is both a v50 (.1476) and a v17 (.290).
      if (typ == 'v' || typ == 'g') {
        // v17, g18
        if (exists(db + "_0101.290")) {
          pick(db, kDatabase290);
          return true;
        }
      }
    }

    if (fov > 20 && exists("w08_0101.001")) {
      pick("w08", kDatabaseWideField);
      return true;
    }

    if (fov > 6 && exists("g05_0101.290")) {
      pick("g05", kDatabase290); // preference for G05 for a large FOV
    } else if (fov > 6 && exists("v05_0101.290")) {
      pick("v05", kDatabase290);
    } else if (fov > 6 && exists("v17_0101.290")) {
      pick("v17", kDatabase290);
      old_database = true;
    } else if (exists("d80_0101.1476")) {
      pick("d80", kDatabase1476); // for a tiny field of view
    } else if (exists("d50_0101.1476")) {
      pick("d50", kDatabase1476);
    } else if (exists("v50_0101.1476")) {
      pick("v50", kDatabase1476); // photometry database
    } else if (exists("d20_0101.1476")) {
      pick("d20", kDatabase1476);
    } else if (exists("d05_0101.1476")) {
      pick("d05", kDatabase1476);
    } else if (exists("g05_0101.290")) {
      pick("g05", kDatabase290);
    } else if (exists("v05_0101.290")) {
      pick("v05", kDatabase290);
    } else if (exists("h18_0101.1476")) {
      pick("h18", kDatabase1476);
      old_database = true;
    } else if (exists("g18_0101.290")) {
      pick("g18", kDatabase290);
      old_database = true;
    } else if (exists("h17_0101.1476")) {
      pick("h17", kDatabase1476);
      old_database = true;
    } else if (exists("v17_0101.290")) {
      pick("v17", kDatabase290);
      old_database = true;
    } else if (exists("g17_0101.290")) {
      pick("g17", kDatabase290);
      old_database = true;
    } else {
      return false;
    }

    if (old_database && warning) *warning = "Old database!";
    return true;
  }

  void StarDatabase::close() {
    if (file_open_) {
      file_.close();
      file_open_ = false;
    }
  }

  bool StarDatabase::open_area(double telescope_dec, int area) {
    cos_telescope_dec_ = std::cos(telescope_dec); // here to save CPU time

    if (area != old_area_ || !file_open_) {
      close();

      const std::string namefile = name_database_ + "_" + area_filename(database_type_, area);
      file_.open(database_path_ + namefile, std::ios::binary);
      if (!file_.is_open()) return false;
      file_open_ = true;

      cache_valid_pos_ = 0;

      char header[110];
      file_.read(header, 110); // 10x11 is 110 bytes
      if (file_.gcount() != 110) {
        close();
        return false;
      }
      record_size_ = header[109] == ' ' ? 11 : static_cast<unsigned char>(header[109]);
      database_version_ = static_cast<unsigned char>(header[108]);

      file_.seekg(0, std::ios::end);
      const std::streamoff size = file_.tellg();
      file_.seekg(110, std::ios::beg);
      cache_size_ = static_cast<size_t>(size) - 110;

      if (cache_size_ > cache_.size()) cache_.resize(cache_size_); // only grow, resizing costs time

      old_area_ = area;
    }
    // else: re-use the data already in the cache

    cache_position_ = 0;
    return true;
  }

  bool StarDatabase::read_star(double telescope_ra, double telescope_dec, double field_diameter,
                               double &ra, double &dec, double &mag, double &b_v) {
    constexpr size_t kBlockSize = 5 * 6 * 4 * 1024; // a multiple of the record sizes 5 and 6

    double delta_ra;
    do {
      bool header_record;
      do {
        if (cache_position_ >= cache_size_) return false; // end of file, no more data

        if (cache_position_ >= cache_valid_pos_) {
          size_t block_to_read = std::min(cache_size_, kBlockSize);
          block_to_read = std::min(block_to_read, cache_size_ - cache_valid_pos_);
          file_.read(reinterpret_cast<char *>(cache_.data() + cache_valid_pos_),
                     static_cast<std::streamsize>(block_to_read));
          cache_valid_pos_ += block_to_read;
        }

        const uint8_t *rec = cache_.data() + cache_position_;
        cache_position_ += static_cast<size_t>(record_size_);

        header_record = false;
        if (record_size_ == 5 || record_size_ == 6) {
          const uint32_t ra_raw = static_cast<uint32_t>(rec[0]) |
                                  (static_cast<uint32_t>(rec[1]) << 8) |
                                  (static_cast<uint32_t>(rec[2]) << 16);
          if (ra_raw == 0xFFFFFF) {
            // special magnitude record
            mag = static_cast<double>(rec[4]) - 16; // shifted by 16 to make Sirius positive
            // The magnitude stays valid until the next header record.
            dec9_storage_ = static_cast<int8_t>(static_cast<int>(rec[3]) - 128);
            header_record = true;
          } else {
            ra = ra_raw * (kPi * 2 / ((256.0 * 256 * 256) - 1));
            const int32_t dec_raw = (static_cast<int32_t>(dec9_storage_) << 16) +
                                    (static_cast<int32_t>(rec[4]) << 8) + static_cast<int32_t>(rec[3]);
            dec = dec_raw * (kPi * 0.5 / ((128.0 * 256 * 256) - 1));
            if (record_size_ == 6)
              b_v = database_version_ == 2
                      ? static_cast<double>(static_cast<int8_t>(rec[5]))
                      : -999; // old type, mark as not available
          }
        } else {
          return false; // only the record sizes 5 and 6 are supported, like astap_cli
        }
      } while (header_record);

      // Skip stars too far from the centre of the field of interest.
      delta_ra = std::fabs(ra - telescope_ra);
      if (delta_ra > kPi) delta_ra = kPi * 2 - delta_ra;
    } while (!(delta_ra * cos_telescope_dec_ < field_diameter / 2 &&
               std::fabs(dec - telescope_dec) < field_diameter / 2));
    return true;
  }

  bool StarDatabase::read_stars_wide_field(const std::string &database_path) {
    if (wide_database_ == name_database_ && !wide_field_stars_.empty()) return true;

    std::ifstream fs(database_path + name_database_ + "_0101.001", std::ios::binary);
    if (!fs.is_open()) return false;

    int32_t isize = 0;
    fs.read(reinterpret_cast<char *>(&isize), sizeof(isize));
    if (!fs || isize <= 0) return false;

    wide_field_stars_.resize(static_cast<size_t>(isize) * 3);
    fs.read(reinterpret_cast<char *>(wide_field_stars_.data()),
            static_cast<std::streamsize>(isize) * 3 * static_cast<std::streamsize>(sizeof(float)));
    if (!fs) {
      wide_field_stars_.clear();
      return false;
    }

    wide_database_ = name_database_; // remember which database is in memory
    return true;
  }
} // namespace astap
