#include "astap/star_database.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <list>
#include <mutex>
#include <unordered_map>

#include "astap/astro_math.h"

#ifdef _WIN32
#include <windows.h>
#elif defined(__APPLE__)
#include <mach-o/dyld.h>
#else
#include <unistd.h>
#endif

namespace astap {
  // Exactly the scale factors the record decoder applied.
  static constexpr double kRaScale = kPi * 2 / ((256.0 * 256 * 256) - 1);
  static constexpr double kDecScale = kPi * 0.5 / ((128.0 * 256 * 256) - 1);

  std::string executable_directory() {
    std::string path;
#ifdef _WIN32
    char buf[4096];
    const DWORD n = GetModuleFileNameA(nullptr, buf, sizeof(buf));
    if (n) path.assign(buf, n);
#elif defined(__APPLE__)
    char buf[4096];
    uint32_t size = sizeof(buf);
    if (_NSGetExecutablePath(buf, &size) == 0) path = buf;
#else
    char buf[4096];
    const ssize_t n = readlink("/proc/self/exe", buf, sizeof(buf) - 1);
    if (n > 0) path.assign(buf, static_cast<size_t>(n));
#endif
    if (path.empty()) return path;
    const size_t slash = path.find_last_of("/\\");
    return slash == std::string::npos ? std::string() : path.substr(0, slash + 1);
  }

  std::vector<std::string> default_database_directories() {
    std::vector<std::string> dirs;
    const std::string here = executable_directory();
    if (!here.empty()) dirs.push_back(here);

#ifdef _WIN32
    // The installer's default, and the same path for the 32 bit build that some
    // machines still have.
    dirs.push_back("C:/Program Files/astap/");
    dirs.push_back("C:/Program Files (x86)/astap/");
#elif defined(__APPLE__)
    dirs.push_back("/Applications/ASTAP.app/Contents/MacOS/");
#else
    // Where the .deb and the star database packages put them, and where the
    // distribution-packaged builds expect them instead.
    dirs.push_back("/opt/astap/");
    dirs.push_back("/usr/share/astap/data/");
    dirs.push_back("/usr/share/astap/");
    dirs.push_back("/usr/local/share/astap/");
#endif
    return dirs;
  }

  namespace {
    bool file_exists(const std::string &name) {
      std::error_code ec;
      return std::filesystem::exists(name, ec);
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

  namespace {
    // The area files of a star database, mapped once and shared by everything
    // that reads them.
    //
    // A solve returns to an area over and over: four tiles cover most fields and
    // the spiral walks back and forth across them. Opening one of those files
    // costs about 35 microseconds on Windows, out of the 64 a whole 500 star
    // request takes, so the open is over half of it. Keeping the mapping alive
    // between visits removes the open, and the read with it.
    //
    // Process wide, because each search thread has its own StarDatabase and they
    // would otherwise map the same file once per thread. Bounded, so that a long
    // lived process working across several databases does not accumulate
    // mappings without limit; the bound is above the area count of either
    // database, so in the ordinary case nothing is ever evicted.
    class AreaMaps {
    public:
      std::shared_ptr<const MappedFile> get(const std::string &path) {
        std::lock_guard<std::mutex> lock(mutex_);
        const auto it = index_.find(path);
        if (it != index_.end()) {
          order_.splice(order_.begin(), order_, it->second.where);
          return it->second.map;
        }

        auto map = std::make_shared<MappedFile>();
        if (!map->open(path)) return nullptr;
        // Deliberately not advise_random(): records run bright to faint from the
        // start of the file and a request stops once it has enough stars, so the
        // access is a sequential prefix and read ahead is working with us. That
        // hint belongs to the index cache, which is queried scattered.
        order_.push_front(path);
        index_.emplace(path, Entry{map, order_.begin()});

        // Evicting drops this cache's reference and nothing else: a reader part
        // way through a mapping holds its own, so its pages stay under it.
        while (index_.size() > kLimit) {
          index_.erase(order_.back());
          order_.pop_back();
        }
        return map;
      }

    private:
      // Above the 1476 of the largest database, so a whole sky search maps every
      // area once and keeps them all.
      static constexpr size_t kLimit = 2048;

      struct Entry {
        std::shared_ptr<const MappedFile> map;
        std::list<std::string>::iterator where;
      };

      std::mutex mutex_;
      std::list<std::string> order_; // most recently used first
      std::unordered_map<std::string, Entry> index_;
    };

    AreaMaps &area_maps() {
      static AreaMaps maps;
      return maps;
    }
  } // namespace

  void StarDatabase::close() {
    area_map_.reset();
    records_ = nullptr;
    records_size_ = 0;
    position_ = 0;
  }

  bool StarDatabase::open_area(double telescope_dec, int area) {
    cos_telescope_dec_ = std::cos(telescope_dec); // here to save CPU time

    if (area != old_area_ || !area_map_) {
      close();

      const std::string namefile = name_database_ + "_" + area_filename(database_type_, area);
      area_map_ = area_maps().get(database_path_ + namefile);
      if (!area_map_) return false;
      if (area_map_->size() < 110) {
        // 10x11 is 110 bytes
        close();
        return false;
      }

      const uint8_t *header = area_map_->data();
      record_size_ = header[109] == ' ' ? 11 : header[109];
      database_version_ = header[108];

      records_ = header + 110;
      records_size_ = area_map_->size() - 110;
      old_area_ = area;
    }
    // else: this reader already holds the mapping for that area

    position_ = 0;
    return true;
  }

  bool StarDatabase::read_star(double telescope_ra, double telescope_dec, double field_diameter,
                               double &ra, double &dec, double &mag, double &b_v) {
    double delta_ra;
    do {
      bool header_record;
      do {
        // The whole record has to be there, not just its first byte. A file
        // whose length is not a whole number of records would otherwise be read
        // past its end, which used to run off a heap buffer quietly and would
        // now run off a mapping and fault.
        if (position_ + static_cast<size_t>(record_size_) > records_size_)
          return false; // end of file, no more data

        const uint8_t *rec = records_ + position_;
        position_ += static_cast<size_t>(record_size_);

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
            ra = ra_raw * kRaScale;
            const int32_t dec_raw = (static_cast<int32_t>(dec9_storage_) << 16) +
                                    (static_cast<int32_t>(rec[4]) << 8) + static_cast<int32_t>(rec[3]);
            dec = dec_raw * kDecScale;
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

  int StarDatabase::warm() {
    if (database_type_ == kDatabaseWideField) return 0;

    // What one request ever consumes: records run bright to faint and the read
    // stops as soon as it has enough stars, so the tail of a file is never
    // reached. Touching a page of every 4096 is what pulls it in; reading the
    // whole file would cost several times as much for nothing.
    constexpr size_t kWanted = 5 * 6 * 4 * 1024; // as much as read_star can use
    constexpr size_t kPage = 4096;

    int opened = 0;
    volatile uint8_t sink = 0;
    for (int area = 1; area <= database_type_; area++) {
      if (!open_area(0, area)) continue;
      opened++;
      const size_t n = std::min(kWanted, records_size_);
      for (size_t i = 0; i < n; i += kPage) sink = static_cast<uint8_t>(sink + records_[i]);
    }
    (void) sink;
    close();
    old_area_ = 9999999;
    return opened;
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
