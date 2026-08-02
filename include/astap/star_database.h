// Reader for the HNSKY .290 / .1476 star databases and the .001 wide field
// database. Ported from unit_command_line_star_database.pas and
// unit_command_line_stars_wide_field.pas.
//
// The sky is divided into 290 respectively 1476 areas of nearly equal surface,
// one file per area. Each file starts with a 110 byte header of which the last
// byte holds the record size. Stars are sorted from bright to faint in 0.1
// magnitude steps; every group is preceded by a special record with RA=$FFFFFF
// that carries the magnitude and the high declination byte.
//
// A star record of 5 bytes:
//   ra7 ra8 ra9   RA as a 3 byte word,   resolution 360*3600/(2^24-1) = 0.077"
//   dec7 dec8     the low two bytes of the 3 byte two's complement declination,
//                 the high byte comes from the preceding header record,
//                 resolution 90*3600/(2^23-1) = 0.039"
// A record of 6 bytes adds one shortint with the Gaia colour (B-V)*50.

#pragma once

#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

namespace astap {
  enum DatabaseType : int {
    kDatabaseWideField = 1, // .001 files
    kDatabase290 = 290, // .290 files
    kDatabase1476 = 1476, // .1476 files
  };

  // The directory holding the running executable, with a trailing separator.
  // Empty when the platform will not say.
  std::string executable_directory();

  // Where to look for a star database when nobody said where it is.
  //
  // A database is gigabytes and is downloaded once, so the overwhelmingly likely
  // situation is that the machine already has one because it already has ASTAP.
  // Looking where ASTAP puts it costs nothing and saves the user from repeating
  // themselves; a caller who names a directory is still taken at their word and
  // nothing below is consulted.
  //
  // In order: beside this executable, which is where a self-contained
  // distribution puts both; then the places ASTAP's own installers use.
  std::vector<std::string> default_database_directories();

  class StarDatabase {
  public:
    // Select a star database, report false when none is found.
    // `database` is either "auto" or an abbreviation such as d80, v50, g05, w08.
    bool select(const std::string &database_path, const std::string &database, double fov,
                std::string *warning = nullptr);

    int database_type() const { return database_type_; }
    const std::string &name() const { return name_database_; }
    const std::string &path() const { return database_path_; }

    // Adopt an already made selection, so several reader instances can share one
    // choice while keeping their own file handle and cache.
    void configure(const std::string &database_path, const std::string &name, int type) {
      database_path_ = database_path;
      name_database_ = name;
      database_type_ = type;
    }

    // Find up to four database areas covering the square image. `fov` is in
    // radians and must stay below the tile size (9.53 degrees for .290 files,
    // 5.14 degrees for .1476 files), otherwise a tile beyond the next one could
    // be selected.
    void find_areas(double ra1, double dec1, double fov, int &area1, int &area2, int &area3,
                    int &area4, double &frac1, double &frac2, double &frac3, double &frac4) const;

    // Open the database file for one area (1 based). Returns false when the file
    // could not be opened.
    bool open_area(double telescope_dec, int area);

    // Read the next star within `field_diameter` (radians) of the telescope
    // position. Returns false when the file holds no more data.
    // `mag` is the magnitude * 10, `b_v` the colour information (unused by the
    // solver, -999 when not available).
    bool read_star(double telescope_ra, double telescope_dec, double field_diameter, double &ra,
                   double &dec, double &mag, double &b_v);

    // Load the wide field (.001) database: an int32 star count followed by
    // count*3 floats holding magnitude, ra and dec.
    bool read_stars_wide_field(const std::string &database_path);

    const std::vector<float> &wide_field_stars() const { return wide_field_stars_; }

    void close();

  private:
    static std::string area_filename(int database_type, int area);

    void area_and_boundaries(double ra1, double dec1, int &area_nr, double &spaceE, double &spaceW,
                             double &spaceN, double &spaceS) const;

    std::string database_path_;
    std::string name_database_;
    int database_type_ = kDatabase1476;
    int database_version_ = 0;

    std::ifstream file_;
    bool file_open_ = false;
    int old_area_ = 9999999;

    int record_size_ = 11;
    int8_t dec9_storage_ = 0;
    double cos_telescope_dec_ = 1.0;

    // The whole file is cached: re-using it is about 35% faster for a field of
    // view of 0.5 degrees, where reading the database is 60% of the total time.
    std::vector<uint8_t> cache_;
    size_t cache_size_ = 0;
    size_t cache_valid_pos_ = 0;
    size_t cache_position_ = 0;

    std::vector<float> wide_field_stars_;
    std::string wide_database_;
  };
} // namespace astap
