// The ASTAP astrometric plate solving method by Han Kleijn.
//
//      => Image <=                                    |  => Star database <=
// 1) Find background, noise and star level            |
// 2) Find stars and their CCD x, y position           | Extract the same amount of stars
//    (standard coordinates)                           | (area corrected) from the area of
//                                                     | interest and convert the alpha, delta
//                                                     | equatorial coordinates into standard
//                                                     | coordinates (rigid method)
// 3) Construct the smallest irregular tetrahedrons     | Idem
//    (quads) of four stars, calculate the six          |
//    distances and the mean x,y position of the quad   |
// 4) Sort the six distances, d1 longest, d6 shortest   | Idem
// 5) Scale them as (d1, d2/d1 .. d6/d1): the hash code | Idem
//
//                        => matching process <=
// 6) Find quad hash code matches where the five ratios match within a small tolerance.
// 7) Calculate the median of the d1_found/d1_reference ratios and remove the quads
//    outside a small tolerance.
// 8) Solve the overdetermined system A * S = X_ref for the six plate constants (LSQ_FIT),
//    then derive the image centre position, the pixel scale and the rotation.
//
// Ported from solve_image() in unit_command_line_solving.pas.

#pragma once

#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "astap/matching.h"
#include "astap/star_database.h"
#include "astap/star_detection.h"
#include "astap/types.h"

namespace astap {
  struct SolverSettings {
    std::string database_path; // directory holding the star database files
    std::string star_database = "auto"; // "auto" or an abbreviation such as d80, g05, w08
    double search_fov = 0; // image height in degrees, 0 selects auto FOV
    bool fov_specified = false; // true when search_fov overrides the header
    double radius_search = 180; // search radius in degrees
    int max_stars = 500; // maximum number of stars to use
    double quad_tolerance = 0.007; // hash code matching tolerance
    double min_star_size = 1.5; // arcsec, ignores hot pixels
    int downsample = 0; // 0 = auto, else the binning factor
    bool force_oversize = false; // "slow" mode, more overlap while searching
    bool check_pattern_filter = false; // normalise the bayer pattern of raw OSC images
    bool add_sip = false; // add SIP distortion coefficients
    bool show_log = false; // log every search step
  };

  class Solver {
  public:
    explicit Solver(SolverSettings settings) : settings_(std::move(settings)) {
    }

    // Called for every progress message. Defaults to no output.
    void set_log(LogFn log) { log_ = std::move(log); }

    // Solves `img`. `head` supplies the start position (ra0, dec0) and the image
    // dimensions and receives the solution. Returns true when a solution is found.
    bool solve(ImageArray img, Header &head);

    int errorlevel() const { return errorlevel_; }
    const std::string &warning() const { return warning_str_; }
    const SipCoefficients &sip() const { return sip_; }
    // Number of matching quads used for the solution, and the number found.
    int nr_references() const { return match_.nr_references; }
    int nr_references2() const { return match_.nr_references2; }
    const SolutionVector &solution_x() const { return match_.solution_vector_x; }
    const SolutionVector &solution_y() const { return match_.solution_vector_y; }
    double mount_ra() const { return ra_mount_; }
    double mount_dec() const { return dec_mount_; }

    void set_mount(double ra, double dec) {
      ra_mount_ = ra;
      dec_mount_ = dec;
    }

    double solved_seconds() const { return solved_seconds_; }
    double search_offset() const { return sep_search_; } // radians
    double magnitude_limit() const { return mag2_ / 10; }

    // Where the time went. The three `_cpu` entries are summed over all worker
    // threads, so they exceed the wall clock; the two `_wall` entries do not
    // overlap and do add up to the run time.
    struct Timing {
      double image_wall = 0;     // binning, background and star detection, per FOV attempt
      double spiral_wall = 0;    // the whole squared spiral search
      double read_stars_cpu = 0; // star database I/O and record decoding
      double quads_cpu = 0;      // database quad construction
      double match_cpu = 0;      // hash matching and the least squares fit
    };

    const Timing &timing() const { return timing_; }

  private:
    // Make mono, bin and crop.
    static void bin_mono_and_crop(int &binning, double crop, const ImageArray &img, ImageArray &img2,
                                  const LogFn &log);

    static void convert_mono(ImageArray &img);

    // Normalise the bayer pattern, avoids colour shifts for raw OSC images.
    static void apply_check_pattern_filter(ImageArray &img, const LogFn &log);

    void bin_and_find_stars(const ImageArray &img, Header &head, int binfactor, double cropping,
                            double hfd_min, int max_stars, RowList &starlist,
                            std::string &short_warning);

    int report_binning_astrometric(double height, double arcsec_per_px) const;

    // Read the stars of the area of interest from the database and convert them
    // to standard coordinates relative to telescope_ra, telescope_dec.
    // `db` and `mag_out` are passed explicitly so that several search threads can
    // each read from their own file handle and tile cache.
    bool read_stars(StarDatabase &db, double telescope_ra, double telescope_dec, double search_field,
                    int nrstars_required, RowList &starlist, double &mag_out);

    bool add_sip_coefficients(const Header &head, double ra_database, double dec_database);

    // Per-thread state for the parallel spiral search. Each worker needs its own
    // database file handle and tile cache, its own star list and its own match
    // state, so that the positions of a batch do not interfere.
    struct SearchWorker {
      StarDatabase database;
      MatchState match;
      double t_read = 0, t_quads = 0, t_match = 0;  // CPU time, this worker only
      RowList starlist;
      double mag2 = 0;
    };

    // Creates (or reuses) one worker per thread and points each at the database
    // selected for this solve.
    void prepare_workers();

    void say(const std::string &s) const {
      if (log_) log_(s);
    }

    SolverSettings settings_;
    LogFn log_;
    StarDatabase database_;
    // unique_ptr because a worker holds an ifstream, which is not movable in a
    // way that survives vector growth on every standard library.
    std::vector<std::unique_ptr<SearchWorker> > workers_;
    MatchState match_;
    Histogram histogram_;
    SipCoefficients sip_;

    std::string warning_str_;
    int errorlevel_ = kErrNone;
    double mag2_ = 0; // magnitude * 10 of the faintest star used
    double ra_mount_ = 99999;
    double dec_mount_ = 99999;
    double solved_seconds_ = 0;
    double sep_search_ = 0;
    Timing timing_;
  };

  // Writes the solution in the astap_cli .ini format.
  bool write_ini(const std::string &filename, bool solution, const Header &head,
                 const std::string &cmdline, int errorlevel, const std::string &warning);

  // Updates or inserts a keyword in the header cards, used to build the .wcs file.
  void update_text(std::vector<std::string> &cards, const std::string &key,
                   const std::string &value_and_comment);

  void update_float(std::vector<std::string> &cards, const std::string &key,
                    const std::string &comment, double x);

  void update_integer(std::vector<std::string> &cards, const std::string &key,
                      const std::string &comment, long x);

  void remove_key(std::vector<std::string> &cards, const std::string &key);

  // Fills the header cards with the solution (and the SIP coefficients when
  // present) so they can be written to a .wcs file.
  void update_solution_cards(std::vector<std::string> &cards, const Header &head,
                             const SipCoefficients &sip, bool solution, const std::string &comment);
} // namespace astap
