// The spiral search port with the same lifecycle as SolveService: set up once,
// solve as often as you like.
//
// The port has less to set up than the index solver — there is no ladder to
// build — but the setup it does have is worth keeping between images. Selecting
// the database, starting the worker threads and giving each its own reader all
// happen once here instead of per frame, and with `warm` the area files are
// pulled into the page cache before the first solve rather than during it.
//
// Deliberately the same shape as SolveService, down to the SolveOutcome it
// returns, so that a caller can swap the reference solver for the fast one by
// changing the type it holds. The fields of SolveOutcome that describe the tier
// sweep have no meaning here and keep their defaults.
//
// Not thread safe, for the same reason SolveService is not: the workers and
// their tile caches belong to this instance.

#pragma once

#include <memory>
#include <string>

#include "astap/solve_service.h"
#include "astap/solver.h"
#include "astap/star_database.h"
#include "astap/types.h"

namespace astap {
  struct SpiralServiceSettings {
    // Directory holding the star database. Empty searches
    // default_database_directories(), exactly as SolveService does.
    std::string database_path;
    std::string database = "auto";  // abbreviation such as d80, or "auto"

    // Touch every area file during load(), so the first solve does not pay to
    // fault them in. See StarDatabase::warm.
    bool warm = false;

    // 0 leaves the thread pool as it is, which is one worker per hardware
    // thread. Set it to bound a solver sharing a machine with a capture loop.
    unsigned threads = 0;
  };

  // The per image half. The database and the thread pool are in the settings
  // above because they belong to the service, not to the frame.
  struct SpiralParams {
    // Start position in radians, 99999 leaving the one already in the header.
    // This is where the spiral begins, so a good value is worth a great deal of
    // time and a bad one costs it.
    double ra = 99999;
    double dec = 99999;

    // Where the mount said it was pointing, if it said. This does not steer the
    // search — `ra` and `dec` above do — it only gives the solved offset
    // something to be measured from in the report.
    double mount_ra = 99999;
    double mount_dec = 99999;

    // Image height in degrees, 0 selects the automatic sweep. `fov_specified`
    // is what stops the header's plate scale being believed over this.
    double fov = 0;
    bool fov_specified = false;
    double radius = 180;  // search radius in degrees

    int max_stars = 500;
    double quad_tolerance = 0.007;
    double min_star_size = 1.5;  // arcsec
    int downsample = 0;          // 0 selects the factor automatically
    bool force_oversize = false; // "slow": more overlap while searching
    bool check_pattern_filter = false;
    bool want_sip = false;

    // What to call this image in the messages.
    std::string label;
  };

  struct SpiralRequest {
    std::string filename;
    std::string output_base;  // empty names the outputs after `filename`
    SpiralParams params;
    bool write_wcs = false;
    std::string cmdline;
  };

  class SpiralService {
  public:
    // Selects the star database and gets the workers ready. Returns false when
    // there is no usable database, which is the one failure a later solve cannot
    // work around.
    bool load(const SpiralServiceSettings &s, const LogFn &log = nullptr);

    bool ready() const { return ready_; }

    // Solves an image already in memory, touching no files. `head` supplies the
    // dimensions and the start position, and receives the solution.
    SolveOutcome solve_image(ImageArray img, Header &head, const SpiralParams &p,
                             const LogFn &progress = nullptr);

    // Solves one image file, writing the .ini and, on request, the .wcs, the
    // same way astap_solve does.
    SolveOutcome solve(const SpiralRequest &r, const LogFn &progress = nullptr);

    const std::string &database_name() const { return database_name_; }
    const std::string &database_path() const { return settings_.database_path; }
    int areas_warmed() const { return warmed_; }

  private:
    SpiralServiceSettings settings_;
    std::string database_name_;
    // A pointer because the Solver is what holds the workers, and it is built
    // once load() knows which database they should read.
    std::unique_ptr<Solver> solver_;
    bool ready_ = false;
    int warmed_ = 0;
  };
} // namespace astap
