#include "astap/spiral_service.h"

#include <chrono>
#include <cmath>

#include "astap/astro_math.h"
#include "astap/parallel.h"

namespace astap {
  namespace {
    using Clock = std::chrono::steady_clock;

    double secs(Clock::time_point a, Clock::time_point b) {
      return std::chrono::duration<double>(b - a).count();
    }
  } // namespace

  bool SpiralService::load(const SpiralServiceSettings &s, const LogFn &log) {
    settings_ = s;
    ready_ = false;
    warmed_ = 0;
    auto say = [&](const std::string &m) {
      if (log) log(m);
    };

    if (settings_.threads > 0) set_thread_count(settings_.threads);

    // The same search SolveService does: a named directory is taken at its word,
    // otherwise the places a database is likely to already be are tried in turn.
    StarDatabase probe;
    bool found = false;
    if (!settings_.database_path.empty()) {
      settings_.database_path = with_separator(settings_.database_path);
      found = probe.select(settings_.database_path, settings_.database, 1.0);
    } else {
      for (const std::string &dir : default_database_directories()) {
        if (!probe.select(with_separator(dir), settings_.database, 1.0)) continue;
        settings_.database_path = with_separator(dir);
        found = true;
        break;
      }
    }
    if (!found) {
      say("Error, no star database found. Download and install a star database.");
      return false;
    }
    database_name_ = probe.name();
    say("Star database: " + database_name_ + " in " + settings_.database_path);

    if (settings_.warm) {
      const auto t0 = Clock::now();
      warmed_ = probe.warm();
      if (warmed_ > 0)
        say("Warmed " + std::to_string(warmed_) + " areas of " + database_name_ + " in " +
            float_to_str(secs(t0, Clock::now()), 2) + " sec.");
    }

    // Built once. Its workers, and the mapped area files they hold, then live
    // across every solve rather than being made again per frame.
    SolverSettings base;
    base.database_path = settings_.database_path;
    base.star_database = settings_.database;
    solver_.reset(new Solver(base));
    ready_ = true;
    return true;
  }

  SolveOutcome SpiralService::solve_image(ImageArray img, Header &head, const SpiralParams &p,
                                          const LogFn &progress) {
    SolveOutcome o;
    const auto t_start = Clock::now();
    auto say = [&](const std::string &m) {
      o.messages.push_back(m);
      if (progress) progress(m);
    };
    const std::string what = p.label.empty() ? std::string("the image") : p.label;

    if (!ready_ || !solver_) {
      say("No star database loaded.");
      o.errorlevel = kErrNoStarDatabase;
      o.total_seconds = secs(t_start, Clock::now());
      return o;
    }

    if (head.width <= 0 || head.height <= 0) {
      head.width = img.width();
      head.height = img.height();
    }
    // A caller who named a position means it; otherwise whatever the header
    // carried is left in place for the search to start from.
    if (p.ra < 999) head.ra0 = p.ra;
    if (p.dec < 999) head.dec0 = p.dec;

    SolverSettings s;
    s.database_path = settings_.database_path;
    s.star_database = settings_.database;
    s.search_fov = p.fov;
    s.fov_specified = p.fov_specified;
    s.radius_search = p.radius;
    s.max_stars = p.max_stars;
    s.quad_tolerance = p.quad_tolerance;
    s.min_star_size = p.min_star_size;
    s.downsample = p.downsample;
    s.force_oversize = p.force_oversize;
    s.check_pattern_filter = p.check_pattern_filter;
    s.add_sip = p.want_sip;
    s.show_log = static_cast<bool>(progress);
    solver_->configure(s);

    // The solver's own running commentary is progress, not summary: it goes to
    // the caller's callback but does not fill up the outcome's messages.
    solver_->set_log(progress ? progress : LogFn());
    solver_->set_mount(p.mount_ra, p.mount_dec);

    const bool solved = solver_->solve(std::move(img), head);

    o.solved = solved;
    o.errorlevel = solver_->errorlevel();
    o.head = head;
    o.sip = solver_->sip();
    o.solve_seconds = solver_->solved_seconds();
    o.nr_inliers = solver_->nr_references();
    o.bin = 1; // the port bins internally and reports in original pixels

    if (solved) {
      say("Solution found: " + prepare_ra(head.ra0, ": ") + " " + prepare_dec(head.dec0, "d "));
      say("Solved in " + float_to_str(o.solve_seconds, 3) + " sec. Scale " +
          float_to_str(std::fabs(head.cdelt2) * 3600, 4) + "\"/px, rotation " +
          float_to_str(head.crota2, 2) + "d, " + std::to_string(o.nr_inliers) +
          " quads. Offset was " +
          distance_to_string(solver_->search_offset(), solver_->search_offset()) + "." +
          (o.sip.valid ? " SIP." : ""));
    } else if (o.errorlevel == kErrNotEnoughStars) {
      say("Not enough stars detected in " + what);
    } else {
      say("No solution for " + what);
      if (o.errorlevel == kErrNone) o.errorlevel = kErrNoSolution;
    }

    o.total_seconds = secs(t_start, Clock::now());
    return o;
  }
} // namespace astap
