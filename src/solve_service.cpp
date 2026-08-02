#include "astap/solve_service.h"

#include <chrono>
#include <cmath>
#include <cstdlib>

#include "astap/astro_math.h"
#include "astap/image/fits.h"
#include "astap/image/image_io.h"
#include "astap/solver.h"

namespace astap {
  namespace {
    using Clock = std::chrono::steady_clock;

    double secs(Clock::time_point a, Clock::time_point b) {
      return std::chrono::duration<double>(b - a).count();
    }

    std::string directory_of(const std::string &path) {
      const size_t slash = path.find_last_of("/\\");
      return slash == std::string::npos ? std::string(".") : path.substr(0, slash);
    }

    // Mono, binned copy of the image, plus the binning factor that produced it.
    //
    // Binning past about 1200 px on the long side throws away the stars the
    // solve depends on, so the automatic factor stops there. The index solver
    // works on this copy, and its solution is scaled back to original pixels
    // afterwards.
    int bin_mono(const ImageArray &img, int width, int height, int forced, ImageArray &out) {
      int bin = forced;
      if (bin <= 0) {
        bin = 1;
        while (std::max(width, height) / (bin + 1) >= 1200) bin++;
      }
      out = ImageArray(1, height / bin, width / bin);
      for (int y = 0; y < out.height(); y++)
        for (int x = 0; x < out.width(); x++) {
          double v = 0;
          for (int c = 0; c < img.colours(); c++)
            for (int i = 0; i < bin; i++)
              for (int j = 0; j < bin; j++) v += img.at(c, y * bin + i, x * bin + j);
          out.at(0, y, x) = static_cast<float>(v / (img.colours() * bin * bin));
        }
      return bin;
    }

    // The index solver returns a WCS in binned pixels. Binned pixel xb covers
    // original pixels [xb*bin, xb*bin+bin), so its centre sits at
    // xb*bin + (bin-1)/2 in original zero-based coordinates. Rewriting the
    // solution in those coordinates is a scale of the pixel size and a shift of
    // the reference pixel; the rotation is unchanged.
    void scale_solution_to_original(const IndexSolveResult &r, int bin, Header &head) {
      head.ra0 = r.ra0;
      head.dec0 = r.dec0;
      head.crpix1 = (r.crpix1 - 1) * bin + (bin - 1) / 2.0 + 1;
      head.crpix2 = (r.crpix2 - 1) * bin + (bin - 1) / 2.0 + 1;
      head.cdelt1 = r.cdelt1 / bin;
      head.cdelt2 = r.cdelt2 / bin;
      head.crota1 = r.crota1;
      head.crota2 = r.crota2;
      head.cd1_1 = r.cd1_1 / bin;
      head.cd1_2 = r.cd1_2 / bin;
      head.cd2_1 = r.cd2_1 / bin;
      head.cd2_2 = r.cd2_2 / bin;
    }

    // SIP coefficients from the binned frame, rewritten for original pixels.
    //
    // The polynomial takes an offset from the reference pixel and returns a
    // correction, both in pixels, and offsets scale exactly: u_orig = bin *
    // u_bin, because the reference pixel is mapped by the same relation as every
    // other. A term c*u^p*v^q therefore becomes c * bin^(1-p-q) in the original
    // frame.
    void scale_sip_to_original(SipCoefficients &sip, int bin) {
      if (!sip.valid || bin == 1) return;
      double (*tables[4])[4] = {sip.a, sip.b, sip.ap, sip.bp};
      for (double (*t)[4] : tables)
        for (int p = 0; p < 4; p++)
          for (int q = 0; q < 4; q++)
            if (t[p][q] != 0) t[p][q] *= std::pow(static_cast<double>(bin), 1 - p - q);
    }
  } // namespace

  const std::vector<double> &default_tier_ladder() {
    static const std::vector<double> ladder = {0.5, 1, 2, 4, 8, 16, 32, 60, 125, 250, 500, 900};
    return ladder;
  }

  const std::vector<double> &deep_tier_ladder() {
    static const std::vector<double> ladder = {1800, 3600};
    return ladder;
  }

  std::vector<double> resolve_ladder(const std::string &tiers, double max_tier) {
    if (!tiers.empty()) return parse_density_list(tiers);
    std::vector<double> ladder = default_tier_ladder();
    for (double d : deep_tier_ladder())
      if (d <= max_tier) ladder.push_back(d);
    return ladder;
  }

  std::vector<double> parse_density_list(const std::string &v) {
    std::vector<double> out;
    size_t p = 0;
    while (p <= v.size()) {
      const size_t c = v.find(',', p);
      const std::string tok = v.substr(p, c == std::string::npos ? std::string::npos : c - p);
      if (!tok.empty()) out.push_back(std::atof(tok.c_str()));
      if (c == std::string::npos) break;
      p = c + 1;
    }
    return out;
  }

  std::string change_file_ext(const std::string &path, const std::string &ext) {
    size_t slash = path.find_last_of("/\\");
    size_t dot = path.find_last_of('.');
    if (dot == std::string::npos || (slash != std::string::npos && dot < slash)) return path + ext;
    return path.substr(0, dot) + ext;
  }

  std::string with_separator(std::string dir) {
    if (!dir.empty() && dir.back() != '/' && dir.back() != '\\') dir += '/';
    return dir;
  }

  // Reads whatever rungs of the wanted ladder this file holds into the slots
  // still empty. A file may hold one rung (the cached-per-rung layout) or the
  // whole ladder (a named cache, or one written before the rungs were split
  // apart); both are just a list of tiers to match against what is wanted.
  void SolveService::take_tiers_from(const std::string &path, bool quiet, const LogFn &log) {
    auto say = [&](const std::string &m) {
      if (log && !quiet) log(m);
    };
    std::string err;
    QuadIndexFile info;
    if (!read_index_file_header(path, info, &err)) return;
    // A cache is only usable when it was built the way this run would build it.
    // Anything else is rebuilt rather than silently half-used.
    if (std::fabs(info.quad_tolerance - settings_.quad_tolerance) > 1e-12) {
      say("Index cache " + path + " was built at tolerance " + float_to_str(info.quad_tolerance, 4) +
          ", rebuilding at " + float_to_str(settings_.quad_tolerance, 4) + ".");
      return;
    }
    std::vector<QuadIndex> in;
    if (!load_index_file(path, in, 0, 0, &err)) {
      say("Index cache " + path + " unusable (" + err + "), rebuilding it.");
      return;
    }
    for (QuadIndex &ix : in)
      for (size_t i = 0; i < settings_.ladder.size(); i++)
        if (!have_[i] && std::fabs(ix.settings().star_density - settings_.ladder[i]) < 1e-9) {
          tiers_[i] = std::move(ix);
          have_[i] = true;
          from_cache_++;
          break;
        }
  }

  bool SolveService::load(const SolveServiceSettings &s, const LogFn &log) {
    settings_ = s;
    if (settings_.ladder.empty()) settings_.ladder = default_tier_ladder();
    from_cache_ = 0;
    built_ = 0;

    auto say = [&](const std::string &m) {
      if (log) log(m);
    };

    have_db_ = db_.select(with_separator(settings_.database_path), settings_.database, 1.0);
    const std::string db_name = have_db_ ? db_.name() : "unknown";
    const int db_type = have_db_ ? db_.database_type() : 0;

    // The pre-split layout: one file for the whole ladder. Still read when it is
    // there, and still written when a caller names a cache outright.
    const std::string ladder_file =
        settings_.index_cache.empty()
            ? default_index_cache_path(db_name, db_type, settings_.quad_tolerance)
            : settings_.index_cache;
    cache_path_ = settings_.index_cache.empty() ? directory_of(ladder_file) : ladder_file;

    // A rung at a time. Each is cached in its own file, so the rungs this run
    // shares with a previous one are read back and only the new ones are built:
    // raising the ceiling costs the deep rung, not the whole ladder.
    tiers_.assign(settings_.ladder.size(), QuadIndex());
    have_.assign(settings_.ladder.size(), false);

    const auto t0 = Clock::now();
    if (settings_.use_cache && !settings_.rebuild) {
      if (!settings_.index_cache.empty()) {
        take_tiers_from(ladder_file, false, log);
      } else {
        for (double d : settings_.ladder)
          take_tiers_from(index_tier_cache_path(db_name, db_type, settings_.quad_tolerance, d),
                          true, log);
        // Whatever is still missing may be in a ladder file written before the
        // rungs were split apart.
        if (from_cache_ < settings_.ladder.size()) take_tiers_from(ladder_file, true, log);
      }
    }

    std::vector<double> missing;
    for (size_t i = 0; i < settings_.ladder.size(); i++)
      if (!have_[i]) missing.push_back(settings_.ladder[i]);

    if (missing.empty()) {
      say("Index: " + std::to_string(tiers_.size()) + " tiers read from cache in " +
          float_to_str(secs(t0, Clock::now()), 2) + " sec.");
    } else {
      // Building is the only step that needs the database. Rungs written by an
      // earlier run stand on their own, which is why the check comes here rather
      // than at the top: a machine that has the cache but has since moved its
      // database directory can still solve, it just cannot rebuild.
      if (!have_db_) {
        say("No star database found in " + with_separator(settings_.database_path));
        tiers_.clear();
        return false;
      }
      QuadIndexSettings qs;
      qs.quad_tolerance = settings_.quad_tolerance;
      say("Building " + std::to_string(missing.size()) + " index tier(s) from " + db_.name() +
          "; each one happens once.");
      // One pass over the database covers every missing rung: a tile is read
      // once at the deepest of them and the shallower ones take a prefix.
      std::vector<QuadIndex> made;
      if (!build_tiers(db_, qs, missing, made)) {
        say("Could not read the star database in " + with_separator(settings_.database_path));
        tiers_.clear();
        return false;
      }
      for (QuadIndex &ix : made) {
        const double d = ix.settings().star_density;
        if (settings_.use_cache && settings_.index_cache.empty()) {
          const std::string path =
              index_tier_cache_path(db_.name(), db_.database_type(), settings_.quad_tolerance, d);
          std::vector<QuadIndex> one;
          one.push_back(ix);
          if (!ensure_parent_directory(path) || !save_index_file(path, one))
            say("Warning: could not write the index cache to " + path);
        }
        for (size_t i = 0; i < settings_.ladder.size(); i++)
          if (!have_[i] && std::fabs(d - settings_.ladder[i]) < 1e-9) {
            tiers_[i] = std::move(ix);
            have_[i] = true;
            built_++;
            break;
          }
      }
      if (settings_.use_cache && !settings_.index_cache.empty()) {
        if (ensure_parent_directory(ladder_file) && save_index_file(ladder_file, tiers_))
          say("Cached to " + ladder_file);
        else
          say("Warning: could not write the index cache to " + ladder_file);
      }
      say("Index: " + std::to_string(tiers_.size()) + " tiers (" + std::to_string(from_cache_) +
          " cached, " + std::to_string(built_) + " built), " + std::to_string(quad_count()) +
          " quads, ready in " + float_to_str(secs(t0, Clock::now()), 2) + " sec.");
    }

    // A rung that could not be built at all would otherwise sit in the sweep as
    // an empty index; drop it so the tier count reported is the tier count used.
    for (size_t i = tiers_.size(); i-- > 0;)
      if (!have_[i]) {
        tiers_.erase(tiers_.begin() + static_cast<long>(i));
        have_.erase(have_.begin() + static_cast<long>(i));
      }
    if (tiers_.empty()) {
      say("No usable index tiers.");
      return false;
    }
    return true;
  }

  size_t SolveService::quad_count() const {
    size_t n = 0;
    for (const QuadIndex &ix : tiers_) n += ix.size();
    return n;
  }

  size_t SolveService::bytes() const {
    size_t n = 0;
    for (const QuadIndex &ix : tiers_) n += ix.bytes();
    return n;
  }

  std::vector<double> SolveService::densities() const {
    std::vector<double> d;
    d.reserve(tiers_.size());
    for (const QuadIndex &ix : tiers_) d.push_back(ix.settings().star_density);
    return d;
  }

  SolveOutcome SolveService::solve(const SolveRequest &r, const LogFn &progress) {
    SolveOutcome o;
    const auto t_start = Clock::now();
    auto say = [&](const std::string &m) {
      o.messages.push_back(m);
      if (progress) progress(m);
    };
    const std::string out_base = r.output_base.empty() ? r.filename : r.output_base;
    auto finish = [&](bool solved, int level) {
      o.solved = solved;
      o.errorlevel = level;
      o.total_seconds = secs(t_start, Clock::now());
      return o;
    };

    Header head;
    ImageArray img;
    const ImageLoadResult lr = load_image(r.filename, head, img);
    if (!lr.ok) {
      say(lr.error);
      write_ini(change_file_ext(out_base, ".ini"), false, head, r.cmdline, kErrImageRead, "");
      o.head = head;
      return finish(false, kErrImageRead);
    }
    if (!lr.warning.empty()) say(lr.warning);

    ImageArray small;
    const int bin = bin_mono(img, head.width, head.height, r.downsample, small);
    o.bin = bin;

    // A star size in arcseconds only means something once the plate scale is
    // known, which a blind solve does not have. With a field size it does.
    double hfd_min = 0.8;
    if (r.fov > 0) {
      const double arcsec_per_px = r.fov * 3600 / head.height;
      hfd_min = std::max(0.8, r.min_star_size / (bin * arcsec_per_px));
    }

    Header bhead = head;
    Histogram hist;
    get_background(0, small, bhead, true, true, r.max_stars, hist);
    RowList stars;
    double mean_hfd = 0;
    find_stars(small, bhead, hfd_min, r.max_stars, stars, mean_hfd, progress ? progress : LogFn());
    o.stars = static_cast<int>(stars.count());
    if (stars.count() < 4) {
      say("Not enough stars detected in " + r.filename);
      write_ini(change_file_ext(out_base, ".ini"), false, head, r.cmdline, kErrNotEnoughStars, "");
      o.head = head;
      return finish(false, kErrNotEnoughStars);
    }

    if (progress)
      say("Image: " + std::to_string(stars.count()) + " stars, binning " + std::to_string(bin));

    // The tier sweep is ordered by the image's star density when the field size
    // is known; without it every tier is tried, cheapest first.
    double density_hint = 0;
    if (r.fov > 0) {
      const double w_deg = r.fov * head.width / head.height;
      density_hint = stars.count() / std::max(1e-9, r.fov * w_deg);
    }

    const auto s0 = Clock::now();
    IndexSolveResult res =
        solve_stars_with_tiers(tiers_, stars, small.width(), small.height(), {}, density_hint);
    o.solve_seconds = secs(s0, Clock::now());

    if (!res.solved) {
      say("No solution for " + r.filename + " (" + res.reason + ", " +
          std::to_string(res.nr_matches) + " candidate quads)");
      write_ini(change_file_ext(out_base, ".ini"), false, head, r.cmdline, kErrNone, "");
      o.head = head;
      return finish(false, kErrNoSolution);
    }

    // Second pass: with the position known, read the database once at it and
    // redo the match there. This is what lifts the accuracy and produces enough
    // matched quads for a distortion fit.
    if (!r.refine && r.want_sip)
      say("Note: SIP needs the second pass, which is disabled. No SIP written.");
    if (r.refine && !have_db_ && r.want_sip)
      say("Note: SIP needs the star database, which was not found. No SIP written.");
    if (r.refine && have_db_) {
      const auto p0 = Clock::now();
      IndexRefineResult ref =
          refine_with_database(db_, stars, small.width(), small.height(), res, {}, r.want_sip);
      const double refine_secs = secs(p0, Clock::now());
      if (ref.ok && progress)
        say("Second pass: " + std::to_string(ref.nr_quads) + " quads matched against " +
            std::to_string(ref.nr_candidates) + " database quads in " +
            float_to_str(refine_secs * 1000, 1) + " ms" +
            (ref.residual_before >= 0
                 ? ", residual " + float_to_str(ref.residual_before, 3) + "\" -> " +
                       float_to_str(ref.residual_after, 3) + "\""
                 : "") +
            (ref.kept ? "." : ", discarded: " + ref.reason));
      if (!ref.ok && progress) say("Second pass skipped: " + ref.reason);
      if (ref.sip_valid) {
        o.sip = ref.sip;
        scale_sip_to_original(o.sip, bin);
      } else if (r.want_sip) {
        say("No SIP coefficients: " + ref.reason);
      }
    }

    scale_solution_to_original(res, bin, head);
    o.head = head;
    o.nr_inliers = res.nr_inliers;
    o.tiers_tried = res.tiers_tried;
    o.tier_density = res.tier_density;
    o.many_quads_pass = res.many_quads_pass;
    o.refined = res.refined;
    o.stars_used = res.stars_used;
    o.stars_detected = res.stars_detected;

    say("Solution found: " + prepare_ra(head.ra0, ": ") + " " + prepare_dec(head.dec0, "d "));
    say("Solved in " + float_to_str(o.solve_seconds, 3) + " sec. Scale " +
        float_to_str(std::fabs(head.cdelt2) * 3600, 4) + "\"/px, rotation " +
        float_to_str(head.crota2, 2) + "d, " + std::to_string(res.nr_inliers) +
        " quads at depth tier " + float_to_str(res.tier_density, 1) + " stars/deg^2 (" +
        std::to_string(res.tiers_tried) + " tiers tried" +
        (res.many_quads_pass ? ", larger quad groups" : "") +
        (res.stars_used < res.stars_detected
             ? ", brightest " + std::to_string(res.stars_used) + " of " +
                   std::to_string(res.stars_detected) + " stars"
             : "") +
        (res.refined ? ", refined" : "") + (o.sip.valid ? ", SIP" : "") + ").");

    write_ini(change_file_ext(out_base, ".ini"), true, head, r.cmdline, kErrNone, "");

    if (r.write_wcs) {
      const std::string comment =
          "Solved in " + float_to_str(o.solve_seconds, 3) + " sec by the index method.";
      update_solution_cards(head.cards, head, o.sip, true, comment);
      remove_key(head.cards, "NAXIS1  =");
      remove_key(head.cards, "NAXIS2  =");
      update_integer(head.cards, "NAXIS   =",
                     " / Minimal header                                 ", 0);
      update_integer(head.cards, "BITPIX  =",
                     " /                                                ", 8);
      write_fits_header_file(change_file_ext(out_base, ".wcs"), head.cards);
    }

    return finish(true, kErrNone);
  }
} // namespace astap
