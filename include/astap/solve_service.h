// One image in, one .ini out, with the index ladder held between calls.
//
// The index solver's cost is lopsided: reading a 2.7 GB ladder takes 1.6 s and
// the solve that follows takes 5 ms. A command line run pays the read every
// time, which is the right trade for a batch and the wrong one for an imaging
// application that solves a frame every few minutes. This class is the
// separation of the two: `load` once, `solve` as often as you like.
//
// The whole per image pipeline lives here - binning, star detection, the tier
// sweep, the second pass, the rescale back to original pixels and the .ini and
// .wcs files - so that `astap_index_solve` and the resident server produce the
// same output for the same image by construction rather than by agreement.
//
// Not thread safe. `solve` reads the star database through a file handle and
// tile cache that belong to this instance, so concurrent calls have to be
// serialised by the caller (the server takes a mutex); the ladder itself is
// immutable once loaded and could be shared, but the second pass is what makes
// the difference between a 0.2 px and a 0.14 px solution and it needs the
// database.

#pragma once

#include <string>
#include <vector>

#include "astap/index_solver.h"
#include "astap/quad_index.h"
#include "astap/star_database.h"
#include "astap/star_detection.h"
#include "astap/types.h"

namespace astap {
  // The default depth ladder. One tier reaches about a factor of two in image
  // star density either side of itself, so the 2x steps overlap and cover
  // 0.2 to 2000 stars/deg^2 continuously. Measured over the corpus, an image
  // whose density falls between two rungs is solved by one of them.
  const std::vector<double> &default_tier_ladder();

  // Rungs past the default ceiling, for a caller that asks for a deeper one.
  // They are not in the default ladder because they are most of the cache, and
  // because density matching reaches the same depths for free down to about
  // 0.25 degrees. See the small fields section of README.md.
  const std::vector<double> &deep_tier_ladder();

  // Parses a comma separated density list such as "60,125,250,500".
  std::vector<double> parse_density_list(const std::string &v);

  // The ladder a run should use. `tiers` names one outright when it is not
  // empty; otherwise the default one is taken and the deep rungs up to
  // `max_tier` are added, which is the common case - a narrow field needs
  // deeper rungs and nothing else about the ladder changes.
  std::vector<double> resolve_ladder(const std::string &tiers, double max_tier);

  // Replaces the extension of `path`, or appends when it has none.
  std::string change_file_ext(const std::string &path, const std::string &ext);

  // The database path is concatenated with a bare file name, so it has to end in
  // a separator. Windows takes '/' too, and a path the user already ended with
  // '\' is left as it is.
  std::string with_separator(std::string dir);

  // A header for an image that came from memory rather than from a file: the
  // dimensions the pipeline needs, and no FITS cards, so a .wcs written from it
  // carries the solution and nothing else.
  //
  // The pixel values are expected on the scale the detector produced, the same
  // 0..65535 range the histogram works in. An image normalised to 0..1 has no
  // stars as far as the detection is concerned.
  Header header_for_image(const ImageArray &img);

  struct SolveServiceSettings {
    // Directory holding the .1476 / .290 files. Empty searches
    // default_database_directories(), which is beside the executable and then
    // wherever ASTAP installs one; a directory named here is used on its own,
    // because a caller who says where the database is should not silently get a
    // different one.
    std::string database_path;
    std::string database = "auto"; // abbreviation such as d80, or "auto"

    // Sets the index hash bin width as well as the match tolerance, so a cache
    // belongs to exactly one value of it.
    double quad_tolerance = 0.007;
    std::vector<double> ladder; // empty selects default_tier_ladder()

    // One file holding the whole ladder, instead of the cached-per-rung layout.
    // Empty is the normal case and the one that lets ladders compose.
    std::string index_cache;
    bool use_cache = true; // false builds in memory and writes nothing
    bool rebuild = false; // build even when a usable cache exists
  };

  // What solving one image needs, whatever the pixels arrived in. Everything
  // here is a property of the image or of how hard to look; nothing in it names
  // a file. `SolveRequest` below is this plus a file to read and files to write.
  struct SolveParams {
    // Field diameter in degrees, 0 when unknown. It orders the tier sweep and
    // gives the minimum star size a plate scale to mean something against; it
    // never restricts which tiers are tried, so a wrong value costs time rather
    // than a solution.
    double fov = 0;
    int max_stars = 500;
    double min_star_size = 1.5; // arcsec, applied only when fov is known
    int downsample = 0; // 0 selects the factor automatically

    bool want_sip = false;
    bool refine = true; // the second pass against the database

    // What to call this image in the messages. The file name when it came from
    // one; a frame number or nothing at all when it did not.
    std::string label;
  };

  struct SolveRequest {
    std::string filename;
    std::string output_base; // empty names the outputs after `filename`

    SolveParams params;

    bool write_wcs = false;

    // Recorded verbatim in the .ini, as astap_cli records its own command line.
    std::string cmdline;
  };

  struct SolveOutcome {
    bool solved = false;
    int errorlevel = kErrNone;

    // The solution in original, unbinned pixel coordinates.
    Header head;
    SipCoefficients sip;

    double solve_seconds = 0; // the tier sweep alone
    double total_seconds = 0; // load, detect, sweep, second pass, files

    int stars = 0;
    int bin = 1;
    // Stars the solve actually used, and how many were detected. They differ
    // when density matching had to thin the list to reach a tier.
    size_t stars_used = 0;
    size_t stars_detected = 0;
    int nr_inliers = 0;
    int tiers_tried = 0;
    double tier_density = 0;
    bool many_quads_pass = false;
    bool refined = false;

    // What the command line front end would have printed for this image, in
    // order. The server relays these to its client so a solve through the pipe
    // reads like a solve without it.
    std::vector<std::string> messages;
  };

  class SolveService {
  public:
    // Selects the star database and gets a ladder, from the cache when one was
    // written by a run that would have built the same thing, otherwise by
    // building it and caching the result. Returns false when there is no usable
    // database, which is the one failure a later solve cannot work around.
    bool load(const SolveServiceSettings &s, const LogFn &log = nullptr);

    bool ready() const { return !tiers_.empty(); }

    // Solves an image already in memory, and touches no files at all.
    //
    // This is the whole per image pipeline - binning, star detection, the tier
    // sweep, the second pass, the rescale back to original pixels. `solve`
    // below is exactly this with a file read in front and the .ini and .wcs
    // written after, which is what keeps the two agreeing.
    //
    // `head` supplies the image dimensions and receives the solution. A caller
    // holding pixels and nothing else can get one from `header_for_image`; a
    // caller that read a file passes the header that came with it, so the FITS
    // cards travel into the .wcs.
    SolveOutcome solve_image(const ImageArray &img, Header &head, const SolveParams &p,
                             const LogFn &progress = nullptr);

    // Solves one image file. Writes the .ini always and the .wcs on request,
    // exactly as the command line front end does, and reports what happened.
    // `progress` receives the per step detail that -progress prints; the summary
    // lines land in `SolveOutcome::messages` either way.
    SolveOutcome solve(const SolveRequest &r, const LogFn &progress = nullptr);

    // --- what is resident, for a status report --------------------------------
    const std::string &database_name() const { return db_.name(); }
    const std::string &database_path() const { return db_.path(); }
    // The single cache file when one was named, otherwise the directory the
    // per-rung files live in.
    const std::string &cache_path() const { return cache_path_; }
    size_t tiers_from_cache() const { return from_cache_; }
    size_t tiers_built() const { return built_; }
    size_t tier_count() const { return tiers_.size(); }

    size_t quad_count() const;

    size_t bytes() const;

    std::vector<double> densities() const;

  private:
    // Reads whatever rungs of the wanted ladder `path` holds into the slots that
    // are still empty.
    void take_tiers_from(const std::string &path, bool quiet, const LogFn &log);

    // Selects the star database, searching the default locations when none was
    // named. Records where it was found in settings_.database_path.
    bool select_database(const LogFn &log);

    SolveServiceSettings settings_;
    StarDatabase db_;
    bool have_db_ = false;
    std::vector<QuadIndex> tiers_;
    std::vector<bool> have_;
    std::string cache_path_;
    size_t from_cache_ = 0;
    size_t built_ = 0;
  };
} // namespace astap
