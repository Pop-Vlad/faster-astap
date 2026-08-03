// The file half of SolveService: read an image, solve it, write the .ini and
// the .wcs.
//
// It lives apart from the pipeline it wraps because it is the only part that
// needs an image decoder. Everything solve_image() does is in astap_solver,
// which depends on nothing outside the standard library, so a caller that
// already holds pixels — a Python binding, a camera driver, a test — can solve
// without libpng, libjpeg, libtiff and LibRaw coming along for the ride.

#include <chrono>

#include "astap/astro_math.h"
#include "astap/image/fits.h"
#include "astap/image/image_io.h"
#include "astap/solve_service.h"
#include "astap/solver.h"

namespace astap {
  SolveOutcome SolveService::solve(const SolveRequest &r, const LogFn &progress) {
    using Clock = std::chrono::steady_clock;
    const auto t_start = Clock::now();
    const std::string out_base = r.output_base.empty() ? r.filename : r.output_base;

    // Said before the pipeline runs, so they have to be put in front of the
    // messages it produces rather than appended to them.
    std::vector<std::string> pre;
    auto say = [&](const std::string &m) {
      pre.push_back(m);
      if (progress) progress(m);
    };

    Header head;
    ImageArray img;
    const ImageLoadResult lr = load_image(r.filename, head, img);
    if (!lr.ok) {
      say(lr.error);
      write_ini(change_file_ext(out_base, ".ini"), false, head, r.cmdline, kErrImageRead, "");
      SolveOutcome o;
      o.messages = pre;
      o.head = head;
      o.errorlevel = kErrImageRead;
      o.total_seconds = std::chrono::duration<double>(Clock::now() - t_start).count();
      return o;
    }
    if (!lr.warning.empty()) say(lr.warning);

    SolveParams p = r.params;
    if (p.label.empty()) p.label = r.filename;
    SolveOutcome o = solve_image(img, head, p, progress);
    o.messages.insert(o.messages.begin(), pre.begin(), pre.end());

    // What astap_cli records for each outcome, which is not the same as what the
    // process exits with: a failure to find a solution is a status of its own to
    // the caller, and a plain "not solved" in the file.
    const int ini_level = o.errorlevel == kErrNotEnoughStars ? kErrNotEnoughStars : kErrNone;
    write_ini(change_file_ext(out_base, ".ini"), o.solved, o.head, r.cmdline, ini_level, "");

    if (o.solved && r.write_wcs) {
      const std::string comment =
          "Solved in " + float_to_str(o.solve_seconds, 3) + " sec by the index method.";
      update_solution_cards(o.head.cards, o.head, o.sip, true, comment);
      // A minimal header: the .wcs file carries no image data.
      remove_key(o.head.cards, "NAXIS1  =");
      remove_key(o.head.cards, "NAXIS2  =");
      update_integer(o.head.cards, "NAXIS   =",
                     " / Minimal header                                 ", 0);
      update_integer(o.head.cards, "BITPIX  =",
                     " /                                                ", 8);
      write_fits_header_file(change_file_ext(out_base, ".wcs"), o.head.cards);
    }

    o.total_seconds = std::chrono::duration<double>(Clock::now() - t_start).count();
    return o;
  }
} // namespace astap
