// The file half of SpiralService, matching solve_service_files.cpp: read an
// image, solve it, write the .ini and the .wcs. Split out for the same reason -
// this is the only part that needs an image decoder.

#include <chrono>

#include "astap/astro_math.h"
#include "astap/image/fits.h"
#include "astap/image/image_io.h"
#include "astap/spiral_service.h"

namespace astap {
  SolveOutcome SpiralService::solve(const SpiralRequest &r, const LogFn &progress) {
    using Clock = std::chrono::steady_clock;
    const auto t_start = Clock::now();
    const std::string out_base = r.output_base.empty() ? r.filename : r.output_base;

    std::vector<std::string>
    pre;
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

    SpiralParams p = r.params;
    if (p.label.empty()) p.label = r.filename;
    // The search starts from the position load_image put in the header, which
    // the caller can override through p.ra/p.dec. What the mount reported is a
    // separate thing and is only used to report the offset, so it is passed
    // through as itself rather than folded into the start position.
    if (p.mount_ra >= 999) p.mount_ra = lr.ra_mount;
    if (p.mount_dec >= 999) p.mount_dec = lr.dec_mount;

    SolveOutcome o = solve_image(std::move(img), head, p, progress);
    o.messages.insert(o.messages.begin(), pre.begin(), pre.end());

    write_ini(change_file_ext(out_base, ".ini"), o.solved, o.head, r.cmdline, o.errorlevel, "");

    if (o.solved && r.write_wcs) {
      const std::string comment = "Solved in " + float_to_str(o.solve_seconds, 1) + " sec.";
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
