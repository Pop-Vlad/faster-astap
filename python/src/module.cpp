// The Python extension: both solvers, with the load-once/solve-many lifecycle
// they were built around.
//
// This layer is deliberately thin. It converts an array into an ImageArray,
// releases the GIL, takes the instance's lock and calls the same solve_image()
// the command line front ends call. Everything that could be written in Python
// is, in faster_astap/; what is here is what cannot be.
//
// Neither service is thread safe — the workers, their database handles and the
// tile caches belong to one instance — so every entry point takes a per instance
// mutex. That makes concurrent calls safe rather than parallel; a caller wanting
// parallelism should hold one solver per thread.

#include <memory>
#include <mutex>
#include <string>
#include <vector>

#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "astap/solve_service.h"
#include "astap/spiral_service.h"
#include "astap/star_database.h"

#ifdef ASTAP_PYTHON_IMAGE_IO
#include "astap/image/image_io.h"
#endif

namespace py = pybind11;
using namespace astap;

namespace {
  // Anything array-like, copied into the solver's own layout.
  //
  // forcecast means an integer camera frame arrives without the caller having to
  // convert it; c_style means the copy below can walk it in order. Accepted
  // shapes are (h, w), (c, h, w) and (h, w, c) — the last because that is what
  // an image library hands back and it would be unkind to make the caller
  // transpose it.
  using FloatArray = py::array_t<float, py::array::c_style | py::array::forcecast>;

  ImageArray to_image(const FloatArray &a) {
    const py::buffer_info info = a.request();
    int colours = 1, height = 0, width = 0;
    bool channels_last = false;

    if (info.ndim == 2) {
      height = static_cast<int>(info.shape[0]);
      width = static_cast<int>(info.shape[1]);
    } else if (info.ndim == 3) {
      // (c, h, w) or (h, w, c). Three or fewer planes at the front is a colour
      // axis; at the back it is one too. An image three pixels tall is not
      // something a plate solver will ever see, so the ambiguity is theoretical.
      if (info.shape[0] <= 4 && info.shape[2] > 4) {
        colours = static_cast<int>(info.shape[0]);
        height = static_cast<int>(info.shape[1]);
        width = static_cast<int>(info.shape[2]);
      } else if (info.shape[2] <= 4) {
        height = static_cast<int>(info.shape[0]);
        width = static_cast<int>(info.shape[1]);
        colours = static_cast<int>(info.shape[2]);
        channels_last = true;
      } else {
        throw py::value_error("a 3 dimensional image must be (colours, height, width) or "
                              "(height, width, colours), with at most 4 colours");
      }
    } else {
      throw py::value_error("expected a 2 or 3 dimensional image, got " +
                            std::to_string(info.ndim) + " dimensions");
    }

    if (width < 8 || height < 8)
      throw py::value_error("the image is too small to detect stars in: " +
                            std::to_string(width) + "x" + std::to_string(height));

    const float *src = static_cast<const float *>(info.ptr);

    // Star detection bins the pixel values into a 0..65535 histogram, so an
    // image scaled to 0..1 has no stars in it as far as this is concerned. That
    // fails as "not enough stars detected", which tells the caller nothing about
    // what they actually did wrong, so it is worth catching here.
    double peak = 0;
    const size_t n = static_cast<size_t>(colours) * height * width;
    for (size_t i = 0; i < n; i++)
      if (src[i] > peak) peak = src[i];
    if (peak > 0 && peak <= 1.0)
      throw py::value_error(
          "the image peaks at " + std::to_string(peak) +
          ", which looks like it was normalised to 0..1. Star detection works on the scale the "
          "detector produced (0..65535); multiply by 65535 if that is what happened.");

    ImageArray img(colours, height, width);
    for (int c = 0; c < colours; c++)
      for (int y = 0; y < height; y++) {
        float *row = img.row(c, y);
        if (channels_last) {
          const float *in = src + (static_cast<size_t>(y) * width) * colours + c;
          for (int x = 0; x < width; x++) row[x] = in[static_cast<size_t>(x) * colours];
        } else {
          const float *in = src + (static_cast<size_t>(c) * height + y) * width;
          for (int x = 0; x < width; x++) row[x] = in[x];
        }
      }
    return img;
  }

  // A log callback that reacquires the GIL, since the solve that calls it has
  // dropped it. Null when the caller passed nothing, which costs the solver
  // nothing at all.
  LogFn make_log(const py::object &fn) {
    if (fn.is_none()) return LogFn();
    // Shared so the lambda can be copied into the solver without copying the
    // handle itself, which would need the GIL each time.
    auto held = std::make_shared<py::object>(fn);
    return [held](const std::string &line) {
      py::gil_scoped_acquire gil;
      try {
        (*held)(line);
      } catch (py::error_already_set &e) {
        // A raising progress callback must not take the solve down with it.
        e.discard_as_unraisable("faster_astap progress callback");
      }
    };
  }

  // Both services get the same wrapper: a lock, and the GIL dropped around
  // anything that takes real time. Named for what it is rather than what it
  // does, because astap::Solver is already a thing and this is not it.
  template <typename Service>
  struct Resident {
    Service service;
    std::mutex mutex;
    bool loaded = false;
  };
} // namespace

PYBIND11_MODULE(_core, m) {
  m.doc() = "Bindings for the faster-astap plate solvers.";

  m.attr("__version__") = "0.1.0";

  m.def("default_database_directories", &default_database_directories,
        "Where a star database is looked for when none is named, in order.");

#ifdef ASTAP_PYTHON_IMAGE_IO
  m.attr("has_image_io") = true;
  m.def("supported_image_extensions", &supported_image_extensions,
        "Image file extensions this build can read.");
#else
  m.attr("has_image_io") = false;
  m.def("supported_image_extensions", []() { return std::string(); },
        "Empty: this build reads no image files, only arrays.");
#endif

  // --- results --------------------------------------------------------------
  py::class_<SipCoefficients>(m, "Sip", "Third order SIP distortion coefficients.")
      .def_property_readonly("valid", [](const SipCoefficients &s) { return s.valid; })
      .def_property_readonly("a",
                             [](const SipCoefficients &s) {
                               return py::array_t<double>({4, 4}, &s.a[0][0]);
                             })
      .def_property_readonly("b",
                             [](const SipCoefficients &s) {
                               return py::array_t<double>({4, 4}, &s.b[0][0]);
                             })
      .def_property_readonly("ap",
                             [](const SipCoefficients &s) {
                               return py::array_t<double>({4, 4}, &s.ap[0][0]);
                             })
      .def_property_readonly("bp", [](const SipCoefficients &s) {
        return py::array_t<double>({4, 4}, &s.bp[0][0]);
      });

  py::class_<Header>(m, "Header", "The WCS the solver produced, in original pixels.")
      .def_readonly("width", &Header::width)
      .def_readonly("height", &Header::height)
      .def_readonly("ra0", &Header::ra0, "Image centre right ascension, radians.")
      .def_readonly("dec0", &Header::dec0, "Image centre declination, radians.")
      .def_readonly("crpix1", &Header::crpix1)
      .def_readonly("crpix2", &Header::crpix2)
      .def_readonly("cdelt1", &Header::cdelt1, "Pixel width, degrees.")
      .def_readonly("cdelt2", &Header::cdelt2, "Pixel height, degrees.")
      .def_readonly("crota1", &Header::crota1)
      .def_readonly("crota2", &Header::crota2, "Rotation at the centre, degrees.")
      .def_readonly("cd1_1", &Header::cd1_1)
      .def_readonly("cd1_2", &Header::cd1_2)
      .def_readonly("cd2_1", &Header::cd2_1)
      .def_readonly("cd2_2", &Header::cd2_2)
      .def_property_readonly("cards", [](const Header &h) { return h.cards; },
                            "The FITS header cards, when the image came from a file.");

  py::class_<SolveOutcome>(m, "Outcome", "What a solve produced.")
      .def_readonly("solved", &SolveOutcome::solved)
      .def_readonly("errorlevel", &SolveOutcome::errorlevel)
      .def_readonly("head", &SolveOutcome::head)
      .def_readonly("sip", &SolveOutcome::sip)
      .def_readonly("solve_seconds", &SolveOutcome::solve_seconds)
      .def_readonly("total_seconds", &SolveOutcome::total_seconds)
      .def_readonly("stars", &SolveOutcome::stars)
      .def_readonly("bin", &SolveOutcome::bin)
      .def_readonly("stars_used", &SolveOutcome::stars_used)
      .def_readonly("stars_detected", &SolveOutcome::stars_detected)
      .def_readonly("nr_inliers", &SolveOutcome::nr_inliers)
      .def_readonly("tiers_tried", &SolveOutcome::tiers_tried)
      .def_readonly("tier_density", &SolveOutcome::tier_density)
      .def_readonly("many_quads_pass", &SolveOutcome::many_quads_pass)
      .def_readonly("refined", &SolveOutcome::refined)
      .def_readonly("messages", &SolveOutcome::messages);

  // --- the index solver -----------------------------------------------------
  using Index = Resident<SolveService>;
  py::class_<Index>(m, "IndexSolver")
      .def(py::init<>())
      .def(
          "load",
          [](Index &self, const std::string &database_path, const std::string &database,
             double quad_tolerance, const std::vector<double> &ladder,
             const std::string &index_cache, bool use_cache, bool rebuild, const py::object &log) {
            SolveServiceSettings s;
            s.database_path = database_path;
            s.database = database;
            s.quad_tolerance = quad_tolerance;
            s.ladder = ladder;
            s.index_cache = index_cache;
            s.use_cache = use_cache;
            s.rebuild = rebuild;
            const LogFn fn = make_log(log);
            bool ok;
            {
              py::gil_scoped_release unlock;
              std::lock_guard<std::mutex> lock(self.mutex);
              ok = self.service.load(s, fn);
            }
            self.loaded = ok;
            return ok;
          },
          py::arg("database_path") = "", py::arg("database") = "auto",
          py::arg("quad_tolerance") = 0.007, py::arg("ladder") = std::vector<double>(),
          py::arg("index_cache") = "", py::arg("use_cache") = true, py::arg("rebuild") = false,
          py::arg("log") = py::none())
      .def(
          "solve_array",
          [](Index &self, const FloatArray &a, double fov, int max_stars, double min_star_size,
             int downsample, bool want_sip, bool refine, const std::string &label,
             const py::object &progress) {
            SolveParams p;
            p.fov = fov;
            p.max_stars = max_stars;
            p.min_star_size = min_star_size;
            p.downsample = downsample;
            p.want_sip = want_sip;
            p.refine = refine;
            p.label = label;
            // Converted while the GIL is still held: it reads a Python object.
            ImageArray img = to_image(a);
            const LogFn fn = make_log(progress);
            Header head = header_for_image(img);
            py::gil_scoped_release unlock;
            std::lock_guard<std::mutex> lock(self.mutex);
            return self.service.solve_image(img, head, p, fn);
          },
          py::arg("image"), py::arg("fov") = 0.0, py::arg("max_stars") = 500,
          py::arg("min_star_size") = 1.5, py::arg("downsample") = 0, py::arg("want_sip") = false,
          py::arg("refine") = true, py::arg("label") = "", py::arg("progress") = py::none())
      .def_property_readonly("ready", [](Index &self) { return self.service.ready(); })
      .def_property_readonly("database",
                             [](Index &self) { return self.service.database_name(); })
      .def_property_readonly("database_path",
                             [](Index &self) { return self.service.database_path(); })
      .def_property_readonly("cache_path", [](Index &self) { return self.service.cache_path(); })
      .def_property_readonly("tier_count", [](Index &self) { return self.service.tier_count(); })
      .def_property_readonly("quad_count", [](Index &self) { return self.service.quad_count(); })
      .def_property_readonly("bytes", [](Index &self) { return self.service.bytes(); })
      .def_property_readonly("densities", [](Index &self) { return self.service.densities(); })
      .def_property_readonly("tiers_from_cache",
                             [](Index &self) { return self.service.tiers_from_cache(); })
      .def_property_readonly("tiers_built", [](Index &self) { return self.service.tiers_built(); });

  // --- the spiral port ------------------------------------------------------
  using Spiral = Resident<SpiralService>;
  py::class_<Spiral>(m, "SpiralSolver")
      .def(py::init<>())
      .def(
          "load",
          [](Spiral &self, const std::string &database_path, const std::string &database, bool warm,
             unsigned threads, const py::object &log) {
            SpiralServiceSettings s;
            s.database_path = database_path;
            s.database = database;
            s.warm = warm;
            s.threads = threads;
            const LogFn fn = make_log(log);
            bool ok;
            {
              py::gil_scoped_release unlock;
              std::lock_guard<std::mutex> lock(self.mutex);
              ok = self.service.load(s, fn);
            }
            self.loaded = ok;
            return ok;
          },
          py::arg("database_path") = "", py::arg("database") = "auto", py::arg("warm") = false,
          py::arg("threads") = 0u, py::arg("log") = py::none())
      .def(
          "solve_array",
          [](Spiral &self, const FloatArray &a, double ra, double dec, double fov,
             bool fov_specified, double radius, int max_stars, double quad_tolerance,
             double min_star_size, int downsample, bool force_oversize, bool check_pattern_filter,
             bool want_sip, const std::string &label, const py::object &progress) {
            SpiralParams p;
            p.ra = ra;
            p.dec = dec;
            p.fov = fov;
            p.fov_specified = fov_specified;
            p.radius = radius;
            p.max_stars = max_stars;
            p.quad_tolerance = quad_tolerance;
            p.min_star_size = min_star_size;
            p.downsample = downsample;
            p.force_oversize = force_oversize;
            p.check_pattern_filter = check_pattern_filter;
            p.want_sip = want_sip;
            p.label = label;
            ImageArray img = to_image(a);
            const LogFn fn = make_log(progress);
            Header head = header_for_image(img);
            py::gil_scoped_release unlock;
            std::lock_guard<std::mutex> lock(self.mutex);
            return self.service.solve_image(std::move(img), head, p, fn);
          },
          py::arg("image"), py::arg("ra") = 99999.0, py::arg("dec") = 99999.0, py::arg("fov") = 0.0,
          py::arg("fov_specified") = false, py::arg("radius") = 180.0, py::arg("max_stars") = 500,
          py::arg("quad_tolerance") = 0.007, py::arg("min_star_size") = 1.5,
          py::arg("downsample") = 0, py::arg("force_oversize") = false,
          py::arg("check_pattern_filter") = false, py::arg("want_sip") = false,
          py::arg("label") = "", py::arg("progress") = py::none())
      .def_property_readonly("ready", [](Spiral &self) { return self.service.ready(); })
      .def_property_readonly("database",
                             [](Spiral &self) { return self.service.database_name(); })
      .def_property_readonly("database_path",
                             [](Spiral &self) { return self.service.database_path(); })
      .def_property_readonly("areas_warmed",
                             [](Spiral &self) { return self.service.areas_warmed(); });

#ifdef ASTAP_PYTHON_IMAGE_IO
  m.def(
      "load_image_file",
      [](const std::string &path) {
        Header head;
        ImageArray img;
        ImageLoadResult lr;
        {
          py::gil_scoped_release unlock;
          lr = load_image(path, head, img);
        }
        if (!lr.ok) throw py::value_error(lr.error);
        // (colours, height, width), which is the solver's own order.
        py::array_t<float> out({img.colours(), img.height(), img.width()});
        auto view = out.mutable_unchecked<3>();
        for (int c = 0; c < img.colours(); c++)
          for (int y = 0; y < img.height(); y++)
            for (int x = 0; x < img.width(); x++) view(c, y, x) = img.at(c, y, x);
        py::dict meta;
        meta["ra_mount"] = lr.ra_mount;
        meta["dec_mount"] = lr.dec_mount;
        meta["focallen"] = lr.focallen;
        meta["warning"] = lr.warning;
        meta["cards"] = head.cards;
        meta["cdelt2"] = head.cdelt2;
        meta["ra0"] = head.ra0;
        meta["dec0"] = head.dec0;
        return py::make_tuple(out, meta);
      },
      py::arg("path"),
      "Reads an image file into a (colours, height, width) array of the pixel values, plus what "
      "the header said. Raises ValueError when the file cannot be read.");
#endif
}
