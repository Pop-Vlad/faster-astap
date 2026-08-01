#include "astap/image/image_io.h"

#include <algorithm>
#include <cctype>
#include <initializer_list>

#include "astap/astro_math.h"
#include "astap/image/fits.h"
#include "image_loaders.h"

namespace astap {
  namespace {
    std::string lower_extension(const std::string &path) {
      const size_t slash = path.find_last_of("/\\");
      const size_t dot = path.find_last_of('.');
      if (dot == std::string::npos || (slash != std::string::npos && dot < slash))
        return std::string();
      std::string ext = path.substr(dot);
      for (char &c: ext) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
      return ext;
    }

    bool one_of(const std::string &ext, std::initializer_list<const char *> list) {
      for (const char *e: list)
        if (ext == e) return true;
      return false;
    }

    // The raw camera formats LibRaw handles. Not part of astap_cli, which has no
    // raw support at all; the GUI shells out to dcraw or LibRaw's unprocessed_raw
    // and reads back the resulting PGM, which is what this reproduces in process.
    bool is_raw_extension(const std::string &ext) {
      return one_of(ext, {".cr2", ".cr3", ".crw", ".nef", ".nrw", ".arw", ".srf", ".sr2",
                          ".orf", ".rw2", ".raf", ".dng", ".pef", ".raw", ".3fr", ".fff",
                          ".iiq", ".mos", ".mef", ".mrw", ".erf", ".kdc", ".dcr", ".srw",
                          ".x3f", ".rwl", ".dcs", ".cap", ".bay"});
    }
  } // namespace

  ImageLoadResult load_image(const std::string &filename, Header &head, ImageArray &img) {
    const std::string ext = lower_extension(filename);

    if (one_of(ext, {".fit", ".fits", ".fts", ".fz", ".new"})) return load_fits(filename, head, img);
    if (one_of(ext, {".ppm", ".pgm", ".pfm", ".pbm"})) return imageio::load_pnm(filename, head, img);
    if (ext == ".bmp") return imageio::load_bmp(filename, head, img);
    if (ext == ".png") return imageio::load_png(filename, head, img);
    if (one_of(ext, {".jpg", ".jpeg"})) return imageio::load_jpeg(filename, head, img);
    if (one_of(ext, {".tif", ".tiff"})) return imageio::load_tiff(filename, head, img);
    if (is_raw_extension(ext)) return imageio::load_raw(filename, head, img);

    ImageLoadResult res;
    res.error = "Error, unsupported file type '" + (ext.empty() ? filename : ext) +
                "'. This build reads: " + supported_image_extensions();
    return res;
  }

  std::string supported_image_extensions() {
    std::string s = ".fit .fits .fts .fz .new .ppm .pgm .pfm .bmp";
#ifdef ASTAP_WITH_PNG
    s += " .png";
#endif
#ifdef ASTAP_WITH_JPEG
    s += " .jpg .jpeg";
#endif
#ifdef ASTAP_WITH_TIFF
    s += " .tif .tiff";
#endif
#ifdef ASTAP_WITH_LIBRAW
    s += " and raw camera files (.cr2 .cr3 .nef .arw .dng ...)";
#endif
    return s;
  }

  namespace imageio {
    void synthesise_header(Header &head, int bitpix, int naxis3, double datamax) {
      head.bitpix = bitpix;
      head.naxis3 = naxis3;
      head.naxis = naxis3 == 1 ? 2 : 3;
      head.datamax_org = datamax;

      // Same keywords in the same order as head1[0..10] in the original, with
      // NAXIS3 skipped for a mono image. The solver only ever appends to this,
      // so the cards are written out in order rather than through the update
      // helpers of the solver module.
      head.cards.clear();
      auto card = [&head](const std::string &key, const std::string &value_and_comment) {
        std::string c = key;
        c.resize(9, ' ');
        c += value_and_comment;
        c.resize(80, ' ');
        head.cards.push_back(c);
      };
      auto integer_card = [&card](const std::string &key, long value, const std::string &comment) {
        std::string v = std::to_string(value);
        if (v.size() < 20) v = std::string(20 - v.size(), ' ') + v;
        card(key, v + comment);
      };

      card("SIMPLE  =", "                    T / FITS header                                    ");
      integer_card("BITPIX  =", bitpix, " / Bits per entry                                 ");
      integer_card("NAXIS   =", head.naxis, " / Number of dimensions                           ");
      integer_card("NAXIS1  =", head.width, " / length of x axis                               ");
      integer_card("NAXIS2  =", head.height, " / length of y axis                               ");
      if (naxis3 != 1)
        integer_card("NAXIS3  =", naxis3, " / length of z axis (mostly colors)               ");
      card("EQUINOX =", "               2000.0 / Equinox of coordinates                         ");
      integer_card("DATAMIN =", 0, " / Minimum data value                             ");
      integer_card("DATAMAX =", static_cast<long>(pround(datamax)),
                   " / Maximum data value                             ");
      card("BZERO   =", "                  0.0 / Physical_value = BZERO + BSCALE * array_value  ");
      card("BSCALE  =", "                  1.0 / Physical_value = BZERO + BSCALE * array_value  ");
      card("COMMENT 1", "  Written by ASTAP, Astrometric STAcking Program. www.hnsky.org        ");
      card("END", "");
    }

    ImageLoadResult missing_support(const std::string &format, const std::string &library) {
      ImageLoadResult res;
      res.error = "Error, this build reads no " + format + " files. Rebuild with " + library +
                  " available, or convert the image to FITS.";
      return res;
    }
  } // namespace imageio
} // namespace astap
