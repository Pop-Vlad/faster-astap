// TIFF reader, the equivalent of read_tiff() in unit_tiff_unthreaded.pas: grey
// scale or RGB, 8 or 16 bit integer or 32 bit float samples, rows stored bottom
// up in the FITS convention, an 8 bit sample scaled by 257 and a float sample by
// 65535. libtiff takes care of the compression schemes the original handles
// itself (none, LZW, Deflate) and of the rest besides.
//
// Anything outside that shape (tiled, planar, palette, unusual sample widths)
// goes through libtiff's RGBA reader, which normalises it to 8 bit RGB.
//
// Only compiled with a decoder when the build found libtiff.

#include <cstdint>
#include <string>
#include <vector>

#include "image_loaders.h"

#ifdef ASTAP_WITH_TIFF
#include <cstdarg>

#include <tiffio.h>
#endif

namespace astap {
  namespace imageio {
#ifdef ASTAP_WITH_TIFF
    namespace {
      void silence(const char *, const char *, va_list) {
      }
    } // namespace

    ImageLoadResult load_tiff(const std::string &filename, Header &head, ImageArray &img) {
      ImageLoadResult res;

      // libtiff writes its complaints to stderr by default, which would end up
      // in the middle of the solver output.
      TIFFSetErrorHandler(silence);
      TIFFSetWarningHandler(silence);

      TIFF *t = TIFFOpen(filename.c_str(), "r");
      if (!t) {
        res.error = "Error, accessing the file!";
        return res;
      }

      head = Header();
      img = ImageArray();

      uint32_t width = 0, height = 0;
      uint16_t samples = 1, bits = 8, format = SAMPLEFORMAT_UINT, planar = PLANARCONFIG_CONTIG;
      TIFFGetField(t, TIFFTAG_IMAGEWIDTH, &width);
      TIFFGetField(t, TIFFTAG_IMAGELENGTH, &height);
      TIFFGetFieldDefaulted(t, TIFFTAG_SAMPLESPERPIXEL, &samples);
      TIFFGetFieldDefaulted(t, TIFFTAG_BITSPERSAMPLE, &bits);
      TIFFGetFieldDefaulted(t, TIFFTAG_SAMPLEFORMAT, &format);
      TIFFGetFieldDefaulted(t, TIFFTAG_PLANARCONFIG, &planar);

      if (width == 0 || height == 0) {
        TIFFClose(t);
        res.error = "Error, invalid TIFF dimensions!";
        return res;
      }

      const bool is_float = bits == 32 && format == SAMPLEFORMAT_IEEEFP;
      const bool direct = !TIFFIsTiled(t) && planar == PLANARCONFIG_CONTIG &&
                          (samples == 1 || samples == 3) && (bits == 8 || bits == 16 || is_float);

      head.width = static_cast<int>(width);
      head.height = static_cast<int>(height);

      if (direct) {
        const int naxis3 = samples == 3 ? 3 : 1;
        img.resize(naxis3, head.height, head.width);

        std::vector<uint8_t> row(static_cast<size_t>(TIFFScanlineSize(t)));
        for (uint32_t y = 0; y < height; y++) {
          if (TIFFReadScanline(t, row.data(), y) < 0) {
            TIFFClose(t);
            res.error = "Error, reading the TIFF image data!";
            return res;
          }
          const int dst_y = head.height - 1 - static_cast<int>(y); // FITS is bottom up
          for (int k = 0; k < naxis3; k++) {
            float *dst = img.row(k, dst_y);
            for (uint32_t x = 0; x < width; x++) {
              const size_t i = static_cast<size_t>(x) * samples + k;
              if (bits == 8) {
                dst[x] = row[i] * 257;
              } else if (bits == 16) {
                dst[x] = reinterpret_cast<const uint16_t *>(row.data())[i];
              } else {
                dst[x] = reinterpret_cast<const float *>(row.data())[i] * 65535;
              }
            }
          }
        }
        TIFFClose(t);
        synthesise_header(head, is_float ? -32 : 16, naxis3, 65535);
        res.ok = true;
        return res;
      }

      // Everything else through the RGBA reader, which hands back 8 bit pixels
      // with the first row at the bottom, the order the image array wants.
      std::vector<uint32_t> raster(static_cast<size_t>(width) * height);
      if (!TIFFReadRGBAImageOriented(t, width, height, raster.data(), ORIENTATION_BOTLEFT, 0)) {
        TIFFClose(t);
        res.error = "Error, reading the TIFF image data!";
        return res;
      }
      TIFFClose(t);

      img.resize(3, head.height, head.width);
      for (uint32_t y = 0; y < height; y++) {
        const uint32_t *src = raster.data() + static_cast<size_t>(y) * width;
        float *r = img.row(0, static_cast<int>(y));
        float *g = img.row(1, static_cast<int>(y));
        float *b = img.row(2, static_cast<int>(y));
        for (uint32_t x = 0; x < width; x++) {
          r[x] = TIFFGetR(src[x]) * 257;
          g[x] = TIFFGetG(src[x]) * 257;
          b[x] = TIFFGetB(src[x]) * 257;
        }
      }

      synthesise_header(head, 16, 3, 65535);
      res.warning = "TIFF read through the RGBA path, reduced to 8 bits per channel.";
      res.ok = true;
      return res;
    }
#else
    ImageLoadResult load_tiff(const std::string &, Header &, ImageArray &) {
      return missing_support("TIFF", "libtiff");
    }
#endif
  } // namespace imageio
} // namespace astap
