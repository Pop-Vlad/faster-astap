// Raw camera files (CR2, CR3, NEF, ARW, DNG and the rest LibRaw knows).
//
// astap_cli has no raw support at all. The ASTAP GUI gets it by running dcraw
// or LibRaw's unprocessed_raw over the file and reading back the 16 bit PGM
// they write: the undemosaiced sensor frame, margins included, no white
// balance, no gamma. That is exactly what this does in process, so the pixels
// the solver sees are the same ones the GUI works with, without the detour
// through a temporary file.
//
// The Bayer mosaic is handed to the solver as it comes off the sensor. Star
// detection works on it directly: a star covers many pixels, so the mosaic
// only costs a little sensitivity. -check applies the check pattern filter for
// raw OSC images, the same as in the original.
//
// Only compiled with a decoder when the build found LibRaw.

#include <cstdint>
#include <string>

#include "image_loaders.h"

#ifdef ASTAP_WITH_LIBRAW
#include <libraw/libraw.h>
#endif

namespace astap {
  namespace imageio {
#ifdef ASTAP_WITH_LIBRAW
    ImageLoadResult load_raw(const std::string &filename, Header &head, ImageArray &img) {
      ImageLoadResult res;

      LibRaw raw;
      int status = raw.open_file(filename.c_str());
      if (status != LIBRAW_SUCCESS) {
        res.error = std::string("Error, accessing the file! ") + libraw_strerror(status);
        return res;
      }
      status = raw.unpack();
      if (status != LIBRAW_SUCCESS) {
        res.error = std::string("Error, unpacking the raw file! ") + libraw_strerror(status);
        return res;
      }

      const libraw_rawdata_t &data = raw.imgdata.rawdata;
      const libraw_image_sizes_t &sizes = raw.imgdata.sizes;
      if (data.raw_image == nullptr) {
        res.error = "Error, this raw file holds no single channel sensor image "
            "(Foveon and already demosaiced files are not supported).";
        return res;
      }

      head = Header();
      img = ImageArray();

      const int width = sizes.raw_width;
      const int height = sizes.raw_height;
      if (width <= 0 || height <= 0) {
        res.error = "Error, invalid raw image dimensions!";
        return res;
      }
      // raw_pitch is in bytes, the image itself is 16 bit.
      const size_t stride = sizes.raw_pitch ? sizes.raw_pitch / 2 : static_cast<size_t>(width);

      head.width = width;
      head.height = height;
      img.resize(1, height, width);

      // Row order as in the PGM the GUI's converters write, which the PPM
      // loader also stores unflipped.
      float measured_max = 0;
      for (int y = 0; y < height; y++) {
        const uint16_t *src = data.raw_image + static_cast<size_t>(y) * stride;
        float *dst = img.row(0, y);
        for (int x = 0; x < width; x++) {
          dst[x] = src[x];
          if (dst[x] > measured_max) measured_max = dst[x];
        }
      }

      synthesise_header(head, 16, 1, measured_max > 255 ? 65535 : 255);
      res.ok = true;
      return res;
    }
#else
    ImageLoadResult load_raw(const std::string &, Header &, ImageArray &) {
      return missing_support("raw camera", "LibRaw");
    }
#endif
  } // namespace imageio
} // namespace astap
