// Per format loaders behind load_image(). Internal to the library.
//
// Every loader fills `img` as [colour][y][x] with row 0 at the bottom of the
// image (the FITS convention) and calls synthesise_header() so that the .ini
// and .wcs output has a header to work from, exactly as astap_cli builds one
// for a non FITS file.

#pragma once

#include <string>

#include "astap/image/image_io.h"
#include "astap/types.h"

namespace astap {
  namespace imageio {
    ImageLoadResult load_pnm(const std::string &filename, Header &head, ImageArray &img);

    ImageLoadResult load_bmp(const std::string &filename, Header &head, ImageArray &img);

    ImageLoadResult load_png(const std::string &filename, Header &head, ImageArray &img);

    ImageLoadResult load_jpeg(const std::string &filename, Header &head, ImageArray &img);

    ImageLoadResult load_tiff(const std::string &filename, Header &head, ImageArray &img);

    ImageLoadResult load_raw(const std::string &filename, Header &head, ImageArray &img);

    // The minimal FITS header astap_cli synthesises for an image that carries
    // none of its own, in the same keyword order.
    void synthesise_header(Header &head, int bitpix, int naxis3, double datamax);

    // Error for a format this build has no library for.
    ImageLoadResult missing_support(const std::string &format, const std::string &library);
  } // namespace imageio
} // namespace astap
