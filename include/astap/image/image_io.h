// Image loading, the equivalent of load_image() in unit_command_line_general.pas.
//
// The file type follows from the extension, exactly as in astap_cli:
//
//   .fit .fits .fts .new   uncompressed FITS
//   .fz                    Rice compressed FITS (ZCMPTYPE = 'RICE_1')
//   .ppm .pgm .pfm         P5 / P6 / PF / Pf binary Netpbm and Portable Float Map
//   .bmp                   Windows bitmap
//   .png .jpg .jpeg        via libpng / libjpeg, when the build found them
//   .tif .tiff             via libtiff, when the build found it
//   plus the raw camera formats below via LibRaw, when the build found it
//
// Which of the optional formats a binary supports is a build time property;
// supported_image_extensions() reports it and the loader says so when asked for
// a format it was built without.

#pragma once

#include <string>

#include "astap/types.h"

namespace astap {
  struct ImageLoadResult {
    bool ok = false;
    std::string error;
    std::string warning; // non fatal remark, e.g. an unsupported tile encoding
    // Values picked up from the header that the solver uses as a starting point.
    double ra_mount = 99999; // >= 999 means "not specified"
    double dec_mount = 99999;
    double focallen = 0;
    double equinox = 2000;
  };

  // Loads `filename` into `img` (indexed [colour][y][x]) and fills `head`,
  // including the 80 character header cards the .wcs output is built from.
  ImageLoadResult load_image(const std::string &filename, Header &head, ImageArray &img);

  // Extensions this build accepts, space separated, for the usage text.
  std::string supported_image_extensions();
} // namespace astap
