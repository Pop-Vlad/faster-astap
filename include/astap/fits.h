// Minimal FITS reader for the plate solver.
// Ported from load_fits() in unit_command_line_general.pas, restricted to
// uncompressed images: BITPIX 8, 16, 32, -32 and -64, NAXIS 2 or 3 (and the
// RGB-in-NAXIS1 variant), with BSCALE / BZERO applied.
//
// Rice compressed (.fz) files and the TIFF/PNG/JPEG loaders of the original are
// not part of this port.

#pragma once

#include <string>

#include "astap/types.h"

namespace astap {
  struct FitsLoadResult {
    bool ok = false;
    std::string error;
    // Values picked up from the header that the solver uses as a starting point.
    double ra_mount = 99999; // >= 999 means "not specified"
    double dec_mount = 99999;
    double focallen = 0;
    double equinox = 2000;
  };

  // Loads a FITS file into `img` (indexed [colour][y][x], y=0 is the first row
  // stored in the file) and fills `head`. The raw 80 character header cards are
  // kept in head.cards so the .wcs output can be built from them.
  FitsLoadResult load_fits(const std::string &filename, Header &head, ImageArray &img);

  // Writes the header cards as a FITS style file (80 character records padded to
  // a multiple of 2880 bytes), the ".wcs" output of astap_cli -wcs.
  bool write_fits_header_file(const std::string &filename, const std::vector<std::string> &cards);
} // namespace astap
