// FITS reader for the plate solver.
// Ported from load_fits() in unit_command_line_general.pas: BITPIX 8, 16, 32,
// -32 and -64, NAXIS 2 or 3 (and the RGB-in-NAXIS1 variant), with BSCALE /
// BZERO applied, plus Rice compressed images (.fz, ZCMPTYPE = 'RICE_1') stored
// in a BINTABLE extension.
//
// When the primary HDU holds no image (NAXIS = 0) the first extension is used,
// which is how .fz and multi extension files are laid out.

#pragma once

#include <string>

#include "astap/image/image_io.h"
#include "astap/types.h"

namespace astap {
  // Loads a FITS file into `img` (indexed [colour][y][x], y=0 is the first row
  // stored in the file) and fills `head`. The raw 80 character header cards are
  // kept in head.cards so the .wcs output can be built from them.
  ImageLoadResult load_fits(const std::string &filename, Header &head, ImageArray &img);

  // Writes the header cards as a FITS style file (80 character records padded to
  // a multiple of 2880 bytes), the ".wcs" output of astap_cli -wcs.
  bool write_fits_header_file(const std::string &filename, const std::vector<std::string> &cards);
} // namespace astap
