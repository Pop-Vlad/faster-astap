// Rice decompression for FITS tiled image compression (ZCMPTYPE = 'RICE_1'),
// the format fpack writes and astap_cli reads from .fz files.
//
// Ported from unit_ricecomp_unthreaded.pas, which is itself a conversion of the
// decompression half of ricecomp.c from NASA's CFITSIO. The original code was
// written by Richard White at STScI and made available for use in CFITSIO in
// July 1999, under the CFITSIO licence:
//
//   Copyright (Unpublished--all rights reserved under the copyright laws of
//   the United States), U.S. Government as represented by the Administrator
//   of the National Aeronautics and Space Administration. No copyright is
//   claimed in the United States under Title 17, U.S. Code.
//
//   Permission to freely use, copy, modify, and distribute this software
//   and its documentation without fee is hereby granted, provided that this
//   copyright notice and disclaimer of warranty appears in all copies.
//
// The algorithm depends on unsigned 32 bit wraparound, which uint32_t provides.

#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "astap/types.h"

namespace astap {
  // How the tiles of the primary data column are encoded.
  enum class TileCodec {
    rice,  // ZCMPTYPE = 'RICE_1'
    gzip1, // ZCMPTYPE = 'GZIP_1', the byte stream of the values
    gzip2, // ZCMPTYPE = 'GZIP_2', the same with the value bytes shuffled
  };

  // True when this build can inflate GZIP tiles. Even a RICE_1 file needs it
  // for the tiles fpack could not compress, which it stores in a second column
  // as GZIP instead. astap_cli skips those tiles and reports a warning.
  bool gzip_tiles_available();

  // Everything the tile loop needs, filled by load_fits from the BINTABLE
  // header. Mirrors Trice_decode_params.
  struct RiceDecodeParams {
    TileCodec codec = TileCodec::rice;
    // Source buffers, owned by the caller.
    const uint8_t *table_buffer = nullptr; // the NAXIS1 * NAXIS2 table rows
    const uint8_t *heap_buffer = nullptr;  // the PCOUNT byte heap, may be null
    long long heap_size = 0;
    int table_rowwidth = 0; // NAXIS1, bytes per row
    int table_rows = 0;     // NAXIS2, number of rows = number of tiles

    // Tile grid and image geometry.
    int tiles_x = 0, tiles_y = 0, tiles_z = 0;
    int total_tiles = 0;
    int ztile1 = 0, ztile2 = 0, ztile3 = 0;
    int znaxis1 = 0, znaxis2 = 0, znaxis3 = 1;
    int img_width = 0, img_height = 0, img_naxis3 = 1;

    // Byte offsets of the columns within one table row, < 0 when absent.
    int off_comp = -1, off_gzip = -1, off_zscale = -1, off_zzero = -1, off_zblank = -1;

    // Decoder and scaling parameters.
    int zbitpix = 0; // data type of the uncompressed image
    int bytepix = 0, blocksize = 32;
    bool zquantiz_is_none = true;
    bool dither_active = false;
    bool dither_is_2 = false;
    int zdither0 = 0;
    double zscale = 1, zzero = 0;
    int zblank = 0;
    bool zblank_present = false;
    double bscale = 1, bzero = 0;
    bool fastpath_possible = false;

    // CFITSIO random table, 10000 entries, or null when not dithering.
    const float *dither_table = nullptr;
  };

  struct RiceDecodeStatus {
    float measured_max = 0;
    float measured_min = 0;
    bool err_gzip = false;   // a GZIP tile was met and this build cannot inflate
    bool err_decode = false; // rice_decode failed on a tile
    bool err_range = false;  // a heap descriptor pointed outside the heap
    int err_tile = -1;       // first tile with a decode or range problem
    std::string err_msg;
  };

  // CFITSIO's 10000 entry random value table (fits_init_randoms), a Park-Miller
  // minimal standard generator seeded with 1. Needed for the dithered
  // quantisation modes.
  void build_dither_table(std::vector<float> &table);

  // Decodes and places every tile of a RICE_1 compressed image into `img`,
  // which must already be sized [naxis3][height][width].
  RiceDecodeStatus rice_decode_tiles(ImageArray &img, const RiceDecodeParams &p);
} // namespace astap
