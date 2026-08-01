#include "astap/image/rice.h"

#include <algorithm>
#include <cmath>
#include <cstring>

#ifdef ASTAP_WITH_ZLIB
#include <zlib.h>
#endif

namespace astap {
  namespace {
    // How many bytes past the end of a compressed tile the decoder may touch.
    // The heap is allocated with this much zero padding, so the hard limit
    // below is always inside the caller's allocation.
    constexpr int kInputPadding = 16;

    // Number of bits in an 8 bit value, leading zeros excluded.
    const int *nonzero_count_table() {
      static int t[256];
      static bool init = false;
      if (!init) {
        t[0] = 0;
        for (int i = 1; i < 256; i++) {
          int n = 0;
          for (int v = i; v; v >>= 1) n++;
          t[i] = n;
        }
        init = true;
      }
      return t;
    }

    // (1 << n) - 1, defined for n up to 32.
    inline uint32_t low_mask(int n) {
      return n >= 32 ? 0xFFFFFFFFu : ((static_cast<uint32_t>(1) << n) - 1);
    }

    // fits_rdecomp / fits_rdecomp_short / fits_rdecomp_byte in one function: the
    // three differ only in the constants derived from the pixel size and in the
    // width the running value is truncated to. Returns 0 on success, 1 on
    // failure, exactly like the C original.
    //
    // Unlike CFITSIO this checks every read against the end of the buffer.
    // CFITSIO relies on the caller over-allocating the input; the checks cost
    // nothing measurable and stop a corrupt stream from running away.
    template<typename T>
    int rdecomp(const uint8_t *c, int clen, T *arr, int nx, int nblock) {
      const int bsize = static_cast<int>(sizeof(T));
      const int fsbits = bsize == 1 ? 3 : (bsize == 2 ? 4 : 5);
      const int fsmax = bsize == 1 ? 6 : (bsize == 2 ? 14 : 25);
      const int bbits = 1 << fsbits;
      const int *nonzero_count = nonzero_count_table();

      if (clen < bsize) return 1;

      // The first pixel is stored raw, big endian, and is not encoded.
      uint32_t lastpix = 0;
      for (int i = 0; i < bsize; i++) lastpix = (lastpix << 8) | c[i];

      const uint8_t *const cstart = c;
      c += bsize;
      const uint8_t *const cend = cstart + clen;     // CFITSIO's per block check
      const uint8_t *const chard = cstart + clen + kInputPadding; // hard limit

      uint32_t b = *c++; // bit buffer
      int nbits = 8;     // bits remaining in b
      int i = 0;
      while (i < nx) {
        // The FS value is in the first fsbits bits.
        nbits -= fsbits;
        while (nbits < 0) {
          if (c >= chard) return 1;
          b = (b << 8) | *c++;
          nbits += 8;
        }
        const int fs = static_cast<int>(b >> nbits) - 1;

        b &= low_mask(nbits);
        int imax = i + nblock;
        if (imax > nx) imax = nx;

        if (fs < 0) {
          // Low entropy case, all differences are zero.
          while (i < imax) arr[i++] = static_cast<T>(lastpix);
        } else if (fs == fsmax) {
          // High entropy case, the pixel values are coded directly.
          while (i < imax) {
            int k = bbits - nbits;
            uint32_t diff = k >= 32 ? 0u : (b << k); // b is 0 when k is 32
            k -= 8;
            while (k >= 0) {
              if (c >= chard) return 1;
              b = *c++;
              diff |= b << k;
              k -= 8;
            }
            if (nbits > 0) {
              if (c >= chard) return 1;
              b = *c++;
              diff |= b >> (-k);
              b &= low_mask(nbits);
            } else {
              b = 0;
            }

            // Undo the mapping and the differencing. These operations overflow
            // the unsigned arithmetic on purpose.
            diff = (diff >> 1) ^ (0u - (diff & 1));
            arr[i] = static_cast<T>(diff + lastpix);
            lastpix = arr[i];
            i++;
          }
        } else {
          // Normal case, Rice coding.
          while (i < imax) {
            while (b == 0) { // count the leading zeros
              if (c >= chard) return 1;
              nbits += 8;
              b = *c++;
            }
            const int nzero = nbits - nonzero_count[b];
            nbits -= nzero + 1;
            b ^= static_cast<uint32_t>(1) << nbits; // flip the leading one bit
            nbits -= fs;                            // the fs trailing bits
            while (nbits < 0) {
              if (c >= chard) return 1;
              b = (b << 8) | *c++;
              nbits += 8;
            }
            uint32_t diff = (static_cast<uint32_t>(nzero) << fs) | (b >> nbits);
            b &= low_mask(nbits);

            diff = (diff >> 1) ^ (0u - (diff & 1));
            arr[i] = static_cast<T>(diff + lastpix);
            lastpix = arr[i];
            i++;
          }
        }
        if (c > cend) return 1; // hit the end of the compressed byte stream
      }
      return 0;
    }

    // Big endian readers for the table columns.
    inline uint32_t be32(const uint8_t *p) {
      return (static_cast<uint32_t>(p[0]) << 24) | (static_cast<uint32_t>(p[1]) << 16) |
             (static_cast<uint32_t>(p[2]) << 8) | static_cast<uint32_t>(p[3]);
    }

    inline double be_double(const uint8_t *p) {
      uint64_t q = 0;
      for (int i = 0; i < 8; i++) q = (q << 8) | p[i];
      double d;
      std::memcpy(&d, &q, 8);
      return d;
    }

    // Inflates one GZIP tile into exactly `want` bytes. The stream is a gzip
    // member, the format CFITSIO writes for both GZIP_1 tiles and the GZIP
    // fallback of a Rice compressed file.
    bool inflate_tile(const uint8_t *src, size_t len, uint8_t *dst, size_t want) {
#ifdef ASTAP_WITH_ZLIB
      z_stream zs;
      std::memset(&zs, 0, sizeof(zs));
      // 15 window bits plus 32, which lets zlib accept both gzip and zlib wrappers.
      if (inflateInit2(&zs, 15 + 32) != Z_OK) return false;
      zs.next_in = const_cast<Bytef *>(reinterpret_cast<const Bytef *>(src));
      zs.avail_in = static_cast<uInt>(len);
      zs.next_out = reinterpret_cast<Bytef *>(dst);
      zs.avail_out = static_cast<uInt>(want);
      const int rc = inflate(&zs, Z_FINISH);
      const uLong produced = zs.total_out;
      inflateEnd(&zs);
      return (rc == Z_STREAM_END || rc == Z_OK) && produced == want;
#else
      (void) src;
      (void) len;
      (void) dst;
      (void) want;
      return false;
#endif
    }

    // A GZIP tile holds the values as they would sit in an uncompressed FITS:
    // big endian, `bytepix` bytes each. GZIP_2 additionally stores all the most
    // significant bytes first, then all the second ones, and so on.
    //
    // The values are rewritten in place in the byte order of this machine,
    // which is what the placement loop reads them back as.
    void unpack_gzip_tile(uint8_t *buf, size_t pixels, int bytepix, bool shuffled,
                          std::vector<uint8_t> &spare) {
      if (bytepix == 1) return;
      const size_t bytes = pixels * static_cast<size_t>(bytepix);
      spare.assign(buf, buf + bytes);
      for (size_t i = 0; i < pixels; i++) {
        uint64_t v = 0;
        for (int b = 0; b < bytepix; b++) // b = 0 is the most significant byte
          v = (v << 8) |
              (shuffled ? spare[static_cast<size_t>(b) * pixels + i] : spare[i * bytepix + b]);
        switch (bytepix) {
          case 2: {
            const uint16_t x = static_cast<uint16_t>(v);
            std::memcpy(buf + i * 2, &x, 2);
            break;
          }
          case 4: {
            const uint32_t x = static_cast<uint32_t>(v);
            std::memcpy(buf + i * 4, &x, 4);
            break;
          }
          default: {
            std::memcpy(buf + i * 8, &v, 8);
            break;
          }
        }
      }
    }
  } // namespace

  bool gzip_tiles_available() {
#ifdef ASTAP_WITH_ZLIB
    return true;
#else
    return false;
#endif
  }

  void build_dither_table(std::vector<float> &table) {
    // Park-Miller minimal standard generator seeded with 1, values in [0, 1).
    const int ia = 16807, im = 2147483647, iq = 127773, ir = 2836;
    table.resize(10000);
    int seed = 1;
    for (int i = 0; i < 10000; i++) {
      const int kk = seed / iq;
      seed = ia * (seed - kk * iq) - ir * kk;
      if (seed < 0) seed += im;
      double temp = static_cast<double>(seed) / im;
      if (temp > 0.9999999) temp = 0.9999999;
      table[i] = static_cast<float>(temp);
    }
  }

  RiceDecodeStatus rice_decode_tiles(ImageArray &img, const RiceDecodeParams &p) {
    RiceDecodeStatus st;

    // The tile grid may be marginally larger than the number of table rows.
    int ntiles = p.total_tiles;
    if (ntiles > p.table_rows) ntiles = p.table_rows;
    if (ntiles <= 0) return st;

    long long max_tile_pixels =
        static_cast<long long>(p.ztile1) * p.ztile2 * p.ztile3;
    if (max_tile_pixels < 1) max_tile_pixels = 1;
    // A tile fpack failed to compress is stored as GZIP of the *original*
    // values, so the scratch buffer has to fit those as well as the encoded
    // ones, which for a quantised float image are narrower.
    const int original_bytes = std::max(1, std::abs(p.zbitpix) / 8);
    const int widest = std::max(p.bytepix, original_bytes);
    std::vector<uint8_t> scratch(static_cast<size_t>(max_tile_pixels) * widest + 16);
    std::vector<uint8_t> spare; // used when a GZIP tile has to be reordered

    for (int tile_index = 0; tile_index < ntiles; tile_index++) {
      const int tx = tile_index % p.tiles_x;
      const int ty = (tile_index / p.tiles_x) % p.tiles_y;
      const int tz = tile_index / (p.tiles_x * p.tiles_y);

      const uint8_t *row =
          p.table_buffer + static_cast<size_t>(tile_index) * p.table_rowwidth;

      // Variable length array descriptor: element count then heap offset.
      const uint8_t *desc = row + p.off_comp;
      const int compressed_len = static_cast<int>(be32(desc));
      const int heap_offset = static_cast<int>(be32(desc + 4));

      const int tile_w = std::min(p.ztile1, p.znaxis1 - tx * p.ztile1);
      const int tile_h = std::min(p.ztile2, p.znaxis2 - ty * p.ztile2);
      const int tile_d = std::min(p.ztile3, p.znaxis3 - tz * p.ztile3);
      const long long tile_pixels = static_cast<long long>(tile_w) * tile_h * tile_d;
      if (tile_pixels <= 0) continue;

      // A tile fpack could not compress is stored in a second column as GZIP
      // instead, and a whole file may be GZIP encoded to begin with.
      bool gzip_tile = p.codec != TileCodec::rice;
      bool fallback_tile = false;
      int length = compressed_len;
      int offset = heap_offset;
      if (compressed_len <= 0) {
        if (p.off_gzip < 0) continue; // an empty tile, nothing to place
        const uint8_t *desc_gzip = row + p.off_gzip;
        length = static_cast<int>(be32(desc_gzip));
        offset = static_cast<int>(be32(desc_gzip + 4));
        if (length <= 0) continue;
        gzip_tile = true;
        fallback_tile = true;
        if (!gzip_tiles_available()) {
          st.err_gzip = true;
          continue;
        }
      }

      // The fallback column holds the values as they were before quantisation,
      // in the type of the uncompressed image; every other tile holds what the
      // encoder produced.
      const int elem = fallback_tile ? original_bytes : p.bytepix;
      const bool tile_floats = p.zbitpix < 0 && (p.zquantiz_is_none || fallback_tile);

      if (offset < 0 || static_cast<long long>(offset) + length > p.heap_size) {
        if (!st.err_range) {
          st.err_range = true;
          st.err_tile = tile_index;
        }
        continue;
      }

      double tile_scale = p.zscale;
      double tile_zero = p.zzero;
      int tile_blank = p.zblank;
      bool tile_has_blank = p.zblank_present;
      if (p.off_zscale >= 0) tile_scale = be_double(row + p.off_zscale);
      if (p.off_zzero >= 0) tile_zero = be_double(row + p.off_zzero);
      if (p.off_zblank >= 0) {
        tile_blank = static_cast<int>(be32(row + p.off_zblank));
        tile_has_blank = true;
      }

      const uint8_t *compressed = p.heap_buffer + offset;
      const int nx = static_cast<int>(tile_pixels);
      int status = 1;
      if (gzip_tile) {
        const size_t want = static_cast<size_t>(nx) * elem;
        if (inflate_tile(compressed, static_cast<size_t>(length), scratch.data(), want)) {
          // Only the primary column of a GZIP_2 file is byte shuffled.
          unpack_gzip_tile(scratch.data(), static_cast<size_t>(nx), elem,
                           p.codec == TileCodec::gzip2 && !fallback_tile, spare);
          status = 0;
        }
      } else {
        switch (p.bytepix) {
          case 1:
            status = rdecomp<uint8_t>(compressed, length, scratch.data(), nx, p.blocksize);
            break;
          case 2:
            status = rdecomp<uint16_t>(compressed, length,
                                       reinterpret_cast<uint16_t *>(scratch.data()), nx,
                                       p.blocksize);
            break;
          case 4:
            status = rdecomp<uint32_t>(compressed, length,
                                       reinterpret_cast<uint32_t *>(scratch.data()), nx,
                                       p.blocksize);
            break;
          default:
            break;
        }
      }
      if (status != 0) {
        if (!st.err_decode) {
          st.err_decode = true;
          st.err_tile = tile_index;
          st.err_msg = p.bytepix == 1 || p.bytepix == 2 || p.bytepix == 4
                         ? "corrupt or truncated compressed byte stream"
                         : "BYTEPIX must be 1, 2 or 4 bytes";
        }
        continue;
      }

      // Fast path: lossless 16 bit tile spanning a full image row.
      if (p.fastpath_possible && tile_w == p.img_width && ty < p.img_height) {
        const uint16_t *src = reinterpret_cast<const uint16_t *>(scratch.data());
        float *dst = img.row(0, ty); // znaxis3 is 1 on this path
        for (int x = 0; x < p.img_width; x++) {
          const float v =
              static_cast<float>(static_cast<int16_t>(src[x]) * p.bscale + p.bzero);
          dst[x] = v;
          if (v > st.measured_max) st.measured_max = v;
          else if (v < st.measured_min) st.measured_min = v;
        }
        continue;
      }

      int dither_iseed = 0, dither_next = 0;
      if (p.dither_active) {
        dither_iseed = (tile_index + p.zdither0 - 1) % 10000;
        dither_next = static_cast<int>(p.dither_table[dither_iseed] * 500);
      }

      for (int pz = 0; pz < tile_d; pz++)
        for (int py = 0; py < tile_h; py++)
          for (int px = 0; px < tile_w; px++) {
            const size_t idx = static_cast<size_t>(px) + static_cast<size_t>(py) * tile_w +
                               static_cast<size_t>(pz) * tile_w * tile_h;
            const int img_x = tx * p.ztile1 + px;
            const int img_y = ty * p.ztile2 + py;
            const int img_z = tz * p.ztile3 + pz;
            const bool store =
                img_z < p.img_naxis3 && img_y < p.img_height && img_x < p.img_width;

            long long value = 0;
            double float_value = 0;
            if (tile_floats) {
              // The tile holds the pixels themselves rather than quantised
              // integers.
              if (elem == 8) {
                std::memcpy(&float_value, scratch.data() + idx * 8, 8);
              } else {
                float f;
                std::memcpy(&f, scratch.data() + idx * 4, 4);
                float_value = f;
              }
            } else {
              switch (elem) {
                case 1:
                  value = scratch[idx];
                  break;
                case 2:
                  value =
                      static_cast<int16_t>(reinterpret_cast<const uint16_t *>(scratch.data())[idx]);
                  break;
                default: // 4
                  value =
                      static_cast<int32_t>(reinterpret_cast<const uint32_t *>(scratch.data())[idx]);
                  break;
              }
            }

            double col;
            if (tile_floats) {
              col = float_value * p.bscale + p.bzero;
              if (std::isnan(col) || std::isinf(col)) col = 0;
            } else if (!p.zquantiz_is_none) {
              if (value == -2147483647LL)
                col = 0;
              else if (p.dither_is_2 && value == -2147483646LL)
                col = 0.0;
              else if (p.dither_active)
                // Both casts matter. The dither table is single precision, and
                // C++ would evaluate `value - table[i]` in float, which drops
                // the low bits of an integer near 2^31. SUBTRACTIVE_DITHER_2
                // puts every value there, so in float the image comes back
                // quantised to steps of 128 * ZSCALE. FPC evaluates the same
                // expression in extended precision, hence the doubles here.
                col = (static_cast<double>(value) - static_cast<double>(p.dither_table[dither_next]) +
                       0.5) *
                        tile_scale +
                      tile_zero;
              else
                col = value * tile_scale + tile_zero;
              if (std::isnan(col) || std::isinf(col)) col = 0;
            } else {
              if (tile_has_blank && value == tile_blank)
                col = 0;
              else
                col = value * p.bscale + p.bzero;
            }

            if (store) {
              const float v = static_cast<float>(col);
              img.at(img_z, img_y, img_x) = v;
              if (v > st.measured_max) st.measured_max = v;
              else if (v < st.measured_min) st.measured_min = v;
            }

            if (p.dither_active) {
              dither_next++;
              if (dither_next >= 10000) {
                dither_iseed++;
                if (dither_iseed >= 10000) dither_iseed = 0;
                dither_next = static_cast<int>(p.dither_table[dither_iseed] * 500);
              }
            }
          }
    }
    return st;
  }
} // namespace astap
