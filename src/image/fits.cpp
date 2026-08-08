#include "astap/image/fits.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <vector>

#include "astap/astro_math.h"
#include "astap/parallel.h"
#include "astap/image/rice.h"

namespace astap {
  namespace {
    // Value field of a fixed format card, positions 11..30 (1 based).
    bool card_value_double(const std::string &card, double &out) {
      if (card.size() < 11) return false;
      std::string t;
      for (size_t r = 10; r < card.size() && r < 31; r++) {
        // The '/' check is not strictly necessary but safer. CFITSIO may write
        // negative values up to position 31, a violation of the FITS standard.
        if (card[r] == '/') break;
        if (card[r] != ' ') t += card[r];
      }
      if (t.empty()) return false;
      try {
        size_t used = 0;
        out = std::stod(t, &used);
        return used == t.size();
      } catch (...) {
        return false;
      }
    }

    // Quoted string value of a card. Single quotes start at position 11 for fixed
    // format cards (FITS standard 4.0, chapter 4.2.1.1).
    std::string card_value_string(const std::string &card) {
      std::string r;
      for (size_t i = 11; i < card.size() && i < 79; i++) {
        if (card[i] == '\'') break;
        r += card[i];
      }
      // trim
      size_t b = r.find_first_not_of(' ');
      size_t e = r.find_last_not_of(' ');
      if (b == std::string::npos) return std::string();
      return r.substr(b, e - b + 1);
    }

    bool key_is(const std::string &card, const char *key) {
      size_t n = std::strlen(key);
      return card.size() >= n && card.compare(0, n, key) == 0;
    }

    // Big endian readers. FITS data is always stored most significant byte first.
    inline uint16_t be16(const uint8_t *p) {
      return static_cast<uint16_t>((static_cast<uint16_t>(p[0]) << 8) | p[1]);
    }

    inline uint32_t be32(const uint8_t *p) {
      return (static_cast<uint32_t>(p[0]) << 24) | (static_cast<uint32_t>(p[1]) << 16) |
             (static_cast<uint32_t>(p[2]) << 8) | static_cast<uint32_t>(p[3]);
    }

    inline uint64_t be64(const uint8_t *p) {
      uint64_t v = 0;
      for (int i = 0; i < 8; i++) v = (v << 8) | p[i];
      return v;
    }

    // Width in bytes of one TFORM entry of a BINTABLE. A variable length array
    // descriptor ('1PB(1234)' and friends) is 8 bytes, 16 for the Q form.
    int tform_width(const std::string &tf) {
      size_t b = tf.find_first_not_of(' ');
      if (b == std::string::npos) return 0;
      size_t e = tf.find_last_not_of(' ');
      const std::string t = tf.substr(b, e - b + 1);

      size_t n = 0;
      while (n < t.size() && t[n] >= '0' && t[n] <= '9') n++;
      int rep = 1;
      if (n > 0) {
        try {
          rep = std::stoi(t.substr(0, n));
        } catch (...) {
          rep = 1;
        }
      }
      if (rep < 0) rep = 0;
      if (n >= t.size()) return 0;

      switch (std::toupper(static_cast<unsigned char>(t[n]))) {
        case 'P': return 8;
        case 'Q': return 16;
        case 'L': return 1 * rep;
        case 'X': return (rep + 7) / 8;
        case 'B': return 1 * rep;
        case 'I': return 2 * rep;
        case 'J': return 4 * rep;
        case 'K': return 8 * rep;
        case 'A': return 1 * rep;
        case 'E': return 4 * rep;
        case 'D': return 8 * rep;
        case 'C': return 8 * rep;
        case 'M': return 16 * rep;
        default: return rep;
      }
    }

    // Index of a TFORMnnn / TTYPEnnn card, 1 based as written in the card, or 0.
    int keyword_index(const std::string &card, size_t first_digit) {
      std::string digits;
      for (size_t i = first_digit; i < card.size() && i < first_digit + 3; i++)
        if (card[i] != ' ') digits += card[i];
      if (digits.empty()) return 0;
      try {
        return std::stoi(digits);
      } catch (...) {
        return 0;
      }
    }
  } // namespace

  ImageLoadResult load_fits(const std::string &filename, Header &head, ImageArray &img) {
    ImageLoadResult res;

    std::ifstream f(filename, std::ios::binary);
    if (!f.is_open()) {
      res.error = "Error, accessing the file!";
      return res;
    }

    head = Header();
    img = ImageArray();

    double bscale = 1;
    double bzero = 0;
    int naxis1 = 0;
    std::string objctra, objctdec;

    bool simple = false;
    bool end_found = false;
    char block[2880];

    // Tiled image compression (.fz) keywords. The image lives in a BINTABLE
    // extension whose header carries the geometry of the uncompressed image.
    bool zimage = false;
    std::string zcmptype, zquantiz = "NONE";
    int zbitpix = 0, znaxis1 = 0, znaxis2 = 0, znaxis3 = 1;
    int ztile1 = 0, ztile2 = 0, ztile3 = 0;
    int zdither0 = 0;
    double zscale = 1, zzero = 0;
    int zblank = 0;
    bool zblank_present = false;
    long long pcount = 0;
    int blocksize = 32, bytepix = 0;
    int tfields = 0;
    std::vector<std::string> tform, ttype;
    std::string zname[4];

    // Read the header, 2880 byte blocks of 36 cards. When the primary HDU holds
    // no image (NAXIS = 0) the first extension is used instead, which is how
    // compressed .fz files and multi extension files are laid out.
    while (true) {
      end_found = false;
      while (!end_found) {
        f.read(block, 2880);
        if (f.gcount() != 2880) {
          res.error = "Error, unexpected end of file in the header!";
          return res;
        }
        if (!simple) {
          if (std::strncmp(block, "SIMPLE  =", 9) != 0) {
            res.error = "Error, not a FITS file (no SIMPLE keyword)!";
            return res;
          }
          simple = true;
        }

        for (int c = 0; c < 36 && !end_found; c++) {
          std::string card(block + c * 80, 80);
          head.cards.push_back(card);

          double v;
          if (key_is(card, "END ")) {
            end_found = true;
          } else if (key_is(card, "NAXIS   ")) {
            if (card_value_double(card, v)) head.naxis = static_cast<int>(pround(v));
          } else if (key_is(card, "NAXIS1  ")) {
            if (card_value_double(card, v)) {
              naxis1 = static_cast<int>(pround(v));
              head.width = naxis1;
            }
          } else if (key_is(card, "NAXIS2  ")) {
            if (card_value_double(card, v)) head.height = static_cast<int>(pround(v));
          } else if (key_is(card, "NAXIS3  ")) {
            if (card_value_double(card, v)) head.naxis3 = static_cast<int>(pround(v));
          } else if (key_is(card, "BITPIX  ")) {
            if (card_value_double(card, v)) head.bitpix = static_cast<int>(pround(v));
          } else if (key_is(card, "BZERO   ")) {
            if (card_value_double(card, v)) bzero = v;
          } else if (key_is(card, "BSCALE  ")) {
            if (card_value_double(card, v)) bscale = v; // rarely used, normally 1
          } else if (key_is(card, "CDELT1  ")) {
            if (card_value_double(card, v)) head.cdelt1 = v; // deg/pixel for RA
          } else if (key_is(card, "CDELT2  ")) {
            if (card_value_double(card, v)) head.cdelt2 = v; // deg/pixel for DEC
          } else if (key_is(card, "CD1_1   ")) {
            if (card_value_double(card, v)) head.cd1_1 = v;
          } else if (key_is(card, "CD1_2   ")) {
            if (card_value_double(card, v)) head.cd1_2 = v;
          } else if (key_is(card, "CD2_1   ")) {
            if (card_value_double(card, v)) head.cd2_1 = v;
          } else if (key_is(card, "CD2_2   ")) {
            if (card_value_double(card, v)) head.cd2_2 = v;
          } else if (key_is(card, "CRVAL1  ")) {
            if (card_value_double(card, v)) head.ra0 = v * kPi / 180;
          } else if (key_is(card, "CRVAL2  ")) {
            if (card_value_double(card, v)) head.dec0 = v * kPi / 180;
          } else if (key_is(card, "SECPIX") || key_is(card, "SCALE   ") || key_is(card, "PIXSCALE")) {
            // No CDELT found yet, use the alternative.
            if (head.cdelt2 == 0 && card_value_double(card, v)) {
              head.cdelt2 = v / 3600;
              head.cdelt1 = head.cdelt2;
            }
          } else if (key_is(card, "EQUINOX ")) {
            if (card_value_double(card, v)) res.equinox = v;
          } else if (key_is(card, "FOCALLEN")) {
            if (card_value_double(card, v)) res.focallen = v; // in mm, MaximDL keyword
          } else if (key_is(card, "DEC     ")) {
            if (card_value_double(card, v)) {
              res.dec_mount = v * kPi / 180;
              if (head.dec0 == 0) head.dec0 = res.dec_mount; // only when CRVAL2 is missing
            } else {
              objctdec = card_value_string(card);
            }
          } else if (key_is(card, "RA      ")) {
            if (card_value_double(card, v)) {
              res.ra_mount = v * kPi / 180;
              if (head.ra0 == 0) head.ra0 = res.ra_mount; // only when CRVAL1 is missing
            } else {
              objctra = card_value_string(card);
            }
          } else if (key_is(card, "OBJCTRA ")) {
            if (res.ra_mount >= 999) objctra = card_value_string(card); // preference for RA
          } else if (key_is(card, "OBJCTDEC")) {
            if (res.dec_mount >= 999) objctdec = card_value_string(card); // preference for DEC
          } else if (key_is(card, "XPIXSZ  ")) {
            if (card_value_double(card, v)) head.xpixsz = v; // microns, after binning
          } else if (key_is(card, "YPIXSZ  ")) {
            if (card_value_double(card, v)) head.ypixsz = v;
          } else if (key_is(card, "XBINNING")) {
            if (card_value_double(card, v)) head.xbinning = pround(v);
          } else if (key_is(card, "YBINNING")) {
            if (card_value_double(card, v)) head.ybinning = pround(v);
          } else if (key_is(card, "PCOUNT  ")) {
            if (card_value_double(card, v)) pcount = static_cast<long long>(pround(v));
          } else if (key_is(card, "TFIELDS ")) {
            if (card_value_double(card, v)) {
              tfields = static_cast<int>(pround(v));
              if (tfields > 0) {
                tform.assign(static_cast<size_t>(tfields), std::string());
                ttype.assign(static_cast<size_t>(tfields), std::string());
              }
            }
          } else if (key_is(card, "TFORM") && tfields > 0) {
            const int n = keyword_index(card, 5);
            if (n >= 1 && n <= tfields) tform[static_cast<size_t>(n - 1)] = card_value_string(card);
          } else if (key_is(card, "TTYPE") && tfields > 0) {
            const int n = keyword_index(card, 5);
            if (n >= 1 && n <= tfields) ttype[static_cast<size_t>(n - 1)] = card_value_string(card);
          } else if (key_is(card, "ZIMAGE  ")) {
            // Logical card, T when the extension holds a compressed image.
            zimage = card.find('T', 10) != std::string::npos && card.find('T', 10) < 31;
          } else if (key_is(card, "ZBITPIX ")) {
            if (card_value_double(card, v)) zbitpix = static_cast<int>(pround(v));
          } else if (key_is(card, "ZNAXIS1 ")) {
            if (card_value_double(card, v)) znaxis1 = static_cast<int>(pround(v));
          } else if (key_is(card, "ZNAXIS2 ")) {
            if (card_value_double(card, v)) znaxis2 = static_cast<int>(pround(v));
          } else if (key_is(card, "ZNAXIS3 ")) {
            if (card_value_double(card, v)) znaxis3 = static_cast<int>(pround(v));
          } else if (key_is(card, "ZTILE1  ")) {
            if (card_value_double(card, v)) ztile1 = static_cast<int>(pround(v));
          } else if (key_is(card, "ZTILE2  ")) {
            if (card_value_double(card, v)) ztile2 = static_cast<int>(pround(v));
          } else if (key_is(card, "ZTILE3  ")) {
            if (card_value_double(card, v)) ztile3 = static_cast<int>(pround(v));
          } else if (key_is(card, "ZCMPTYPE")) {
            zcmptype = card_value_string(card);
          } else if (key_is(card, "ZQUANTIZ")) {
            zquantiz = card_value_string(card);
          } else if (key_is(card, "ZDITHER0")) {
            if (card_value_double(card, v)) zdither0 = static_cast<int>(pround(v));
          } else if (key_is(card, "ZSCALE  ")) {
            card_value_double(card, zscale); // header level default, may be per tile
          } else if (key_is(card, "ZZERO   ")) {
            card_value_double(card, zzero);
          } else if (key_is(card, "ZBLANK  ")) {
            if (card_value_double(card, v)) {
              zblank = static_cast<int>(pround(v));
              zblank_present = true;
            }
          } else if (key_is(card, "ZNAME")) {
            const int n = keyword_index(card, 5);
            if (n >= 1 && n <= 4) zname[n - 1] = card_value_string(card);
          } else if (key_is(card, "ZVAL")) {
            // ZNAMEn / ZVALn pairs carry BLOCKSIZE and BYTEPIX.
            const int n = keyword_index(card, 4);
            if (n >= 1 && n <= 4 && card_value_double(card, v)) {
              if (zname[n - 1] == "BLOCKSIZE") blocksize = static_cast<int>(pround(v));
              else if (zname[n - 1] == "BYTEPIX") bytepix = static_cast<int>(pround(v));
            }
          }
        }
      }

      // An empty primary HDU means the image is in the next one.
      if (head.naxis != 0 || f.peek() == std::char_traits<char>::eof()) break;

      head.cards.clear();
      head.naxis3 = 1;
      naxis1 = 0;
      head.width = 0;
      head.height = 0;
      bscale = 1;
      bzero = 0;
      zimage = false;
      zcmptype.clear();
      zquantiz = "NONE";
      zbitpix = 0;
      znaxis1 = znaxis2 = 0;
      znaxis3 = 1;
      ztile1 = ztile2 = ztile3 = 0;
      zdither0 = 0;
      zscale = 1;
      zzero = 0;
      zblank = 0;
      zblank_present = false;
      pcount = 0;
      blocksize = 32;
      bytepix = 0;
      tfields = 0;
      tform.clear();
      ttype.clear();
      for (std::string &z: zname) z.clear();
    }

    if (head.naxis < 2) {
      res.error = "Error, the file contains no image!";
      return res;
    }

    // RGB FITS with NAXIS1=3 is treated as 24 bit coded pixels in 2 dimensions.
    bool rgb24 = (head.naxis == 3 && naxis1 == 3);
    if (rgb24) {
      head.width = head.height;
      head.height = head.naxis3;
      head.naxis3 = 3;
    }
    if (head.naxis == 2) head.naxis3 = 1;
    if (head.naxis3 > 3) head.naxis3 = 1; // more than three colours, use only the first

    if (head.ra0 != 0 || head.dec0 != 0) {
      if (res.equinox != 2000) {
        // e.g. SharpCap
        precession_jnow_to_j2000(res.equinox, head.ra0, head.dec0);
        if (res.dec_mount < 999) precession_jnow_to_j2000(res.equinox, res.ra_mount, res.dec_mount);
      }
    } else if (!objctra.empty()) {
      ra_text_to_radians(objctra, head.ra0);
      dec_text_to_radians(objctdec, head.dec0);
    }

    if (head.cdelt2 == 0) {
      if (head.cd1_1 == 0) {
        // no scale, try to fix it
        // XPIXSZ is including binning (MaximDL keyword).
        if (res.focallen != 0 && head.xpixsz != 0)
          head.cdelt2 = 180 / (kPi * 1000) * head.xpixsz / res.focallen;
      } else {
        head.cdelt2 = std::sqrt(sqr(head.cd1_2) + sqr(head.cd2_2));
      }
    }

    // Tiled image compression (.fz). The BINTABLE holds one compressed tile per
    // row; the pixels are in the heap that follows the rows.
    TileCodec codec = TileCodec::rice;
    bool compressed_image = zimage && zcmptype == "RICE_1";
    if (zimage && (zcmptype == "GZIP_1" || zcmptype == "GZIP_2")) {
      // astap_cli refuses these; this build reads them when it has zlib, which
      // is the same code path the GZIP fallback tiles of a Rice file need.
      if (!gzip_tiles_available()) {
        res.error = "GZIP compression not supported by this build!";
        return res;
      }
      codec = zcmptype == "GZIP_2" ? TileCodec::gzip2 : TileCodec::gzip1;
      compressed_image = true;
    }
    if (zimage && !compressed_image && !zcmptype.empty()) {
      res.error = zcmptype + " compression not supported! Only RICE_1 and GZIP are implemented.";
      return res;
    }

    if (compressed_image) {
      // Geometry of the table itself, before head is given the dimensions of
      // the uncompressed image.
      const int table_rowwidth = naxis1;
      const int table_rows = head.height;

      head.bitpix = zbitpix;
      head.width = znaxis1;
      head.height = znaxis2;
      head.naxis3 = znaxis3;
      head.naxis = head.naxis3 > 1 ? 3 : 2;
      if (head.naxis3 > 3) {
        head.naxis3 = 1; // more than three colours, use only the first
        res.warning = "More than three colours in the compressed image.";
      }
      if (ztile1 <= 0) ztile1 = znaxis1;
      if (ztile2 <= 0) ztile2 = 1; // row by row is the FITS default
      if (ztile3 <= 0) ztile3 = 1;
      if (ztile1 <= 0) ztile1 = 1;

      if (head.naxis < 2 || head.width <= 0 || head.height <= 0) {
        res.error = "Error, invalid dimensions in compressed FITS.";
        return res;
      }
      if (table_rows <= 0 || table_rowwidth <= 0) {
        res.error = "Error, compressed BINTABLE has no rows.";
        return res;
      }

      // Byte offset of every column within a table row, and the columns used.
      int col_comp = -1, col_gzip = -1, col_zscale = -1, col_zzero = -1, col_zblank = -1;
      int run_off = 0;
      std::vector<int> col_offsets(static_cast<size_t>(tfields), 0);
      for (int k = 0; k < tfields; k++) {
        col_offsets[static_cast<size_t>(k)] = run_off;
        run_off += tform_width(tform[static_cast<size_t>(k)]);
        const std::string &name = ttype[static_cast<size_t>(k)];
        if (name == "COMPRESSED_DATA") col_comp = k;
        else if (name == "GZIP_COMPRESSED_DATA") col_gzip = k;
        else if (name == "ZSCALE") col_zscale = k;
        else if (name == "ZZERO") col_zzero = k;
        else if (name == "ZBLANK") col_zblank = k;
      }
      if (run_off != table_rowwidth) {
        res.error = "Error, computed BINTABLE row width (" + std::to_string(run_off) +
                    ") does not match NAXIS1 (" + std::to_string(table_rowwidth) + ").";
        return res;
      }
      if (col_comp < 0) {
        res.error = "Error, COMPRESSED_DATA column not found in compressed FITS.";
        return res;
      }

      if (bytepix == 0) {
        // Not given as a ZNAMEn / ZVALn pair, derive it from ZBITPIX. Rice
        // handles 1, 2 and 4 byte pixels; a GZIP tile may also hold 8 byte
        // doubles, since it stores the values rather than encoding them.
        switch (std::abs(zbitpix)) {
          case 8: bytepix = 1;
            break;
          case 16: bytepix = 2;
            break;
          case 64: bytepix = codec == TileCodec::rice ? 4 : 8;
            break;
          default: bytepix = 4;
            break;
        }
      }
      const bool bytepix_ok =
          bytepix == 1 || bytepix == 2 || bytepix == 4 || (bytepix == 8 && codec != TileCodec::rice);
      if (!bytepix_ok) {
        res.error = "Error, unsupported BYTEPIX " + std::to_string(bytepix) + " in compressed FITS.";
        return res;
      }
      if (blocksize <= 0) blocksize = 32;

      // The table rows first, then the heap. Both follow the header directly.
      const size_t table_size = static_cast<size_t>(table_rowwidth) * table_rows;
      const size_t heap_size = pcount > 0 ? static_cast<size_t>(pcount) : 0;
      std::vector<uint8_t> table(table_size);
      f.read(reinterpret_cast<char *>(table.data()), static_cast<std::streamsize>(table_size));
      if (static_cast<size_t>(f.gcount()) != table_size) {
        res.error = "Error reading compressed table data!";
        return res;
      }
      // The decoder may read a few bytes past the end of the last tile, exactly
      // as the C original does, so the heap keeps some zero padding.
      std::vector<uint8_t> heap(heap_size + 16, 0);
      if (heap_size > 0) {
        f.read(reinterpret_cast<char *>(heap.data()), static_cast<std::streamsize>(heap_size));
        if (static_cast<size_t>(f.gcount()) != heap_size) {
          res.error = "Error reading compressed heap data!";
          return res;
        }
      }

      img.resize(head.naxis3, head.height, head.width);

      RiceDecodeParams p;
      p.codec = codec;
      p.zbitpix = zbitpix;
      p.table_buffer = table.data();
      p.heap_buffer = heap.data();
      p.heap_size = static_cast<long long>(heap_size);
      p.table_rowwidth = table_rowwidth;
      p.table_rows = table_rows;
      p.tiles_x = (znaxis1 + ztile1 - 1) / ztile1;
      p.tiles_y = (znaxis2 + ztile2 - 1) / ztile2;
      p.tiles_z = std::max(1, (znaxis3 + ztile3 - 1) / ztile3);
      p.total_tiles = p.tiles_x * p.tiles_y * p.tiles_z;
      p.ztile1 = ztile1;
      p.ztile2 = ztile2;
      p.ztile3 = ztile3;
      p.znaxis1 = znaxis1;
      p.znaxis2 = znaxis2;
      p.znaxis3 = znaxis3;
      p.img_width = head.width;
      p.img_height = head.height;
      p.img_naxis3 = head.naxis3;
      p.off_comp = col_offsets[static_cast<size_t>(col_comp)];
      p.off_gzip = col_gzip >= 0 ? col_offsets[static_cast<size_t>(col_gzip)] : -1;
      p.off_zscale = col_zscale >= 0 ? col_offsets[static_cast<size_t>(col_zscale)] : -1;
      p.off_zzero = col_zzero >= 0 ? col_offsets[static_cast<size_t>(col_zzero)] : -1;
      p.off_zblank = col_zblank >= 0 ? col_offsets[static_cast<size_t>(col_zblank)] : -1;
      p.bytepix = bytepix;
      p.blocksize = blocksize;
      p.zquantiz_is_none = (zquantiz == "NONE");
      p.dither_is_2 = (zquantiz == "SUBTRACTIVE_DITHER_2");
      p.dither_active =
          (zquantiz == "SUBTRACTIVE_DITHER_1" || p.dither_is_2) && zdither0 > 0;
      p.zdither0 = zdither0;
      p.zscale = zscale;
      p.zzero = zzero;
      p.zblank = zblank;
      p.zblank_present = zblank_present;
      p.bscale = bscale;
      p.bzero = bzero;
      // Full width, single row, lossless 16 bit tiles are the common case and
      // are placed straight into the image without the per pixel branching.
      p.fastpath_possible = bytepix == 2 && p.zquantiz_is_none && !p.dither_active &&
                            !zblank_present && col_zblank < 0 && ztile2 == 1 && ztile3 == 1 &&
                            ztile1 == znaxis1 && znaxis3 == 1;

      std::vector<float> dither_table;
      if (p.dither_active) {
        build_dither_table(dither_table);
        p.dither_table = dither_table.data();
      }

      const RiceDecodeStatus st = rice_decode_tiles(img, p);
      if (st.err_gzip) res.warning = "GZIP tile(s) skipped, this build has no zlib!";
      if (st.err_range)
        res.warning = "Tile " + std::to_string(st.err_tile) + " heap offset out of range, skipped.";
      if (st.err_decode)
        res.warning =
            "Decompression error for tile " + std::to_string(st.err_tile) + ": " + st.err_msg;

      // Same rescaling as the uncompressed path.
      if (zbitpix <= -32 || zbitpix == 32) {
        double scalefactor = 1;
        if (st.measured_max > 0 && (st.measured_max <= 1.5 || st.measured_max > 65535 * 1.5))
          scalefactor = 65535 / st.measured_max;
        if (scalefactor != 1) {
          for (int k = 0; k < head.naxis3; k++)
            for (int j = 0; j < head.height; j++) {
              float *r = img.row(k, j);
              for (int i = 0; i < head.width; i++) r[i] = static_cast<float>(r[i] * scalefactor);
            }
          head.datamax_org = 65535;
        } else {
          head.datamax_org = st.measured_max;
        }
      } else if (zbitpix == 8) {
        head.datamax_org = 255;
      } else {
        head.datamax_org = st.measured_max;
      }

      res.ok = true;
      return res;
    }

    // Read the image data.
    const int bytes_per_pixel = rgb24 ? 3 : std::abs(head.bitpix) / 8;
    if (bytes_per_pixel <= 0) {
      res.error = "Error, unsupported BITPIX!";
      return res;
    }

    img.resize(head.naxis3, head.height, head.width);

    float measured_max = 0;

    if (rgb24) {
      std::vector<uint8_t> line(static_cast<size_t>(head.width) * 3);
      for (int j = 0; j < head.height; j++) {
        f.read(reinterpret_cast<char *>(line.data()), static_cast<std::streamsize>(line.size()));
        if (f.gcount() != static_cast<std::streamsize>(line.size())) break;
        for (int i = 0; i < head.width; i++) {
          img.at(0, j, i) = line[static_cast<size_t>(i) * 3 + 0];
          img.at(1, j, i) = line[static_cast<size_t>(i) * 3 + 1];
          img.at(2, j, i) = line[static_cast<size_t>(i) * 3 + 2];
        }
      }
      head.datamax_org = 255;
      head.bitpix = 8; // already converted to separate colour planes
    } else {
      if (head.bitpix != 8 && head.bitpix != 16 && head.bitpix != 32 && head.bitpix != -32 &&
          head.bitpix != -64) {
        res.error = "Error, unsupported BITPIX!";
        return res;
      }

      // One plane is read in a single call and converted afterwards. Reading row
      // by row and branching on BITPIX per pixel, as the original does, prevents
      // the compiler from vectorising the byte swap and costs roughly half of the
      // total run time of a hinted solve.
      const size_t plane_pixels = static_cast<size_t>(head.width) * head.height;
      const size_t plane_bytes = plane_pixels * bytes_per_pixel;
      std::vector<uint8_t> raw(plane_bytes);

      // Converts the pixels [lo, hi) of `raw` into plane k. Reports the local
      // maximum and whether a NaN was seen. Pure, so chunks may run in any order.
      auto convert_range = [&](const uint8_t *base, float *dst, size_t lo, size_t hi, float &local_max,
                               bool &saw_nan) {
        float m = 0;
        bool nan = false;
        switch (head.bitpix) {
          case 8:
            for (size_t i = lo; i < hi; i++) {
              const double col = base[i] * bscale + bzero;
              dst[i] = static_cast<float>(col);
              if (col > m) m = static_cast<float>(col);
            }
            break;
          case 16:
            for (size_t i = lo; i < hi; i++) {
              const double col = static_cast<int16_t>(be16(base + i * 2)) * bscale + bzero;
              dst[i] = static_cast<float>(col);
              if (col > m) m = static_cast<float>(col);
            }
            break;
          case 32:
            for (size_t i = lo; i < hi; i++) {
              const double col = static_cast<int32_t>(be32(base + i * 4)) * bscale + bzero;
              dst[i] = static_cast<float>(col);
              if (col > m) m = static_cast<float>(col);
            }
            break;
          case -32:
            for (size_t i = lo; i < hi; i++) {
              uint32_t u = be32(base + i * 4);
              float fv;
              std::memcpy(&fv, &u, 4);
              const double col = fv * bscale + bzero;
              if (std::isnan(col)) {
                nan = true;
                continue; // the sequential pass fixes these up
              }
              dst[i] = static_cast<float>(col);
              if (col > m) m = static_cast<float>(col);
            }
            break;
          default: // -64
            for (size_t i = lo; i < hi; i++) {
              uint64_t u = be64(base + i * 8);
              double dv;
              std::memcpy(&dv, &u, 8);
              const double col = dv * bscale + bzero;
              if (std::isnan(col)) {
                nan = true;
                continue;
              }
              dst[i] = static_cast<float>(col);
              if (col > m) m = static_cast<float>(col);
            }
            break;
        }
        local_max = m;
        saw_nan = nan;
      };

      for (int k = 0; k < head.naxis3; k++) {
        // all colours
        f.read(reinterpret_cast<char *>(raw.data()), static_cast<std::streamsize>(plane_bytes));
        const size_t got = static_cast<size_t>(f.gcount());
        const size_t pixels = std::min(plane_pixels, got / bytes_per_pixel);
        float *dst = img.row(k, 0);

        const unsigned chunks = range_chunks(pixels);
        std::vector<float> maxima(std::max(1u, chunks), 0.0f);
        std::vector<char> nans(std::max(1u, chunks), 0);
        const float max_before_plane = measured_max;

        parallel_ranges(0, pixels, [&](size_t lo, size_t hi, unsigned t) {
          bool nan = false;
          convert_range(raw.data(), dst, lo, hi, maxima[t], nan);
          nans[t] = nan ? 1 : 0;
        });

        // max is order independent, so the reduction is exact.
        for (float m: maxima) measured_max = std::max(measured_max, m);

        // A NaN takes the value of the running maximum at that pixel, which does
        // depend on the order. That is rare (very high floating point values in
        // PS1 images), so redo the plane sequentially when it happens.
        if (std::find(nans.begin(), nans.end(), 1) != nans.end()) {
          measured_max = max_before_plane;
          for (size_t i = 0; i < pixels; i++) {
            double col;
            if (head.bitpix == -32) {
              uint32_t u = be32(raw.data() + i * 4);
              float fv;
              std::memcpy(&fv, &u, 4);
              col = fv * bscale + bzero;
            } else {
              uint64_t u = be64(raw.data() + i * 8);
              double dv;
              std::memcpy(&dv, &u, 8);
              col = dv * bscale + bzero;
            }
            if (std::isnan(col)) col = measured_max;
            dst[i] = static_cast<float>(col);
            if (col > measured_max) measured_max = static_cast<float>(col);
          }
        }
      }

      if (head.bitpix <= -32 || head.bitpix == 32) {
        double scalefactor = 1;
        // Rescale a 0..1 range float (GIMP, Astro Pixel Processor, PixInsight) to
        // 0..65535, or scale down values far above 65535.
        if (measured_max > 0 && (measured_max <= 1.0 * 1.5 || measured_max > 65535 * 1.5))
          scalefactor = 65535 / measured_max;

        if (scalefactor != 1) {
          for (int k = 0; k < head.naxis3; k++)
            for (int j = 0; j < head.height; j++) {
              float *row = img.row(k, j);
              for (int i = 0; i < head.width; i++) row[i] = static_cast<float>(row[i] * scalefactor);
            }
          head.datamax_org = 65535;
        } else {
          head.datamax_org = measured_max;
        }
      } else if (head.bitpix == 8) {
        head.datamax_org = 255;
      } else {
        head.datamax_org = measured_max;
      }
    }

    res.ok = true;
    return res;
  }

  bool write_fits_header_file(const std::string &filename, const std::vector<std::string> &cards) {
    std::ofstream f(filename, std::ios::binary);
    if (!f.is_open()) return false;

    size_t written = 0;
    for (const std::string &c: cards) {
      std::string line = c;
      line.resize(80, ' ');
      f.write(line.data(), 80);
      written += 80;
    }
    // Pad to a multiple of 2880 bytes with blanks.
    const std::string blanks(80, ' ');
    while (written % 2880 != 0) {
      f.write(blanks.data(), 80);
      written += 80;
    }
    return static_cast<bool>(f);
  }
} // namespace astap
