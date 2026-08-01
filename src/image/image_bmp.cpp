// Windows bitmap reader, the .bmp branch of load_PNGJPEG() in
// unit_command_line_general.pas. The original hands the file to the fcl-image
// BMP reader, which returns 16 bit channel values (an 8 bit sample becomes
// v * 257), and then stores the rows bottom up in the FITS convention. That is
// what this reproduces, so a .bmp reads back the same way it does in astap_cli.
//
// Uncompressed bitmaps of 1, 4, 8, 16, 24 and 32 bits per pixel are read,
// including BI_BITFIELDS. The RLE encoded variants are not.

#include <cstdint>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

#include "image_loaders.h"

namespace astap {
  namespace imageio {
    namespace {
      inline uint16_t le16(const uint8_t *p) {
        return static_cast<uint16_t>(p[0] | (p[1] << 8));
      }

      inline uint32_t le32(const uint8_t *p) {
        return static_cast<uint32_t>(p[0]) | (static_cast<uint32_t>(p[1]) << 8) |
               (static_cast<uint32_t>(p[2]) << 16) | (static_cast<uint32_t>(p[3]) << 24);
      }

      // Number of trailing zero bits, and the number of set bits, of a channel
      // mask, so an arbitrary bitfield can be expanded to the full 16 bit range.
      void mask_shape(uint32_t mask, int &shift, int &bits) {
        shift = 0;
        bits = 0;
        if (mask == 0) return;
        while (((mask >> shift) & 1) == 0) shift++;
        for (uint32_t m = mask >> shift; m & 1; m >>= 1) bits++;
      }

      inline float expand(uint32_t value, int bits) {
        if (bits <= 0) return 0;
        const uint32_t maxv = (bits >= 32) ? 0xFFFFFFFFu : ((1u << bits) - 1);
        return static_cast<float>(static_cast<double>(value) * 65535 / maxv);
      }
    } // namespace

    ImageLoadResult load_bmp(const std::string &filename, Header &head, ImageArray &img) {
      ImageLoadResult res;

      std::ifstream f(filename, std::ios::binary);
      if (!f.is_open()) {
        res.error = "Error, accessing the file!";
        return res;
      }

      head = Header();
      img = ImageArray();

      uint8_t file_header[14];
      f.read(reinterpret_cast<char *>(file_header), 14);
      if (f.gcount() != 14 || file_header[0] != 'B' || file_header[1] != 'M') {
        res.error = "Error, not a BMP file!";
        return res;
      }
      const uint32_t pixel_offset = le32(file_header + 10);

      uint8_t size_field[4];
      f.read(reinterpret_cast<char *>(size_field), 4);
      if (f.gcount() != 4) {
        res.error = "Error, truncated BMP header!";
        return res;
      }
      const uint32_t dib_size = le32(size_field);
      if (dib_size < 12 || dib_size > 1024) {
        res.error = "Error, unsupported BMP header!";
        return res;
      }
      std::vector<uint8_t> dib(dib_size);
      std::memcpy(dib.data(), size_field, 4);
      f.read(reinterpret_cast<char *>(dib.data()) + 4, static_cast<std::streamsize>(dib_size - 4));
      if (static_cast<uint32_t>(f.gcount()) != dib_size - 4) {
        res.error = "Error, truncated BMP header!";
        return res;
      }

      int width = 0, height = 0, bpp = 0;
      uint32_t compression = 0, palette_entries = 0;
      bool core_header = dib_size == 12;
      if (core_header) {
        width = le16(dib.data() + 4);
        height = le16(dib.data() + 6);
        bpp = le16(dib.data() + 10);
      } else {
        width = static_cast<int32_t>(le32(dib.data() + 4));
        height = static_cast<int32_t>(le32(dib.data() + 8));
        bpp = le16(dib.data() + 14);
        compression = le32(dib.data() + 16);
        palette_entries = le32(dib.data() + 32);
      }

      const bool top_down = height < 0;
      if (top_down) height = -height;
      if (width <= 0 || height <= 0) {
        res.error = "Error, invalid BMP dimensions!";
        return res;
      }
      if (compression != 0 && compression != 3) {
        res.error = "Error, compressed BMP files are not supported!";
        return res;
      }
      if (bpp != 1 && bpp != 4 && bpp != 8 && bpp != 16 && bpp != 24 && bpp != 32) {
        res.error = "Error, unsupported BMP colour depth (" + std::to_string(bpp) + " bits)!";
        return res;
      }

      // Channel masks. BI_RGB has fixed layouts; BI_BITFIELDS stores its masks
      // either in the header (V4 and up) or right behind it.
      uint32_t rmask = 0, gmask = 0, bmask = 0;
      if (bpp == 16 || bpp == 32) {
        if (compression == 3) {
          if (dib_size >= 52) {
            rmask = le32(dib.data() + 40);
            gmask = le32(dib.data() + 44);
            bmask = le32(dib.data() + 48);
          } else {
            uint8_t masks[12];
            f.read(reinterpret_cast<char *>(masks), 12);
            if (f.gcount() != 12) {
              res.error = "Error, truncated BMP bit masks!";
              return res;
            }
            rmask = le32(masks);
            gmask = le32(masks + 4);
            bmask = le32(masks + 8);
          }
        } else if (bpp == 16) {
          rmask = 0x7C00; // 5-5-5
          gmask = 0x03E0;
          bmask = 0x001F;
        } else {
          rmask = 0x00FF0000;
          gmask = 0x0000FF00;
          bmask = 0x000000FF;
        }
      }

      // Palette, for the indexed depths.
      std::vector<uint8_t> palette;
      if (bpp <= 8) {
        uint32_t entries = palette_entries;
        if (entries == 0) entries = 1u << bpp;
        const uint32_t entry_size = core_header ? 3 : 4;
        palette.resize(static_cast<size_t>(entries) * entry_size);
        f.read(reinterpret_cast<char *>(palette.data()),
               static_cast<std::streamsize>(palette.size()));
        if (static_cast<size_t>(f.gcount()) != palette.size()) {
          res.error = "Error, truncated BMP palette!";
          return res;
        }
        // Store as 4 byte BGRA entries so the pixel loop needs one layout only.
        if (core_header) {
          std::vector<uint8_t> wide(static_cast<size_t>(entries) * 4, 0);
          for (uint32_t i = 0; i < entries; i++) std::memcpy(&wide[i * 4], &palette[i * 3], 3);
          palette.swap(wide);
        }
      }

      f.seekg(static_cast<std::streamoff>(pixel_offset), std::ios::beg);
      if (!f) {
        res.error = "Error, invalid BMP pixel offset!";
        return res;
      }

      head.width = width;
      head.height = height;
      img.resize(3, height, width);

      int rshift = 0, rbits = 0, gshift = 0, gbits = 0, bshift = 0, bbits = 0;
      mask_shape(rmask, rshift, rbits);
      mask_shape(gmask, gshift, gbits);
      mask_shape(bmask, bshift, bbits);

      const size_t stride = (static_cast<size_t>(width) * bpp + 31) / 32 * 4;
      std::vector<uint8_t> row(stride);

      for (int y = 0; y < height; y++) {
        f.read(reinterpret_cast<char *>(row.data()), static_cast<std::streamsize>(stride));
        if (static_cast<size_t>(f.gcount()) != stride) {
          res.error = "Error, unexpected end of file in the image data!";
          return res;
        }
        // Bitmap rows run bottom up unless the height was negative, and image
        // row 0 is the bottom row, so only a top down file needs flipping.
        const int dst_y = top_down ? height - 1 - y : y;
        float *r = img.row(0, dst_y);
        float *g = img.row(1, dst_y);
        float *b = img.row(2, dst_y);

        for (int x = 0; x < width; x++) {
          uint8_t red = 0, green = 0, blue = 0;
          switch (bpp) {
            case 1:
            case 4:
            case 8: {
              const int per_byte = 8 / bpp;
              const uint8_t byte = row[static_cast<size_t>(x) / per_byte];
              const int shift = (per_byte - 1 - x % per_byte) * bpp;
              const uint32_t index = (byte >> shift) & ((1u << bpp) - 1);
              const size_t entry = static_cast<size_t>(index) * 4;
              if (entry + 2 < palette.size()) {
                blue = palette[entry];
                green = palette[entry + 1];
                red = palette[entry + 2];
              }
              break;
            }
            case 24:
              blue = row[static_cast<size_t>(x) * 3];
              green = row[static_cast<size_t>(x) * 3 + 1];
              red = row[static_cast<size_t>(x) * 3 + 2];
              break;
            default: { // 16 and 32 bit, through the channel masks
              const uint32_t v = bpp == 16 ? le16(&row[static_cast<size_t>(x) * 2])
                                           : le32(&row[static_cast<size_t>(x) * 4]);
              r[x] = expand((v & rmask) >> rshift, rbits);
              g[x] = expand((v & gmask) >> gshift, gbits);
              b[x] = expand((v & bmask) >> bshift, bbits);
              continue;
            }
          }
          // 8 bit samples become 16 bit ones the same way fcl-image does it.
          r[x] = red * 257;
          g[x] = green * 257;
          b[x] = blue * 257;
        }
      }

      synthesise_header(head, 16, 3, 65535);
      res.ok = true;
      return res;
    }
  } // namespace imageio
} // namespace astap
