// Netpbm and Portable Float Map reader, ported from load_PPM_PGM_PFM() in
// unit_command_line_general.pas. P5 (grey) and P6 (colour) binary Netpbm with
// maxval 255 or 65535, and PF / Pf float maps as Photoshop and the raw
// converters write them.
//
// Like the original, the rows are stored in file order rather than flipped.
// That is right for a PFM, which is stored bottom row first, and leaves a
// PGM/PPM mirrored the same way astap_cli leaves it, so both agree on the
// solution they report for the same file.

#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

#include "image_loaders.h"

namespace astap {
  namespace imageio {
    namespace {
      // Next header token, skipping whitespace and # comments.
      bool next_token(std::ifstream &f, std::string &token) {
        token.clear();
        int c;
        while ((c = f.get()) != EOF) {
          if (c == '#') {
            while ((c = f.get()) != EOF && c != '\n') {
            }
            continue;
          }
          if (std::isspace(c)) {
            if (!token.empty()) return true; // the whitespace byte is consumed
            continue;
          }
          token += static_cast<char>(c);
        }
        return !token.empty();
      }

      inline float be_float(const uint8_t *p) {
        uint32_t u = (static_cast<uint32_t>(p[0]) << 24) | (static_cast<uint32_t>(p[1]) << 16) |
                     (static_cast<uint32_t>(p[2]) << 8) | p[3];
        float v;
        std::memcpy(&v, &u, 4);
        return v;
      }

      inline float le_float(const uint8_t *p) {
        uint32_t u = (static_cast<uint32_t>(p[3]) << 24) | (static_cast<uint32_t>(p[2]) << 16) |
                     (static_cast<uint32_t>(p[1]) << 8) | p[0];
        float v;
        std::memcpy(&v, &u, 4);
        return v;
      }
    } // namespace

    ImageLoadResult load_pnm(const std::string &filename, Header &head, ImageArray &img) {
      ImageLoadResult res;

      std::ifstream f(filename, std::ios::binary);
      if (!f.is_open()) {
        res.error = "Error, accessing the file!";
        return res;
      }

      head = Header();
      img = ImageArray();

      std::string magic;
      if (!next_token(f, magic) ||
          (magic != "P5" && magic != "P6" && magic != "PF" && magic != "Pf")) {
        res.error = "Error, unknown format! Only binary P5, P6, PF and Pf are read.";
        return res;
      }
      const bool pfm = magic == "PF" || magic == "Pf";
      const bool colour = magic == "P6" || magic == "PF";

      std::string w_str, h_str, range_str;
      if (!next_token(f, w_str) || !next_token(f, h_str) || !next_token(f, range_str)) {
        res.error = "Error, truncated header!";
        return res;
      }

      int width = 0, height = 0;
      double range = 0;
      try {
        width = std::stoi(w_str);
        height = std::stoi(h_str);
        range = std::stod(range_str);
      } catch (...) {
        res.error = "Error, unreadable header values!";
        return res;
      }
      if (width <= 0 || height <= 0) {
        res.error = "Error, invalid image dimensions!";
        return res;
      }

      if (pfm) {
        if (range == 0) {
          res.error = "Error, invalid PFM scale factor!";
          return res;
        }
        head.bitpix = -32;
      } else if (range == 65535) {
        head.bitpix = 16;
      } else if (range == 255) {
        head.bitpix = 8;
      } else {
        res.error = "Error, only 255 and 65535 are supported as maximum value!";
        return res;
      }

      const int naxis3 = colour ? 3 : 1;
      const int bytes_per_sample = pfm ? 4 : (head.bitpix == 16 ? 2 : 1);
      const size_t row_bytes =
          static_cast<size_t>(width) * naxis3 * bytes_per_sample;

      head.width = width;
      head.height = height;
      img.resize(naxis3, height, width);

      // A PFM with a negative scale holds little endian floats, a positive one
      // big endian floats; both are scaled to the 0..65535 range the solver uses.
      const bool little_endian_floats = range < 0;
      const double float_scale = pfm ? 65535 / std::fabs(range) : 0;

      std::vector<uint8_t> row(row_bytes);
      for (int y = 0; y < height; y++) {
        f.read(reinterpret_cast<char *>(row.data()), static_cast<std::streamsize>(row_bytes));
        if (static_cast<size_t>(f.gcount()) != row_bytes) {
          res.error = "Error, unexpected end of file in the image data!";
          return res;
        }
        for (int k = 0; k < naxis3; k++) {
          float *dst = img.row(k, y);
          const uint8_t *src = row.data() + static_cast<size_t>(k) * bytes_per_sample;
          const size_t stride = static_cast<size_t>(naxis3) * bytes_per_sample;
          for (int x = 0; x < width; x++, src += stride) {
            if (pfm) {
              const float v = little_endian_floats ? le_float(src) : be_float(src);
              dst[x] = static_cast<float>(v * float_scale);
            } else if (bytes_per_sample == 2) {
              dst[x] = static_cast<float>((static_cast<uint16_t>(src[0]) << 8) | src[1]);
            } else {
              dst[x] = src[0];
            }
          }
        }
      }

      synthesise_header(head, head.bitpix, naxis3, head.bitpix == 8 ? 255 : 65535);
      res.ok = true;
      return res;
    }
  } // namespace imageio
} // namespace astap
