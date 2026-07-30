#include "astap/fits.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <vector>

#include "astap/astro_math.h"

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
  } // namespace

  FitsLoadResult load_fits(const std::string &filename, Header &head, ImageArray &img) {
    FitsLoadResult res;

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

    // Read the primary header, 2880 byte blocks of 36 cards.
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
        }
      }
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

    // Read the image data.
    const int bytes_per_pixel = rgb24 ? 3 : std::abs(head.bitpix) / 8;
    if (bytes_per_pixel <= 0) {
      res.error = "Error, unsupported BITPIX!";
      return res;
    }

    img.resize(head.naxis3, head.height, head.width);

    std::vector<uint8_t> line(static_cast<size_t>(head.width) * bytes_per_pixel);
    float measured_max = 0;

    auto read_line = [&]() -> bool {
      f.read(reinterpret_cast<char *>(line.data()), static_cast<std::streamsize>(line.size()));
      return f.gcount() == static_cast<std::streamsize>(line.size());
    };

    if (rgb24) {
      for (int j = 0; j < head.height; j++) {
        if (!read_line()) break;
        for (int i = 0; i < head.width; i++) {
          img.at(0, j, i) = line[static_cast<size_t>(i) * 3 + 0];
          img.at(1, j, i) = line[static_cast<size_t>(i) * 3 + 1];
          img.at(2, j, i) = line[static_cast<size_t>(i) * 3 + 2];
        }
      }
      head.datamax_org = 255;
      head.bitpix = 8; // already converted to separate colour planes
    } else {
      for (int k = 0; k < head.naxis3; k++) {
        // all colours
        for (int j = 0; j < head.height; j++) {
          if (!read_line()) break;
          for (int i = 0; i < head.width; i++) {
            const uint8_t *p = line.data() + static_cast<size_t>(i) * bytes_per_pixel;
            double col;
            switch (head.bitpix) {
              case 8:
                col = p[0] * bscale + bzero;
                break;
              case 16:
                col = static_cast<int16_t>(be16(p)) * bscale + bzero;
                break;
              case 32:
                col = static_cast<int32_t>(be32(p)) * bscale + bzero;
                break;
              case -32: {
                uint32_t u = be32(p);
                float fv;
                std::memcpy(&fv, &u, 4);
                col = fv * bscale + bzero;
                // Not a number, can happen in PS1 images with very high values.
                if (std::isnan(col)) col = measured_max;
                break;
              }
              case -64: {
                uint64_t u = be64(p);
                double dv;
                std::memcpy(&dv, &u, 8);
                col = dv * bscale + bzero;
                if (std::isnan(col)) col = measured_max;
                break;
              }
              default:
                res.error = "Error, unsupported BITPIX!";
                return res;
            }
            img.at(k, j, i) = static_cast<float>(col);
            // Find the maximum value, needed for images with a 0..1 scale.
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
