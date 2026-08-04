// Checks the image loaders behind load_image(): the pixel values, the number of
// colour planes and the row order every format ends up with.
//
// The reference files in tests/data hold the same 8x6 pattern,
// value(x, y) = 1000*y + 100*x + 7, which is unique per pixel, so a flipped,
// transposed or off by one load cannot pass by accident. The Netpbm and BMP
// files are written by the test itself; the compressed ones were produced by
// astropy (CFITSIO for the Rice tiles) and Pillow, so the decoders are checked
// against other implementations rather than against themselves.

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

#include "astap/image/image_io.h"

using namespace astap;

namespace {
  int failures = 0;

  const int kW = 8, kH = 6;

  // The pattern as stored in the file, row 0 first.
  double expected(int x, int y) { return 1000.0 * y + 100.0 * x + 7; }

  std::string data_dir() { return std::string(ASTAP_TEST_DATA_DIR) + "/"; }
  std::string tmp_dir() { return std::string(ASTAP_TEST_TMP_DIR) + "/"; }

  void check(bool cond, const std::string &what) {
    if (!cond) {
      std::printf("FAIL: %s\n", what.c_str());
      failures++;
    } else {
      std::printf("ok  : %s\n", what.c_str());
    }
  }

  // Loads `file` and compares plane `plane` with the pattern. `flipped` is true
  // for the formats whose first row is the top of the image, which the loaders
  // store last to follow the FITS convention. `scale` and `tol` allow for the
  // conversions a format applies to its samples.
  void check_image(const std::string &file, int want_planes, bool flipped, double scale,
                   double tol, const std::string &what) {
    Header head;
    ImageArray img;
    const ImageLoadResult r = load_image(file, head, img);
    if (!r.ok) {
      std::printf("FAIL: %s (%s)\n", what.c_str(), r.error.c_str());
      failures++;
      return;
    }
    if (head.width != kW || head.height != kH || img.colours() != want_planes) {
      std::printf("FAIL: %s (got %dx%d, %d plane(s))\n", what.c_str(), head.width, head.height,
                  img.colours());
      failures++;
      return;
    }

    double worst = 0;
    for (int y = 0; y < kH; y++)
      for (int x = 0; x < kW; x++) {
        const int row = flipped ? kH - 1 - y : y;
        worst = std::max(worst, std::fabs(img.at(0, row, x) - expected(x, y) * scale));
      }
    if (worst > tol) {
      std::printf("FAIL: %s (largest pixel error %.6g, tolerance %g)\n", what.c_str(), worst, tol);
      failures++;
      return;
    }
    std::printf("ok  : %-46s %dx%d, %d plane(s), max error %.4g\n", what.c_str(), head.width,
                head.height, img.colours(), worst);
  }

  void write_file(const std::string &name, const std::vector<uint8_t> &bytes) {
    std::ofstream f(name, std::ios::binary);
    f.write(reinterpret_cast<const char *>(bytes.data()),
            static_cast<std::streamsize>(bytes.size()));
  }

  void append(std::vector<uint8_t> &v, const std::string &s) {
    v.insert(v.end(), s.begin(), s.end());
  }

  // --- the formats the test writes itself -----------------------------------

  // P5 / P6, 8 or 16 bits. Netpbm rows are stored in file order by the original
  // and by this port, so no flip is expected.
  std::string write_pnm(bool colour, int maxval) {
    std::vector<uint8_t> v;
    append(v, (colour ? std::string("P6\n") : std::string("P5\n")) + "# written by image_io_tests\n" +
                  std::to_string(kW) + " " + std::to_string(kH) + "\n" + std::to_string(maxval) +
                  "\n");
    for (int y = 0; y < kH; y++)
      for (int x = 0; x < kW; x++)
        for (int k = 0; k < (colour ? 3 : 1); k++) {
          const int value = maxval == 255 ? static_cast<int>(expected(x, y)) % 256
                                          : static_cast<int>(expected(x, y));
          if (maxval == 65535) v.push_back(static_cast<uint8_t>(value >> 8));
          v.push_back(static_cast<uint8_t>(value & 0xFF));
        }
    const std::string name =
        tmp_dir() + (colour ? "pattern_p6_" : "pattern_p5_") + std::to_string(maxval) +
        (colour ? ".ppm" : ".pgm");
    write_file(name, v);
    return name;
  }

  // Portable float map, little endian, scale -1 so the values are 0..1 and come
  // back multiplied by 65535.
  std::string write_pfm() {
    std::vector<uint8_t> v;
    append(v, "Pf\n" + std::to_string(kW) + " " + std::to_string(kH) + "\n-1.0\n");
    for (int y = 0; y < kH; y++)
      for (int x = 0; x < kW; x++) {
        const float value = static_cast<float>(expected(x, y) / 65535.0);
        uint8_t b[4];
        std::memcpy(b, &value, 4);
        v.insert(v.end(), b, b + 4);
      }
    const std::string name = tmp_dir() + "pattern.pfm";
    write_file(name, v);
    return name;
  }

  void put32(std::vector<uint8_t> &v, uint32_t x) {
    v.push_back(static_cast<uint8_t>(x));
    v.push_back(static_cast<uint8_t>(x >> 8));
    v.push_back(static_cast<uint8_t>(x >> 16));
    v.push_back(static_cast<uint8_t>(x >> 24));
  }

  void put16(std::vector<uint8_t> &v, uint16_t x) {
    v.push_back(static_cast<uint8_t>(x));
    v.push_back(static_cast<uint8_t>(x >> 8));
  }

  // 24 bit bottom up bitmap holding the pattern modulo 256, so the loader has
  // to flip nothing: the first row in the file is the bottom row of the image.
  std::string write_bmp() {
    const size_t stride = (static_cast<size_t>(kW) * 3 + 3) / 4 * 4;
    std::vector<uint8_t> pixels;
    for (int y = kH - 1; y >= 0; y--) { // bottom row first
      const size_t start = pixels.size();
      for (int x = 0; x < kW; x++) {
        const uint8_t value = static_cast<uint8_t>(static_cast<int>(expected(x, y)) % 256);
        pixels.push_back(value); // blue
        pixels.push_back(value); // green
        pixels.push_back(value); // red
      }
      pixels.resize(start + stride, 0);
    }

    std::vector<uint8_t> v;
    append(v, "BM");
    put32(v, static_cast<uint32_t>(14 + 40 + pixels.size()));
    put32(v, 0);
    put32(v, 14 + 40);
    put32(v, 40);
    put32(v, kW);
    put32(v, kH);
    put16(v, 1);
    put16(v, 24);
    put32(v, 0); // BI_RGB
    put32(v, static_cast<uint32_t>(pixels.size()));
    put32(v, 2835);
    put32(v, 2835);
    put32(v, 0);
    put32(v, 0);
    v.insert(v.end(), pixels.begin(), pixels.end());

    const std::string name = tmp_dir() + "pattern.bmp";
    write_file(name, v);
    return name;
  }

  // --- checks ---------------------------------------------------------------

  void test_netpbm_and_bmp() {
    check_image(write_pnm(false, 65535), 1, false, 1, 0, "P5 16 bit grey");
    check_image(write_pnm(true, 65535), 3, false, 1, 0, "P6 16 bit colour");
    check_image(write_pfm(), 1, false, 1, 1e-3, "Pf 32 bit float");

    // The 8 bit variants hold the pattern modulo 256.
    Header head;
    ImageArray img;
    const std::string p5 = write_pnm(false, 255);
    check(load_image(p5, head, img).ok && img.at(0, 2, 3) == static_cast<int>(expected(3, 2)) % 256,
          "P5 8 bit grey");

    const std::string bmp = write_bmp();
    const ImageLoadResult r = load_image(bmp, head, img);
    const int v = static_cast<int>(expected(3, 2)) % 256;
    check(r.ok && img.colours() == 3 && head.width == kW && head.height == kH &&
              img.at(0, kH - 1 - 2, 3) == v * 257 && img.at(2, kH - 1 - 2, 3) == v * 257,
          "BMP 24 bit, bottom up");
  }

  // Compares two files that must decode to the same pixels.
  void check_same(const std::string &a, const std::string &b, double tol,
                  const std::string &what) {
    Header ha, hb;
    ImageArray ia, ib;
    const ImageLoadResult ra = load_image(a, ha, ia);
    const ImageLoadResult rb = load_image(b, hb, ib);
    if (!ra.ok || !rb.ok) {
      std::printf("FAIL: %s (%s)\n", what.c_str(), (ra.ok ? rb.error : ra.error).c_str());
      failures++;
      return;
    }
    if (ia.width() != ib.width() || ia.height() != ib.height() ||
        ia.colours() != ib.colours()) {
      std::printf("FAIL: %s (different geometry)\n", what.c_str());
      failures++;
      return;
    }
    double worst = 0;
    for (int c = 0; c < ia.colours(); c++)
      for (int y = 0; y < ia.height(); y++)
        for (int x = 0; x < ia.width(); x++)
          worst = std::max(worst, std::fabs(static_cast<double>(ia.at(c, y, x)) - ib.at(c, y, x)));
    if (worst > tol) {
      std::printf("FAIL: %s (largest difference %.6g, tolerance %g)\n", what.c_str(), worst, tol);
      failures++;
      return;
    }
    std::printf("ok  : %-46s %dx%d, largest difference %.4g\n", what.c_str(), ia.width(),
                ia.height(), worst);
  }

  void test_fits() {
    check_image(data_dir() + "plain.fits", 1, false, 1, 0, "uncompressed FITS");
    check_image(data_dir() + "rice_rows.fz", 1, false, 1, 0, "Rice .fz, one tile per row");
    check_image(data_dir() + "rice_tiles.fz", 1, false, 1, 0, "Rice .fz, 3x4 tiles");

    // GZIP tiles. astap_cli refuses these, this port reads them when it has
    // zlib, and the same code covers the tiles a Rice file could not compress.
    for (const char *file: {"gzip1.fz", "gzip2.fz"}) {
      Header head;
      ImageArray img;
      const ImageLoadResult r = load_image(data_dir() + file, head, img);
      if (!r.ok && r.error.find("GZIP") != std::string::npos) {
        std::printf("skip: %-46s built without zlib\n", file);
        continue;
      }
      check_image(data_dir() + file, 1, false, 1, 0, std::string(file) + ", one tile per row");
    }

    // Tiles too small for Rice to win are stored as GZIP of the *original*
    // floats, which is a second path through the same file. Without zlib those
    // tiles are skipped and reported, exactly as astap_cli does.
    {
      Header head;
      ImageArray img;
      const std::string file = data_dir() + "rice_float.fz";
      const ImageLoadResult r = load_image(file, head, img);
      if (r.ok && r.warning.find("zlib") != std::string::npos)
        std::printf("skip: %-46s built without zlib\n", "Rice .fz, GZIP fallback tiles");
      else
        check_image(file, 1, false, 0.1, 0.02, "Rice .fz, float with GZIP fallback tiles");
    }

    // Quantised with subtractive dithering. The reference is the same image
    // decoded by astropy, so this compares against CFITSIO's own arithmetic
    // including the random table.
    //
    // This file has GZIP fallback tiles in it too, so like the one above it can
    // only be compared in full when the build has zlib. Without it those tiles
    // are skipped and the comparison would fail on their contents rather than
    // on the dithering it is meant to be testing.
    {
      Header head;
      ImageArray img;
      const std::string file = data_dir() + "rice_dither.fz";
      const ImageLoadResult r = load_image(file, head, img);
      if (r.ok && r.warning.find("zlib") != std::string::npos)
        std::printf("skip: %-46s built without zlib\n", "Rice .fz, dithered quantisation");
      else
        check_same(file, data_dir() + "rice_dither_ref.fits", 0.001,
                   "Rice .fz, dithered quantisation vs reference decode");
    }
  }

  void test_optional_formats() {
    struct Case {
      const char *file;
      int planes;
      double scale;
      double tol;
      const char *what;
    };
    const Case cases[] = {
        {"grey16.png", 1, 1, 0, "PNG 16 bit grey"},
        {"grey16_lzw.tif", 1, 1, 0, "TIFF 16 bit grey, LZW"},
        {"grey16_deflate.tif", 1, 1, 0, "TIFF 16 bit grey, deflate"},
    };
    for (const Case &c: cases) {
      Header head;
      ImageArray img;
      const std::string file = data_dir() + c.file;
      const ImageLoadResult r = load_image(file, head, img);
      if (!r.ok && r.error.find("this build reads no") != std::string::npos) {
        std::printf("skip: %-46s not built with the decoder\n", c.what);
        continue;
      }
      check_image(file, c.planes, true, c.scale, c.tol, c.what);
    }

    // 8 bit colour: r = 10*x, g = 20*y, b = 100, each scaled to 16 bits.
    for (const char *file: {"rgb8.png", "rgb8.tif"}) {
      Header head;
      ImageArray img;
      const ImageLoadResult r = load_image(data_dir() + file, head, img);
      if (!r.ok && r.error.find("this build reads no") != std::string::npos) {
        std::printf("skip: %-46s not built with the decoder\n", file);
        continue;
      }
      bool ok = r.ok && img.colours() == 3 && head.width == kW && head.height == kH;
      if (ok)
        for (int y = 0; y < kH && ok; y++)
          for (int x = 0; x < kW && ok; x++) {
            const int row = kH - 1 - y;
            ok = img.at(0, row, x) == 10 * x * 257 && img.at(1, row, x) == 20 * y * 257 &&
                 img.at(2, row, x) == 100 * 257;
          }
      check(ok, std::string(file) + " 8 bit colour");
    }

    // Lossy: only the geometry and which half of the frame is bright.
    {
      Header head;
      ImageArray img;
      const ImageLoadResult r = load_image(data_dir() + "half_grey.jpg", head, img);
      if (!r.ok && r.error.find("this build reads no") != std::string::npos) {
        std::printf("skip: %-46s not built with the decoder\n", "JPEG grey");
        return;
      }
      const bool ok = r.ok && head.width == kW && head.height == kH && img.colours() == 1 &&
                      img.at(0, 3, 1) < 5000 && img.at(0, 3, kW - 1) > 55000;
      check(ok, "JPEG 8 bit grey");
    }
  }

  void test_rejections() {
    Header head;
    ImageArray img;
    const ImageLoadResult r = load_image(tmp_dir() + "nothing.xyz", head, img);
    check(!r.ok && r.error.find("unsupported file type") != std::string::npos,
          "unknown extension is refused with the list of what is read");

    const ImageLoadResult m = load_image(data_dir() + "does_not_exist.fits", head, img);
    check(!m.ok, "missing file is refused");
  }
} // namespace

int main() {
  std::printf("-- Netpbm and BMP --\n");
  test_netpbm_and_bmp();
  std::printf("\n-- FITS, plain and Rice compressed --\n");
  test_fits();
  std::printf("\n-- optional decoders --\n");
  test_optional_formats();
  std::printf("\n-- refusals --\n");
  test_rejections();

  if (failures) {
    std::printf("\n%d check(s) FAILED\n", failures);
    return 1;
  }
  std::printf("\nall image loader checks passed\n");
  return 0;
}
