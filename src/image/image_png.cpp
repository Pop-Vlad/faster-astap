// PNG reader. astap_cli hands the file to the fcl-image PNG reader, which
// reports whether the image is grey scale and returns 16 bit channel values
// (an 8 bit sample becomes v * 257), after which the rows are stored bottom up
// in the FITS convention. libpng is used here for the same result.
//
// Only compiled with a decoder when the build found libpng.

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "image_loaders.h"

#ifdef ASTAP_WITH_PNG
#include <png.h>
#endif

namespace astap {
  namespace imageio {
#ifdef ASTAP_WITH_PNG
    ImageLoadResult load_png(const std::string &filename, Header &head, ImageArray &img) {
      ImageLoadResult res;

      FILE *fp = std::fopen(filename.c_str(), "rb");
      if (!fp) {
        res.error = "Error, accessing the file!";
        return res;
      }

      head = Header();
      img = ImageArray();

      png_structp png = png_create_read_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
      png_infop info = png ? png_create_info_struct(png) : nullptr;
      // Everything the error path touches has to survive a longjmp, so no
      // containers and no std::string here.
      png_bytep volatile pixels = nullptr;
      png_bytepp volatile rows = nullptr;
      const char *volatile failure = "Error, reading the PNG file!";

      if (!png || !info || setjmp(png_jmpbuf(png))) {
        if (pixels) std::free(pixels);
        if (rows) std::free(rows);
        if (png) png_destroy_read_struct(&png, info ? &info : nullptr, nullptr);
        std::fclose(fp);
        res.error = failure;
        return res;
      }

      png_init_io(png, fp);
      png_read_info(png, info);

      const png_uint_32 width = png_get_image_width(png, info);
      const png_uint_32 height = png_get_image_height(png, info);
      const int bit_depth = png_get_bit_depth(png, info);
      const int colour_type = png_get_color_type(png, info);

      if (colour_type == PNG_COLOR_TYPE_PALETTE) png_set_palette_to_rgb(png);
      if (colour_type == PNG_COLOR_TYPE_GRAY && bit_depth < 8) png_set_expand_gray_1_2_4_to_8(png);
      if (png_get_valid(png, info, PNG_INFO_tRNS)) png_set_tRNS_to_alpha(png);
      png_set_strip_alpha(png);
      png_set_interlace_handling(png);
      png_read_update_info(png, info);

      const int channels = png_get_channels(png, info);
      const int depth = png_get_bit_depth(png, info);
      const size_t row_bytes = png_get_rowbytes(png, info);

      pixels = static_cast<png_bytep>(std::malloc(row_bytes * height));
      rows = static_cast<png_bytepp>(std::malloc(sizeof(png_bytep) * height));
      if (!pixels || !rows) {
        failure = "Error, not enough memory for the image!";
        png_longjmp(png, 1);
      }
      for (png_uint_32 y = 0; y < height; y++) rows[y] = pixels + y * row_bytes;
      png_read_image(png, rows);

      const int naxis3 = channels >= 3 ? 3 : 1;
      head.width = static_cast<int>(width);
      head.height = static_cast<int>(height);
      img.resize(naxis3, head.height, head.width);

      for (png_uint_32 y = 0; y < height; y++) {
        const png_bytep src = rows[y];
        const int dst_y = head.height - 1 - static_cast<int>(y); // FITS is bottom up
        for (int k = 0; k < naxis3; k++) {
          float *dst = img.row(k, dst_y);
          for (png_uint_32 x = 0; x < width; x++) {
            const png_bytep p = src + (static_cast<size_t>(x) * channels + k) * (depth / 8);
            // 16 bit samples arrive most significant byte first, whatever the
            // byte order of this machine.
            dst[x] = depth == 16 ? static_cast<float>((p[0] << 8) | p[1])
                                 : static_cast<float>(p[0] * 257);
          }
        }
      }

      std::free(pixels);
      std::free(rows);
      png_destroy_read_struct(&png, &info, nullptr);
      std::fclose(fp);

      synthesise_header(head, 16, naxis3, 65535);
      res.ok = true;
      return res;
    }
#else
    ImageLoadResult load_png(const std::string &, Header &, ImageArray &) {
      return missing_support("PNG", "libpng");
    }
#endif
  } // namespace imageio
} // namespace astap
