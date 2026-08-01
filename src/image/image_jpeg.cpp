// JPEG reader. astap_cli uses the fcl-image JPEG reader, which reports whether
// the image is grey scale and returns 16 bit channel values (v * 257 for the 8
// bit samples a JPEG always has), stored bottom up in the FITS convention.
// libjpeg is used here for the same result. Baseline and progressive files are
// both handled by libjpeg itself.
//
// Only compiled with a decoder when the build found libjpeg.

#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

#include "image_loaders.h"

#ifdef ASTAP_WITH_JPEG
#include <csetjmp>

#include <jpeglib.h>
#endif

namespace astap {
  namespace imageio {
#ifdef ASTAP_WITH_JPEG
    namespace {
      struct JpegError {
        jpeg_error_mgr pub;
        jmp_buf escape;
        char message[JMSG_LENGTH_MAX];
      };

      void jpeg_error_exit(j_common_ptr cinfo) {
        JpegError *err = reinterpret_cast<JpegError *>(cinfo->err);
        (*cinfo->err->format_message)(cinfo, err->message);
        std::longjmp(err->escape, 1);
      }

      void jpeg_no_message(j_common_ptr) {
      }

      // The decoding proper lives in its own frame: a longjmp out of libjpeg
      // skips the destructors of everything below it, so no container with an
      // allocation may sit in the frame that owns the setjmp.
      void decode(jpeg_decompress_struct &cinfo, Header &head, ImageArray &img, int naxis3) {
        const int width = static_cast<int>(cinfo.output_width);
        const int height = static_cast<int>(cinfo.output_height);
        const int channels = cinfo.output_components;

        head.width = width;
        head.height = height;
        img.resize(naxis3, height, width);

        std::vector<uint8_t> scanline(static_cast<size_t>(width) * channels);
        while (cinfo.output_scanline < cinfo.output_height) {
          uint8_t *line = scanline.data();
          const int y = static_cast<int>(cinfo.output_scanline);
          jpeg_read_scanlines(&cinfo, &line, 1);
          const int dst_y = height - 1 - y; // FITS is bottom up
          for (int k = 0; k < naxis3; k++) {
            float *dst = img.row(k, dst_y);
            for (int x = 0; x < width; x++)
              dst[x] = scanline[static_cast<size_t>(x) * channels + k] * 257;
          }
        }
      }
    } // namespace

    ImageLoadResult load_jpeg(const std::string &filename, Header &head, ImageArray &img) {
      ImageLoadResult res;

      FILE *fp = std::fopen(filename.c_str(), "rb");
      if (!fp) {
        res.error = "Error, accessing the file!";
        return res;
      }

      head = Header();
      img = ImageArray();

      jpeg_decompress_struct cinfo;
      JpegError err;
      cinfo.err = jpeg_std_error(&err.pub);
      err.pub.error_exit = jpeg_error_exit;
      err.pub.output_message = jpeg_no_message; // warnings go nowhere, like fcl-image
      err.message[0] = '\0';

      int naxis3 = 1;
      if (setjmp(err.escape)) {
        jpeg_destroy_decompress(&cinfo);
        std::fclose(fp);
        res.error = std::string("Error, reading the JPEG file: ") +
                    (err.message[0] ? err.message : "decoding failed");
        return res;
      }

      jpeg_create_decompress(&cinfo);
      jpeg_stdio_src(&cinfo, fp);
      jpeg_read_header(&cinfo, TRUE);

      // Grey scale stays grey scale, everything else is decoded to RGB.
      cinfo.out_color_space = cinfo.jpeg_color_space == JCS_GRAYSCALE ? JCS_GRAYSCALE : JCS_RGB;
      jpeg_start_decompress(&cinfo);

      naxis3 = cinfo.output_components >= 3 ? 3 : 1;
      decode(cinfo, head, img, naxis3);

      jpeg_finish_decompress(&cinfo);
      jpeg_destroy_decompress(&cinfo);
      std::fclose(fp);

      synthesise_header(head, 16, naxis3, 65535);
      res.ok = true;
      return res;
    }
#else
    ImageLoadResult load_jpeg(const std::string &, Header &, ImageArray &) {
      return missing_support("JPEG", "libjpeg");
    }
#endif
  } // namespace imageio
} // namespace astap
