# Writes the reference files image_io_tests reads from tests/data, which are not
# in the repository, so a fresh checkout can run the whole suite.
#
#   pip install astropy pillow numpy
#   python tools/make_test_images.py tests/data
#
# The pattern is the one image_io_tests.cpp documents: an 8x6 frame holding
# value(x, y) = 1000*y + 100*x + 7, unique per pixel, so a flipped or off by one
# load cannot pass by accident. Row 0 is the first row in the file.
#
# The point of the .fz files is to check this port's Rice and GZIP tile decoders
# against a different implementation, so they are written by astropy, which
# compresses through CFITSIO itself. rice_dither_ref.fits is astropy's own
# decode of rice_dither.fz, which is what makes that check a comparison against
# CFITSIO's quantisation arithmetic, random dither table included.
#
# Either half can be skipped: without astropy the FITS files are left alone, and
# without Pillow the PNG, TIFF and JPEG ones are. A check whose file is missing
# reports that rather than passing quietly.

import sys

import numpy as np

W, H = 8, 6

PATTERN = np.array([[1000 * y + 100 * x + 7 for x in range(W)] for y in range(H)])


def write_fits(out):
    from astropy.io import fits

    def comp(data, name, ctype, tile, **kw):
        hdu = fits.CompImageHDU(data=data, compression_type=ctype, tile_shape=tile, **kw)
        fits.HDUList([fits.PrimaryHDU(), hdu]).writeto(f"{out}/{name}", overwrite=True)

    i16 = PATTERN.astype(np.int16)
    fits.PrimaryHDU(data=i16).writeto(f"{out}/plain.fits", overwrite=True)
    comp(i16, "rice_rows.fz", "RICE_1", (1, W))
    comp(i16, "rice_tiles.fz", "RICE_1", (3, 4))
    # astap_cli refuses GZIP tiles; this port reads them when it has zlib.
    comp(i16, "gzip1.fz", "GZIP_1", (1, W))
    comp(i16, "gzip2.fz", "GZIP_2", (1, W))

    # Float, which the test reads back at a tenth of the pattern. Tiles Rice
    # cannot shrink end up as GZIP of the original floats, the second path
    # through the same file.
    f32 = (PATTERN * 0.1).astype(np.float32)
    comp(f32, "rice_float.fz", "RICE_1", (1, W), quantize_level=16)

    # Subtractive dithering, with astropy's own decode alongside as reference.
    comp(f32, "rice_dither.fz", "RICE_1", (1, W), quantize_level=16, quantize_method=2)
    decoded = np.asarray(fits.getdata(f"{out}/rice_dither.fz"), dtype=np.float32)
    fits.PrimaryHDU(data=decoded).writeto(f"{out}/rice_dither_ref.fits", overwrite=True)
    print(f"wrote 8 FITS reference files to {out}")


def write_images(out):
    from PIL import Image

    grey = PATTERN.astype(np.uint16)
    Image.fromarray(grey).save(f"{out}/grey16.png")
    Image.fromarray(grey).save(f"{out}/grey16_lzw.tif", compression="tiff_lzw")
    Image.fromarray(grey).save(f"{out}/grey16_deflate.tif", compression="tiff_deflate")

    # 8 bit colour: r = 10*x, g = 20*y, b = 100.
    rgb = np.zeros((H, W, 3), dtype=np.uint8)
    for y in range(H):
        for x in range(W):
            rgb[y, x] = (10 * x, 20 * y, 100)
    Image.fromarray(rgb, mode="RGB").save(f"{out}/rgb8.png")
    Image.fromarray(rgb, mode="RGB").save(f"{out}/rgb8.tif")

    # Lossy, so the test only checks the geometry and which half is bright.
    half = np.zeros((H, W), dtype=np.uint8)
    half[:, W // 2:] = 255
    Image.fromarray(half, mode="L").save(f"{out}/half_grey.jpg", quality=95)

    print(f"wrote 6 reference images to {out}")


def main(out):
    out = out.rstrip("/\\")
    for step, what in ((write_fits, "astropy"), (write_images, "Pillow")):
        try:
            step(out)
        except ImportError:
            print(f"skipped the {what} files: {what} is not installed")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "tests/data")
