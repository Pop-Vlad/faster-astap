#!/usr/bin/env python3
"""Convert a PNG (or any image Pillow can read) into a 16 bit FITS file.

The C++ port reads FITS only, while the ASTAP GUI and astap_cli also accept
PNG/TIFF/JPEG. Use this helper to feed such an image to astap_solve, for
example the test_img.png in this directory:

    python3 tools/png_to_fits.py test_img.png test_img.fits
    ./build/astap_solve -f test_img.fits -fov 0 -r 180 -d /path/to/star_database

Requires Pillow and numpy.
"""
import numpy as np
import sys
from PIL import Image

Image.MAX_IMAGE_PIXELS = None


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        return 1
    src, dst = sys.argv[1], sys.argv[2]

    im = Image.open(src)
    a = np.asarray(im)
    if a.ndim == 3:
        planes = [a[:, :, k] for k in range(min(3, a.shape[2]))]
    else:
        planes = [a]

    h, w = planes[0].shape
    # FITS stores the first row at the bottom, PNG at the top.
    scale = 257 if planes[0].dtype == np.uint8 else 1  # 8 bit -> full 16 bit range
    planes = [np.flipud(p).astype(np.int32) * scale for p in planes]
    naxis3 = len(planes)

    cards = [
        "SIMPLE  =                    T",
        "BITPIX  =                   16",
        "NAXIS   = %20d" % (3 if naxis3 > 1 else 2),
        "NAXIS1  = %20d" % w,
        "NAXIS2  = %20d" % h,
    ]
    if naxis3 > 1:
        cards.append("NAXIS3  = %20d" % naxis3)
    cards += [
        "BZERO   =              32768.0",  # unsigned data stored as signed 16 bit
        "BSCALE  =                  1.0",
        "EQUINOX =               2000.0",
        "END",
    ]

    with open(dst, "wb") as f:
        blob = b"".join(("%-80s" % c).encode("ascii") for c in cards)
        blob += b" " * ((2880 - len(blob) % 2880) % 2880)
        f.write(blob)
        written = 0
        for p in planes:
            raw = (np.clip(p, 0, 65535) - 32768).astype(">i2")
            f.write(raw.tobytes())
            written += raw.nbytes
        f.write(b"\0" * ((2880 - written % 2880) % 2880))

    print("wrote %s, %dx%dx%d" % (dst, w, h, naxis3))
    return 0


if __name__ == "__main__":
    sys.exit(main())
