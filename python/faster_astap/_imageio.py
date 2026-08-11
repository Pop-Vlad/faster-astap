"""Turning an image file into the array the solvers take.

Which decoder owns which extension is fixed here rather than discovered at
import time. The C++ extension reads the formats that need nothing beyond the
compiler - FITS, Rice compressed .fz, Netpbm and BMP - and Pillow and rawpy read
the rest. A wheel is built with the optional C++ decoders switched off (the
comment above ASTAP_PNG in pyproject.toml says why), so PNG, JPEG, TIFF and the
raw camera formats always come through Python, and always the same way on every
platform.

Nothing here falls back from one decoder to another. Each extension has exactly
one owner, whatever a given build happens to have available, because a loader
that tries C++ first and Python second is a loader whose pixel values depend on
how the machine was built - which is the thing the fixed ownership avoids.

The pixel conventions are the ones the C++ loaders use, since both feed the same
solver:

  - 8 bit samples are scaled by 257 and 32 bit floats by 65535, 16 bit samples
    are taken as they are (src/image/image_png.cpp, src/image/image_tiff.cpp).
  - rows are stored bottom up, the FITS convention (include/astap/image/image_io.h).
  - a raw camera file is the undemosaiced sensor frame, margins included and
    rows in file order rather than flipped (src/image/image_raw.cpp).

Pillow decodes a 16 bit colour PNG or TIFF to 8 bits per channel and says
nothing about it, which the C++ readers do not do - they keep all 16. The
warning in the returned metadata is this module saying it instead. It costs
accuracy in no measurable way: on the corpus, a truncated image that solves at
all agrees with the full precision one to under a tenth of a pixel. What it can
cost is the solve itself, when the background noise sits below ~256 counts and
truncation flattens the faint stars into a single level.
"""

from __future__ import annotations

import importlib.util
import os
from typing import Any, Dict, Tuple

import numpy as np

from . import _core

# The formats the C++ extension reads with nothing linked behind them, and so
# the ones a wheel can rely on having. Kept in step with load_image() in
# src/image/image_io.cpp.
CORE_EXTENSIONS = (".fit", ".fits", ".fts", ".fz", ".new",
                   ".ppm", ".pgm", ".pfm", ".pbm", ".bmp")

PILLOW_EXTENSIONS = (".png", ".jpg", ".jpeg", ".tif", ".tiff")

# The same list is_raw_extension() gives LibRaw in src/image/image_io.cpp.
RAW_EXTENSIONS = (".cr2", ".cr3", ".crw", ".nef", ".nrw", ".arw", ".srf", ".sr2",
                  ".orf", ".rw2", ".raf", ".dng", ".pef", ".raw", ".3fr", ".fff",
                  ".iiq", ".mos", ".mef", ".mrw", ".erf", ".kdc", ".dcr", ".srw",
                  ".x3f", ".rwl", ".dcs", ".cap", ".bay")


class ImageReadError(ValueError):
    """Raised when a file cannot be read.

    A ValueError, because that is what _core.load_image_file has always raised
    and callers catch.
    """


# --- the FITS header every loader synthesises ---------------------------------
# The cards synthesise_header() writes in src/image/image_io.cpp, card for card,
# so a solution built from a Pillow decoded file carries the same header as one
# built from a FITS file.

def _card(key: str, value_and_comment: str) -> str:
    return (key.ljust(9) + value_and_comment).ljust(80)


def _integer_card(key: str, value: float, comment: str) -> str:
    # Python's round() breaks ties to even, the same as FPC's Round() that
    # pround() reproduces on the C++ side.
    return _card(key, f"{round(value)}".rjust(20) + comment)


def _synthesise_header(width: int, height: int, bitpix: int, naxis3: int,
                       datamax: float) -> list:
    naxis = 2 if naxis3 == 1 else 3
    cards = [
        _card("SIMPLE  =", "                    T / FITS header                                    "),
        _integer_card("BITPIX  =", bitpix, " / Bits per entry                                 "),
        _integer_card("NAXIS   =", naxis, " / Number of dimensions                           "),
        _integer_card("NAXIS1  =", width, " / length of x axis                               "),
        _integer_card("NAXIS2  =", height, " / length of y axis                               "),
    ]
    if naxis3 != 1:
        cards.append(
            _integer_card("NAXIS3  =", naxis3, " / length of z axis (mostly colors)               "))
    cards += [
        _card("EQUINOX =", "               2000.0 / Equinox of coordinates                         "),
        _integer_card("DATAMIN =", 0, " / Minimum data value                             "),
        _integer_card("DATAMAX =", datamax, " / Maximum data value                             "),
        _card("BZERO   =", "                  0.0 / Physical_value = BZERO + BSCALE * array_value  "),
        _card("BSCALE  =", "                  1.0 / Physical_value = BZERO + BSCALE * array_value  "),
        _card("COMMENT 1", "  Written by ASTAP, Astrometric STAcking Program. www.hnsky.org        "),
        _card("END", ""),
    ]
    return cards


def _meta(cards: list, warning: str = "") -> Dict[str, Any]:
    """The metadata dict load_image_file() returns, for a file whose format
    carries no pointing information at all. Every one of these has a header, so
    the keys are always present and the solver's own defaults are what they say:
    99999 for a coordinate means "not specified"."""
    return {
        "ra_mount": 99999.0,
        "dec_mount": 99999.0,
        "focallen": 0.0,
        "warning": warning,
        "cards": cards,
        "cdelt2": 0.0,
        "ra0": 0.0,
        "dec0": 0.0,
    }


# --- Pillow: PNG, JPEG, TIFF --------------------------------------------------

# Modes Pillow hands back that map straight onto a sample width, with the
# multiplier that puts them on the solver's 0..65535 scale.
_SCALE_BY_MODE = {
    "L": 257.0,      # 8 bit grey
    "RGB": 257.0,    # 8 bit colour, and where a 16 bit colour file lands too
    "I;16": 1.0,     # 16 bit grey, the one 16 bit mode Pillow keeps intact
    "I;16B": 1.0,
    "I;16L": 1.0,
    "I;16N": 1.0,
    "F": 65535.0,    # 32 bit float, as in image_tiff.cpp
    "I": 1.0,        # 32 bit integer TIFF; already on a large scale of its own
}

# Alpha is dropped and a palette expanded, the same two things image_png.cpp
# asks libpng for with png_set_strip_alpha and png_set_palette_to_rgb.
_CONVERT_TO = {
    "1": "L", "LA": "L", "La": "L",
    "P": "RGB", "PA": "RGB", "RGBA": "RGB", "RGBa": "RGB", "RGBX": "RGB",
    "CMYK": "RGB", "YCbCr": "RGB", "HSV": "RGB", "LAB": "RGB",
}

# The two modes whose values are not bounded by their own sample width, and so
# the ones the rescale below can apply to.
_UNBOUNDED_MODES = ("F", "I")


def _load_pillow(path: str) -> Tuple[np.ndarray, Dict[str, Any]]:
    try:
        from PIL import Image, UnidentifiedImageError
    except ImportError:  # pragma: no cover - only on a partial install
        raise ImageReadError(
            "Error, reading this format needs Pillow, which is not installed. "
            "Install it with: pip install pillow") from None

    try:
        with Image.open(path) as opened:
            # The first frame, for a multi page TIFF; the solver takes one image.
            opened.seek(0)
            source_mode = opened.mode
            truncated = _is_truncated_colour(opened, path)
            image = opened.convert(_CONVERT_TO[source_mode]) \
                if source_mode in _CONVERT_TO else opened
            mode = image.mode
            scale = _SCALE_BY_MODE.get(mode)
            if scale is None:
                raise ImageReadError(
                    f"Error, unsupported {source_mode} image; this is not a grey scale "
                    "or colour image the solver can use.")
            pixels = np.asarray(image)
    except ImageReadError:
        raise
    except (UnidentifiedImageError, OSError, ValueError) as exc:
        raise ImageReadError(f"Error, reading the file! {exc}") from None

    if pixels.ndim == 2:
        pixels = pixels[:, :, np.newaxis]
    # Three colours or one, as the C++ loaders reduce to; anything between is
    # already gone by the conversion above.
    colours = 3 if pixels.shape[2] >= 3 else 1
    data = pixels[:, :, :colours].astype(np.float32) * scale

    # Rows bottom up, the FITS convention, then (colours, height, width) as the
    # solver indexes it. Both are views; the copy is taken once, at the end.
    data = np.ascontiguousarray(data[::-1].transpose(2, 0, 1))

    warning = ""
    if mode in _UNBOUNDED_MODES:
        # read_tiff() in ASTAP's unit_tiff_unthreaded.pas normalises a float
        # image whose values run past the 16 bit range back into it, rather than
        # handing the solver numbers it has no scale for. The port's libtiff
        # reader left this out; it is here because the guard costs one pass and
        # a float TIFF written in 0..1 is not the only kind there is.
        measured_max = float(data.max()) if data.size else 0.0
        if measured_max > 65535.0 * 1.5:
            data *= np.float32(65535.0 / measured_max)
    if truncated:
        warning = ("PNG/TIFF read through Pillow, reduced to 8 bits per channel; "
                   "a 16 bit colour image loses its low byte.")

    bitpix = -32 if mode == "F" else 16
    cards = _synthesise_header(data.shape[2], data.shape[1], bitpix, colours, 65535)
    return data, _meta(cards, warning)


def _is_truncated_colour(image, path: str) -> bool:
    """Whether Pillow is about to drop the low byte of a 16 bit colour image.

    Pillow decodes 16 bit colour to 8 bit and reports it as plain RGB, with
    nothing on the image to tell it from a file that was 8 bit to begin with, so
    the sample width has to come from the format's own metadata. 16 bit grey
    survives as I;16 and is not affected.
    """
    if image.mode not in ("RGB", "RGBA"):
        return False
    if image.format == "TIFF":
        bits = image.tag_v2.get(258) if hasattr(image, "tag_v2") else None
        if isinstance(bits, (tuple, list)):
            return any(int(b) > 8 for b in bits)
        return bits is not None and int(bits) > 8
    if image.format == "PNG":
        # Straight from the IHDR rather than from Pillow: the bit depth is the
        # 25th byte of every PNG there is, and the attribute that would carry it
        # is private and cleared once the image has been loaded.
        try:
            with open(path, "rb") as handle:
                header = handle.read(26)
        except OSError:
            return False
        return len(header) >= 26 and header[12:16] == b"IHDR" and header[24] == 16
    return False


# --- rawpy: the raw camera formats --------------------------------------------

def _load_rawpy(path: str) -> Tuple[np.ndarray, Dict[str, Any]]:
    try:
        import rawpy
    except ImportError:  # pragma: no cover - only on a partial install
        raise ImageReadError(
            "Error, reading raw camera files needs rawpy, which is not installed. "
            "Install it with: pip install rawpy") from None

    try:
        with rawpy.imread(path) as raw:
            sensor = raw.raw_image
            if sensor is None or sensor.ndim != 2:
                raise ImageReadError(
                    "Error, this raw file holds no single channel sensor image "
                    "(Foveon and already demosaiced files are not supported).")
            # raw_image includes the masked margins, which is what LibRaw's
            # unprocessed_raw writes and what image_raw.cpp reads; the slice only
            # trims the row padding rawpy may carry past raw_width.
            sizes = raw.sizes
            sensor = sensor[:sizes.raw_height, :sizes.raw_width]
            if sensor.size == 0:
                raise ImageReadError("Error, invalid raw image dimensions!")
            # Copied, not viewed: the buffer belongs to LibRaw and is freed here.
            data = np.array(sensor, dtype=np.float32)[np.newaxis, :, :]
    except ImageReadError:
        raise
    except Exception as exc:  # rawpy raises a family of LibRawError subclasses
        raise ImageReadError(f"Error, accessing the file! {exc}") from None

    # Row order as it comes off the sensor, not flipped: image_raw.cpp stores it
    # the way the PGM that ASTAP's own converters write is stored.
    measured_max = float(data.max()) if data.size else 0.0
    cards = _synthesise_header(data.shape[2], data.shape[1], 16, 1,
                               65535 if measured_max > 255 else 255)
    return data, _meta(cards)


# --- what the rest of the package calls ---------------------------------------

def _extension(path: str) -> str:
    return os.path.splitext(path)[1].lower()


def _installed(module: str) -> bool:
    """Whether a decoder is importable, without paying to import it."""
    try:
        return importlib.util.find_spec(module) is not None
    except (ImportError, ValueError):  # a broken or partial installation
        return False


def supported_extensions() -> str:
    """The extensions this installation can read, space separated."""
    parts = []
    if _core.has_image_io:
        parts.extend(CORE_EXTENSIONS)
    if _installed("PIL"):
        parts.extend(PILLOW_EXTENSIONS)
    if _installed("rawpy"):
        parts.append("and raw camera files (.cr2 .cr3 .nef .arw .dng ...)")
    return " ".join(parts)


def load(path) -> Tuple[np.ndarray, Dict[str, Any]]:
    """Reads an image file.

    Returns the pixel values as a (colours, height, width) float32 array and the
    metadata dict `_core.load_image_file` returns, whichever decoder was used.
    Raises ImageReadError, a ValueError, when the file cannot be read.
    """
    path = os.fspath(path)
    extension = _extension(path)

    if extension in CORE_EXTENSIONS:
        if not _core.has_image_io:
            raise ImageReadError(
                f"Error, this build reads no {extension} files; it was built with "
                "ASTAP_PYTHON_IMAGE_IO off. Pass an array instead.")
        return _core.load_image_file(path)
    if extension in PILLOW_EXTENSIONS:
        return _load_pillow(path)
    if extension in RAW_EXTENSIONS:
        return _load_rawpy(path)

    raise ImageReadError(
        f"Error, unsupported file type '{extension or path}'. This installation reads: "
        + supported_extensions())
