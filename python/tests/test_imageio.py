"""The image loader, which decides what a file's pixels mean.

The C++ tests already check the C++ decoders against reference files. What is
worth testing here is the half of the loader that is Python: the conventions
Pillow and rawpy output has to be put into before the solver sees it - the 257
and 65535 scale factors, the bottom up row order, the colour reduction - and the
dispatch that decides which decoder owns an extension in the first place. Those
are what make a PNG solve to the same coordinates as the FITS it came from.

The fixtures in tests/data are the same image in several formats, so a Pillow
decoded file can be compared against the C++ decoded one directly rather than
against numbers written down here.
"""

from __future__ import annotations

import numpy as np
import pytest
import struct
import zlib
from pathlib import Path

from faster_astap import _core, _imageio

DATA = Path(__file__).resolve().parents[2] / "tests" / "data"

pytest.importorskip("PIL", reason="Pillow is required to read PNG, JPEG and TIFF")


def _data(name: str) -> str:
    path = DATA / name
    if not path.is_file():
        pytest.skip(f"{name} is missing from tests/data")
    return str(path)


def _write_png16(path, array):
    """A 16 bit PNG, which Pillow can read but cannot write.

    Only needed to produce the one case Pillow decodes lossily, so that the
    warning about it can be tested at all.
    """
    height, width, colours = array.shape

    def chunk(kind, payload):
        return (struct.pack(">I", len(payload)) + kind + payload
                + struct.pack(">I", zlib.crc32(kind + payload) & 0xFFFFFFFF))

    raw = b"".join(b"\x00" + array[y].astype(">u2").tobytes() for y in range(height))
    Path(path).write_bytes(
        b"\x89PNG\r\n\x1a\n"
        + chunk(b"IHDR", struct.pack(">IIBBBBB", width, height, 16,
                                     {1: 0, 3: 2}[colours], 0, 0, 0))
        + chunk(b"IDAT", zlib.compress(raw))
        + chunk(b"IEND", b""))


# --- who owns which extension -------------------------------------------------
def test_extension_ownership_does_not_overlap():
    """One decoder per extension, so a file reads the same on every machine."""
    core = set(_imageio.CORE_EXTENSIONS)
    assert not core & set(_imageio.PILLOW_EXTENSIONS)
    assert not core & set(_imageio.RAW_EXTENSIONS)
    assert not set(_imageio.PILLOW_EXTENSIONS) & set(_imageio.RAW_EXTENSIONS)


def test_unsupported_extension_names_what_is_readable(tmp_path):
    nonsense = tmp_path / "frame.xyz"
    nonsense.write_bytes(b"")
    with pytest.raises(ValueError, match="unsupported file type"):
        _imageio.load(nonsense)


def test_supported_extensions_lists_the_installed_decoders():
    listed = _imageio.supported_extensions()
    assert ".png" in listed and ".jpg" in listed  # Pillow is a dependency
    if _core.has_image_io:
        assert ".fits" in listed


# --- the conventions Pillow output has to be put into -------------------------
def test_eight_bit_is_scaled_by_257_and_stored_bottom_up(tmp_path):
    from PIL import Image

    top_down = np.array([[0, 1], [2, 255]], dtype=np.uint8)
    path = tmp_path / "grey.png"
    Image.fromarray(top_down, mode="L").save(path)

    data, meta = _imageio.load(path)
    assert data.shape == (1, 2, 2)
    assert data.dtype == np.float32
    # Row order reversed, and every sample multiplied by 257 as image_png.cpp does.
    np.testing.assert_array_equal(data[0], top_down[::-1].astype(np.float32) * 257)
    assert not meta["warning"]


def test_sixteen_bit_grey_is_taken_as_it_is(tmp_path):
    values = np.array([[7, 1000], [40000, 65535]], dtype=np.uint16)
    path = tmp_path / "grey16.png"
    _write_png16(path, values[:, :, np.newaxis])

    data, meta = _imageio.load(path)
    np.testing.assert_array_equal(data[0], values[::-1].astype(np.float32))
    assert not meta["warning"]


def test_alpha_is_dropped_and_colour_reduced_to_three(tmp_path):
    from PIL import Image

    rgba = np.zeros((2, 2, 4), dtype=np.uint8)
    rgba[..., :3] = 40
    rgba[..., 3] = 128
    path = tmp_path / "rgba.png"
    Image.fromarray(rgba, mode="RGBA").save(path)

    data, _meta = _imageio.load(path)
    assert data.shape == (3, 2, 2)
    assert np.all(data == 40 * 257)


def test_sixteen_bit_colour_truncation_is_reported(tmp_path):
    """Pillow drops the low byte of 16 bit colour and says nothing.

    The C++ readers keep all 16, so this is a real difference and the caller is
    told about it rather than left to find out from a solve that finds no stars.
    """
    values = np.array([[[1000, 2000, 3000], [40000, 50000, 60000]],
                       [[65535, 0, 12345], [7, 8, 9]]], dtype=np.uint16)
    path = tmp_path / "rgb16.png"
    _write_png16(path, values)

    data, meta = _imageio.load(path)
    assert data.shape == (3, 2, 2)
    assert "8 bits per channel" in meta["warning"]
    # Truncated, not scaled: 1000 >> 8 is 3, and 3 * 257 is what comes back.
    assert data[0, 1, 0] == (1000 >> 8) * 257


def test_eight_bit_colour_is_not_falsely_flagged(tmp_path):
    from PIL import Image

    path = tmp_path / "rgb8.png"
    Image.fromarray(np.full((2, 2, 3), 17, dtype=np.uint8), mode="RGB").save(path)
    _data_, meta = _imageio.load(path)
    assert not meta["warning"]


# --- the float guard read_tiff() applies in ASTAP -----------------------------
def _write_float_tiff(path, array):
    """A minimal single strip float TIFF, since this is about sample width."""
    height, width = array.shape
    pixels = array.astype("<f4").tobytes()
    ntags, header = 11, 8
    strip = header + 2 + ntags * 12 + 4

    def tag(code, kind, count, value):
        return struct.pack("<HHII", code, kind, count, value)

    tags = [tag(256, 3, 1, width), tag(257, 3, 1, height), tag(258, 3, 1, 32),
            tag(259, 3, 1, 1), tag(262, 3, 1, 1), tag(273, 4, 1, strip),
            tag(277, 3, 1, 1), tag(278, 3, 1, height),
            tag(279, 4, 1, len(pixels)), tag(284, 3, 1, 1), tag(339, 3, 1, 3)]
    Path(path).write_bytes(b"II" + struct.pack("<HI", 42, header)
                           + struct.pack("<H", ntags) + b"".join(tags)
                           + struct.pack("<I", 0) + pixels)


def test_floats_within_range_are_only_scaled(tmp_path):
    path = tmp_path / "small.tif"
    _write_float_tiff(path, np.linspace(0, 1, 4, dtype=np.float32).reshape(2, 2))
    data, _meta = _imageio.load(path)
    assert data.max() == pytest.approx(65535.0, abs=0.5)


def test_floats_past_the_range_are_rescaled_not_clipped(tmp_path):
    """read_tiff() in unit_tiff_unthreaded.pas normalises a float image whose
    values run past 65535 * 1.5, rather than handing the solver numbers it has
    no scale for. A float TIFF written in 0..1 is not the only kind there is."""
    path = tmp_path / "big.tif"
    _write_float_tiff(path, np.linspace(0, 4, 4, dtype=np.float32).reshape(2, 2))
    data, _meta = _imageio.load(path)
    assert data.max() == pytest.approx(65535.0, abs=0.5)
    # Linear, so the ratios inside the image survive.
    positive = data[data > 0]
    assert float(positive.max() / positive.min()) == pytest.approx(3.0, rel=1e-4)


def test_floats_below_the_guard_are_left_alone(tmp_path):
    path = tmp_path / "edge.tif"
    _write_float_tiff(path, np.full((2, 2), 1.4999, dtype=np.float32))
    data, _meta = _imageio.load(path)
    assert data.max() > 65535.0  # under 65535 * 1.5, so untouched


# --- the header every loader synthesises --------------------------------------
def test_metadata_matches_what_the_core_loader_returns(tmp_path):
    from PIL import Image

    path = tmp_path / "grey.png"
    Image.fromarray(np.zeros((4, 4), dtype=np.uint8), mode="L").save(path)
    _pixels, meta = _imageio.load(path)

    for key in ("ra_mount", "dec_mount", "focallen", "warning", "cards",
                "cdelt2", "ra0", "dec0"):
        assert key in meta
    # 99999 is the solver's own "not specified"; a PNG carries no pointing.
    assert meta["ra_mount"] >= 999 and meta["dec_mount"] >= 999
    assert all(len(card) == 80 for card in meta["cards"])
    assert meta["cards"][0].startswith("SIMPLE  =")
    assert meta["cards"][-1].startswith("END")


def test_colour_image_header_carries_naxis3():
    _pixels, meta = _imageio.load(_data("rgb8.png"))
    assert any(card.startswith("NAXIS3") for card in meta["cards"])


# --- the point of all of it: the formats have to agree ------------------------
@pytest.mark.skipif(not _core.has_image_io, reason="this build reads no image files")
def test_pillow_png_agrees_with_the_core_fits_reader():
    """The Pillow reader has to put pixels where the C++ readers put them.

    plain.fits and grey16.png hold the same array written twice by
    tools/make_test_images.py, which makes them vertical mirrors of each other on
    disk: FITS numbers its rows from the bottom and PNG from the top. Both
    readers store bottom up, so the mirror survives into the arrays, and that is
    the point being checked - the values must agree exactly, and the row order
    must differ by exactly one flip. Get either half wrong and a PNG exported
    from a FITS solves to a mirrored sky.
    """
    from_fits, _ = _imageio.load(_data("plain.fits"))
    from_png, _ = _imageio.load(_data("grey16.png"))
    assert from_fits.shape == from_png.shape
    np.testing.assert_array_equal(from_fits, from_png[:, ::-1, :])


def test_png_and_tiff_agree_with_each_other():
    from_png, _ = _imageio.load(_data("rgb8.png"))
    from_tiff, _ = _imageio.load(_data("rgb8.tif"))
    np.testing.assert_array_equal(from_png, from_tiff)


@pytest.mark.parametrize("name", ["grey16_lzw.tif", "grey16_deflate.tif"])
def test_compressed_tiffs_read_the_same_as_the_uncompressed_png(name):
    reference, _ = _imageio.load(_data("grey16.png"))
    compressed, _ = _imageio.load(_data(name))
    np.testing.assert_array_equal(reference, compressed)


# --- raw camera files ---------------------------------------------------------
def test_raw_extensions_are_routed_to_rawpy(tmp_path):
    """No raw file ships with the package, so what is checked here is that a raw
    extension reaches rawpy at all rather than falling through to Pillow."""
    pytest.importorskip("rawpy", reason="rawpy is required to read raw camera files")
    truncated = tmp_path / "frame.cr2"
    truncated.write_bytes(b"not a raw file")
    with pytest.raises(ValueError, match="accessing the file"):
        _imageio.load(truncated)
