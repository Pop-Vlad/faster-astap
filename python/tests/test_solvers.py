"""Both solvers, through the Python API.

The C++ tests already check that solving an array agrees with solving the file
it came from. What is worth testing here is the layer above that: the argument
handling, the conversions, the Solution, and the lifecycle — the things a caller
touches and the C++ tests never see.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import faster_astap as fa


# --- things that need nothing at all -----------------------------------------
def test_find_databases_returns_paths():
    dirs = fa.find_databases()
    assert isinstance(dirs, list)
    assert all(isinstance(d, str) for d in dirs)


def test_unloaded_solver_refuses_to_solve():
    solver = fa.IndexSolver()
    assert not solver.ready
    with pytest.raises(fa.SolveError, match="not loaded"):
        solver.solve(np.zeros((64, 64), dtype=np.uint16))


def test_missing_database_raises_no_database_error(tmp_path):
    with pytest.raises(fa.NoDatabaseError):
        fa.IndexSolver(database_path=tmp_path).load()


# --- image conversion, which is this layer's own work -------------------------
def test_normalised_image_is_refused(index_solver):
    """0..1 input would find no stars at all; that deserves a real message."""
    with pytest.raises(ValueError, match="normalised"):
        index_solver.solve(np.linspace(0, 1, 64 * 64).reshape(64, 64).astype(np.float32))


def test_wrong_dimensions_are_refused(index_solver):
    with pytest.raises(ValueError, match="2 or 3 dimensional"):
        index_solver.solve(np.zeros((4, 4, 4, 4), dtype=np.uint16))


def test_tiny_image_is_refused(index_solver):
    with pytest.raises(ValueError, match="too small"):
        index_solver.solve(np.zeros((4, 4), dtype=np.uint16))


@pytest.mark.parametrize("dtype", [np.uint8, np.uint16, np.int32, np.float32, np.float64])
def test_every_sane_dtype_reaches_the_solver(index_solver, pixels, dtype):
    """An integer camera frame should not need converting by the caller."""
    mono = pixels[0]
    scaled = mono if dtype != np.uint8 else mono / mono.max() * 255
    solution = index_solver.solve(np.ascontiguousarray(scaled, dtype=dtype))
    assert solution.solved


def test_colour_axis_either_way_round(index_solver, pixels):
    """(c, h, w) and (h, w, c) must give the same answer."""
    chw = np.repeat(pixels, 3, axis=0)
    hwc = np.ascontiguousarray(np.transpose(chw, (1, 2, 0)))
    a, b = index_solver.solve(chw), index_solver.solve(hwc)
    assert a.solved and b.solved
    assert a.ra == pytest.approx(b.ra, abs=1e-9)
    assert a.dec == pytest.approx(b.dec, abs=1e-9)


# --- the Solution ------------------------------------------------------------
def test_solution_units_and_truthiness(index_solver, pixels):
    solution = index_solver.solve(pixels)
    assert solution.solved and bool(solution)

    assert solution.ra_deg == pytest.approx(math.degrees(solution.ra))
    assert solution.ra_hours == pytest.approx(solution.ra_deg / 15)
    assert solution.dec_deg == pytest.approx(math.degrees(solution.dec))
    assert solution.scale_arcsec == pytest.approx(abs(solution.cdelt2) * 3600)
    assert solution.scale_arcsec > 0

    width_deg, height_deg = solution.fov_deg
    assert 0 < width_deg < 180 and 0 < height_deg < 180
    assert solution.nr_inliers >= 3
    assert solution.stars > 0
    assert "Solution found" in " ".join(solution.messages)


def test_unsolved_solution_is_falsey(index_solver):
    """Noise has stars in it, but no sky matches them."""
    rng = np.random.default_rng(7)
    noise = rng.integers(0, 5000, size=(512, 512), dtype=np.uint16)
    solution = index_solver.solve(noise)
    assert not solution.solved
    assert not solution
    assert "unsolved" in repr(solution)


# --- lifecycle ---------------------------------------------------------------
def test_solve_many_reuses_one_load(index_solver, pixels):
    """The point of the whole design: load once, solve repeatedly."""
    first = index_solver.solve(pixels)
    second = index_solver.solve(pixels)
    assert first.solved and second.solved
    assert first.ra == pytest.approx(second.ra, abs=1e-12)


def test_context_manager_loads_and_closes(database, pixels):
    with fa.IndexSolver(database_path=database) as solver:
        assert solver.ready
        assert solver.solve(pixels).solved
    assert not solver.ready


def test_progress_callback_is_called(index_solver, pixels):
    lines = []
    index_solver.solve(pixels, progress=lines.append)
    assert lines and all(isinstance(line, str) for line in lines)


def test_a_raising_progress_callback_does_not_break_the_solve(index_solver, pixels):
    def explode(_line):
        raise RuntimeError("callbacks should not be able to take a solve down")

    assert index_solver.solve(pixels, progress=explode).solved


# --- the port ----------------------------------------------------------------
def test_spiral_solver_agrees_with_the_index(database, image_path, pixels):
    """Two independent methods, one sky: they have to land in the same place."""
    index = fa.IndexSolver(database_path=database).load()
    reference = index.solve(image_path)
    assert reference.solved

    port = fa.SpiralSolver(database_path=database).load()
    got = port.solve(image_path)
    assert got.solved

    separation = math.degrees(math.hypot(
        (got.ra - reference.ra) * math.cos(reference.dec), got.dec - reference.dec))
    assert separation * 3600 < 5.0  # arcseconds
    assert got.scale_arcsec == pytest.approx(reference.scale_arcsec, rel=1e-3)
