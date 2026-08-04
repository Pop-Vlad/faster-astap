"""Shared fixtures.

Everything here needs a star database and a real image, neither of which ships
with the package: a database is gigabytes and comes from www.hnsky.org, and the
test images are fetched with tools/fetch_skyview_corpus.py. Without them the
tests skip rather than fail, so a checkout with neither still reports green for
the things it can actually check.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

import faster_astap as fa

REPO = Path(__file__).resolve().parents[2]


def _first_database() -> str | None:
    """A directory holding a usable star database, or None."""
    override = os.environ.get("FASTER_ASTAP_DATABASE")
    candidates = [override] if override else fa.find_databases()
    for directory in candidates:
        if not directory:
            continue
        path = Path(directory)
        if path.is_dir() and any(path.glob("*.1476")) or any(path.glob("*.290")):
            return str(path)
    return None


@pytest.fixture(scope="session")
def database() -> str:
    found = _first_database()
    if not found:
        pytest.skip("no star database found; set FASTER_ASTAP_DATABASE to one")
    return found


@pytest.fixture(scope="session")
def image_path() -> str:
    """A corpus image around a degree across, which both solvers handle well.

    The corpus deliberately contains sizes nothing can solve — 0.05 degree
    cutouts are down to a handful of stars — so picking by name matters.
    """
    for directory in (REPO / "corpus", REPO / "build" / "mini-corpus"):
        if not directory.is_dir():
            continue
        images = sorted(
            directory.glob("*.fits"),
            key=lambda p: abs(_fov_from_name(p.name) - 1.0),
        )
        if images:
            return str(images[0])
    pytest.skip("no corpus images; fetch them with tools/fetch_skyview_corpus.py")


def _fov_from_name(name: str) -> float:
    """Field size out of <field>_<survey>_<fov>deg.fits, or 0 when absent."""
    marker = name.rfind("deg.fits")
    if marker < 0:
        return 0.0
    start = name.rfind("_", 0, marker)
    try:
        return float(name[start + 1 : marker])
    except ValueError:
        return 0.0


@pytest.fixture(scope="session")
def index_solver(database):
    solver = fa.IndexSolver(database_path=database)
    solver.load()
    return solver


@pytest.fixture(scope="session")
def pixels(image_path):
    """The test image as an array, so tests do not each re-read it."""
    from faster_astap import _core

    if not _core.has_image_io:
        pytest.skip("this build reads no image files")
    data, _meta = _core.load_image_file(image_path)
    return data
