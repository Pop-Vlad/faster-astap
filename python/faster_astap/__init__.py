"""Fast astrometric plate solving, with the index loaded once and kept.

Two solvers, the same lifecycle for both::

    import faster_astap as fa

    solver = fa.IndexSolver(database_path=r"D:\\astap\\d80")
    solver.load()                      # seconds once; solves are milliseconds
    sol = solver.solve(frame, fov=1.2)
    print(sol.ra_deg, sol.dec_deg, sol.scale_arcsec)

`IndexSolver` looks an image's quads up in a whole-sky table, so it needs no
idea where the telescope was pointing. `SpiralSolver` is the faithful port of
ASTAP's own search: slower and it wants a starting position, but it reproduces
`astap_cli`'s answer digit for digit. They take the same shape of call and
return the same `Solution`, so swapping one for the other is a one line change.

Pixel values are expected on the scale the detector produced - the same
0..65535 range a FITS frame carries. An array normalised to 0..1 is rejected
rather than quietly failing to find any stars in it.

The star database is not part of this package: it is gigabytes, and distributed
separately from www.hnsky.org. `find_databases()` reports where one was found.
"""

from __future__ import annotations

from .results import Solution
from .solvers import IndexSolver, SpiralSolver, find_databases, supported_image_extensions

__all__ = [
    "IndexSolver",
    "SpiralSolver",
    "Solution",
    "find_databases",
    "supported_image_extensions",
    "NoDatabaseError",
    "SolveError",
    "__version__",
]

from ._core import __version__  # noqa: E402
from .solvers import NoDatabaseError, SolveError  # noqa: E402
