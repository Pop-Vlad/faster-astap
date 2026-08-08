"""The two solvers, with the setup separated from the solving.

Both classes are built the same way: construct with everything that describes
the *installation* - where the star database is, how deep an index to build -
then `load()` once, then `solve()` per image. `load()` is the expensive call and
the reason these are objects rather than functions.
"""

from __future__ import annotations

import os
from typing import Any, Callable, Optional, Sequence, Union

from . import _core
from .results import Solution, from_outcome

PathLike = Union[str, "os.PathLike[str]"]


class SolveError(RuntimeError):
    """Raised when a solve could not be attempted at all."""


class NoDatabaseError(SolveError):
    """Raised when no star database could be found or read."""


def find_databases() -> list:
    """The directories a star database is looked for in, in order.

    Only says where it looks, not what it found; a solver reports which database
    it actually selected once loaded.
    """
    return list(_core.default_database_directories())


def supported_image_extensions() -> str:
    """Image file extensions this build can read, space separated.

    Empty when the extension was built without image support, in which case
    `solve()` takes arrays only.
    """
    return _core.supported_image_extensions()


class _Solver:
    """Shared lifecycle. Subclasses supply the backend and the per image call."""

    _impl_factory: Callable[[], Any]

    def __init__(self, *, database_path: Optional[PathLike] = None, database: str = "auto"):
        self._impl = self._impl_factory()
        self._database_path = "" if database_path is None else os.fspath(database_path)
        self._database = database
        self._loaded = False

    # --- lifecycle ------------------------------------------------------------
    @property
    def ready(self) -> bool:
        return self._loaded and self._impl.ready

    def load(self, *, log: Optional[Callable[[str], None]] = None) -> "_Solver":
        """Gets everything resident. Returns self, so it can be chained."""
        raise NotImplementedError

    def close(self) -> None:
        """Drops what `load` made resident."""
        self._impl = self._impl_factory()
        self._loaded = False

    def __enter__(self) -> "_Solver":
        if not self.ready:
            self.load()
        return self

    def __exit__(self, *exc) -> None:
        self.close()

    def _require_ready(self) -> None:
        if not self.ready:
            raise SolveError(
                f"{type(self).__name__} is not loaded; call load() first, or use it as a "
                "context manager"
            )

    # --- what is resident -----------------------------------------------------
    @property
    def database(self) -> str:
        """The star database actually selected, e.g. 'd80'."""
        return self._impl.database

    @property
    def database_path(self) -> str:
        return self._impl.database_path

    def _solve_file(self, path: PathLike, **kwargs) -> Solution:
        """Reads a file and solves the array, so both solvers share one path."""
        if not _core.has_image_io:
            raise SolveError(
                "this build reads no image files; pass an array instead "
                "(astropy.io.fits or your own loader)"
            )
        pixels, meta = _core.load_image_file(os.fspath(path))
        kwargs.setdefault("label", os.fspath(path))
        return pixels, meta, kwargs


class IndexSolver(_Solver):
    """Blind plate solving against a pre-built whole-sky quad index.

    Needs no start position and no field size, though a field size makes it
    quicker. `load()` reads the index — seconds the first time it has to build
    one, and effectively nothing afterwards, since the cache is memory mapped.
    """

    _impl_factory = staticmethod(_core.IndexSolver)

    def __init__(
            self,
            *,
            database_path: Optional[PathLike] = None,
            database: str = "auto",
            quad_tolerance: float = 0.007,
            ladder: Optional[Sequence[float]] = None,
            index_cache: Optional[PathLike] = None,
            use_cache: bool = True,
            rebuild: bool = False,
    ):
        super().__init__(database_path=database_path, database=database)
        self._quad_tolerance = quad_tolerance
        self._ladder = list(ladder) if ladder else []
        self._index_cache = "" if index_cache is None else os.fspath(index_cache)
        self._use_cache = use_cache
        self._rebuild = rebuild

    def load(self, *, log: Optional[Callable[[str], None]] = None) -> "IndexSolver":
        ok = self._impl.load(
            database_path=self._database_path,
            database=self._database,
            quad_tolerance=self._quad_tolerance,
            ladder=self._ladder,
            index_cache=self._index_cache,
            use_cache=self._use_cache,
            rebuild=self._rebuild,
            log=log,
        )
        if not ok:
            raise NoDatabaseError(
                "no usable star database or index. Looked in: "
                + (self._database_path or ", ".join(find_databases()))
            )
        self._loaded = True
        return self

    def solve(
            self,
            image,
            *,
            fov: float = 0.0,
            max_stars: int = 500,
            min_star_size: float = 1.5,
            downsample: int = 0,
            want_sip: bool = False,
            refine: bool = True,
            label: str = "",
            progress: Optional[Callable[[str], None]] = None,
    ) -> Solution:
        """Solves one image.

        `image` is an array of pixel values — 2-D, or 3-D as (colours, height,
        width) or (height, width, colours) — or a path, when this build reads
        image files. `fov` is the field diameter in degrees when known: it only
        orders the search, so a wrong value costs time rather than a solution.
        """
        self._require_ready()
        kwargs = dict(
            fov=fov,
            max_stars=max_stars,
            min_star_size=min_star_size,
            downsample=downsample,
            want_sip=want_sip,
            refine=refine,
            label=label,
            progress=progress,
        )
        if isinstance(image, (str, os.PathLike)):
            image, meta, kwargs = self._solve_file(image, **kwargs)
            if not kwargs.get("fov"):
                kwargs["fov"] = _fov_from_meta(meta, image)
        return from_outcome(self._impl.solve_array(image, **kwargs))

    # --- what is resident -----------------------------------------------------
    @property
    def tier_count(self) -> int:
        return self._impl.tier_count

    @property
    def quad_count(self) -> int:
        return self._impl.quad_count

    @property
    def densities(self) -> list:
        """The depth tiers of the ladder, in stars per square degree."""
        return list(self._impl.densities)

    @property
    def bytes(self) -> int:
        """Size of the index it addresses. Mapped, so not all of it is resident."""
        return self._impl.bytes

    @property
    def cache_path(self) -> str:
        return self._impl.cache_path

    def __repr__(self) -> str:
        if not self.ready:
            return "<IndexSolver not loaded>"
        return (
            f"<IndexSolver {self.database} {self.tier_count} tiers "
            f"{self.quad_count / 1e6:.1f}M quads {self.bytes / 2 ** 30:.1f} GiB>"
        )


class SpiralSolver(_Solver):
    """The faithful port of ASTAP's own search.

    Reproduces `astap_cli`'s answer digit for digit, which is what it is for.
    It searches outward from a start position, so it is much quicker when told
    roughly where the telescope was pointing and much slower when not.
    """

    _impl_factory = staticmethod(_core.SpiralSolver)

    def __init__(
            self,
            *,
            database_path: Optional[PathLike] = None,
            database: str = "auto",
            warm: bool = False,
            threads: int = 0,
    ):
        super().__init__(database_path=database_path, database=database)
        self._warm = warm
        self._threads = threads

    def load(self, *, log: Optional[Callable[[str], None]] = None) -> "SpiralSolver":
        ok = self._impl.load(
            database_path=self._database_path,
            database=self._database,
            warm=self._warm,
            threads=self._threads,
            log=log,
        )
        if not ok:
            raise NoDatabaseError(
                "no star database found. Looked in: "
                + (self._database_path or ", ".join(find_databases()))
            )
        self._loaded = True
        return self

    def solve(
            self,
            image,
            *,
            ra: Optional[float] = None,
            dec: Optional[float] = None,
            fov: float = 0.0,
            radius: float = 180.0,
            max_stars: int = 500,
            quad_tolerance: float = 0.007,
            min_star_size: float = 1.5,
            downsample: int = 0,
            slow: bool = False,
            check_pattern_filter: bool = False,
            want_sip: bool = False,
            label: str = "",
            progress: Optional[Callable[[str], None]] = None,
    ) -> Solution:
        """Solves one image, searching outward from (`ra`, `dec`) in radians.

        Without a start position the search begins wherever the file's header
        said, and covers `radius` degrees around it. `fov` is the image height
        in degrees; given, it stops the header's own plate scale being used.
        """
        self._require_ready()
        kwargs = dict(
            ra=99999.0 if ra is None else ra,
            dec=99999.0 if dec is None else dec,
            fov=fov,
            fov_specified=bool(fov),
            radius=radius,
            max_stars=max_stars,
            quad_tolerance=quad_tolerance,
            min_star_size=min_star_size,
            downsample=downsample,
            force_oversize=slow,
            check_pattern_filter=check_pattern_filter,
            want_sip=want_sip,
            label=label,
            progress=progress,
        )
        if isinstance(image, (str, os.PathLike)):
            image, meta, kwargs = self._solve_file(image, **kwargs)
            # Unlike the index solver, the port searches from somewhere, so what
            # the file said about where it was pointing is worth a great deal.
            if ra is None and meta["ra0"]:
                kwargs["ra"] = meta["ra0"]
            if dec is None and meta["dec0"]:
                kwargs["dec"] = meta["dec0"]
        return from_outcome(self._impl.solve_array(image, **kwargs))

    @property
    def areas_warmed(self) -> int:
        return self._impl.areas_warmed

    def __repr__(self) -> str:
        if not self.ready:
            return "<SpiralSolver not loaded>"
        return f"<SpiralSolver {self.database} in {self.database_path}>"


def _fov_from_meta(meta, pixels) -> float:
    """Field height in degrees from a file's own plate scale, when it has one."""
    cdelt2 = meta.get("cdelt2") or 0.0
    if not cdelt2:
        return 0.0
    height = pixels.shape[-2]
    return abs(cdelt2) * height
