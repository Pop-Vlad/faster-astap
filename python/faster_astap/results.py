"""The solution a solve returns, in the units a caller actually wants."""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Sequence


@dataclass(frozen=True)
class Solution:
    """Where the camera was pointed, and what it took to work that out.

    The solver works in radians and degrees-per-pixel because the algorithm
    does; the properties below are the same numbers in the units people write
    down. `head` is the raw WCS if you want it.
    """

    solved: bool
    errorlevel: int
    messages: Sequence[str]

    # Radians, as the solver produces them.
    ra: float = 0.0
    dec: float = 0.0
    # Degrees per pixel, signed as the WCS convention has them.
    cdelt1: float = 0.0
    cdelt2: float = 0.0
    crpix1: float = 0.0
    crpix2: float = 0.0
    crota1: float = 0.0
    crota2: float = 0.0
    cd: tuple = (0.0, 0.0, 0.0, 0.0)
    width: int = 0
    height: int = 0

    stars: int = 0
    stars_used: int = 0
    stars_detected: int = 0
    nr_inliers: int = 0
    binning: int = 1
    tier_density: float = 0.0
    tiers_tried: int = 0
    many_quads_pass: bool = False
    refined: bool = False
    solve_seconds: float = 0.0
    total_seconds: float = 0.0

    sip: Any = None
    cards: Sequence[str] = field(default_factory=tuple)

    # --- the units people write down -----------------------------------------
    @property
    def ra_deg(self) -> float:
        return math.degrees(self.ra)

    @property
    def dec_deg(self) -> float:
        return math.degrees(self.dec)

    @property
    def ra_hours(self) -> float:
        return math.degrees(self.ra) / 15.0

    @property
    def scale_arcsec(self) -> float:
        """Pixel scale in arcseconds, unsigned."""
        return abs(self.cdelt2) * 3600.0

    @property
    def rotation_deg(self) -> float:
        return self.crota2

    @property
    def fov_deg(self) -> tuple:
        """Field width and height in degrees, from the solved scale."""
        return (abs(self.cdelt1) * self.width, abs(self.cdelt2) * self.height)

    def __bool__(self) -> bool:
        return self.solved

    def __repr__(self) -> str:
        if not self.solved:
            why = self.messages[-1] if self.messages else "no solution"
            return f"<Solution unsolved: {why}>"
        return (
            f"<Solution ra={self.ra_deg:.5f}d dec={self.dec_deg:.5f}d "
            f"scale={self.scale_arcsec:.4f}\"/px rot={self.rotation_deg:.2f}d "
            f"quads={self.nr_inliers} in {self.solve_seconds:.3f}s>"
        )

    def wcs_cards(self) -> list:
        """The solution as FITS header cards.

        Present so a caller can hand them to astropy or write them out. Only
        available when the image came from a file, since the cards it is built
        on are the ones that file carried.
        """
        return list(self.cards)

    def to_astropy_wcs(self):
        """The solution as an `astropy.wcs.WCS`.

        astropy is not a dependency of this package; this raises ImportError if
        it is not installed.
        """
        from astropy.wcs import WCS  # imported here so it stays optional

        wcs = WCS(naxis=2)
        wcs.wcs.crpix = [self.crpix1, self.crpix2]
        wcs.wcs.crval = [self.ra_deg, self.dec_deg]
        wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        wcs.wcs.cd = [[self.cd[0], self.cd[1]], [self.cd[2], self.cd[3]]]
        return wcs


def from_outcome(outcome) -> Solution:
    """Builds a Solution from the extension's Outcome."""
    h = outcome.head
    return Solution(
        solved=outcome.solved,
        errorlevel=outcome.errorlevel,
        messages=tuple(outcome.messages),
        ra=h.ra0,
        dec=h.dec0,
        cdelt1=h.cdelt1,
        cdelt2=h.cdelt2,
        crpix1=h.crpix1,
        crpix2=h.crpix2,
        crota1=h.crota1,
        crota2=h.crota2,
        cd=(h.cd1_1, h.cd1_2, h.cd2_1, h.cd2_2),
        width=h.width,
        height=h.height,
        stars=outcome.stars,
        stars_used=outcome.stars_used,
        stars_detected=outcome.stars_detected,
        nr_inliers=outcome.nr_inliers,
        binning=outcome.bin,
        tier_density=outcome.tier_density,
        tiers_tried=outcome.tiers_tried,
        many_quads_pass=outcome.many_quads_pass,
        refined=outcome.refined,
        solve_seconds=outcome.solve_seconds,
        total_seconds=outcome.total_seconds,
        sip=outcome.sip if outcome.sip.valid else None,
        cards=tuple(h.cards),
    )
