#!/usr/bin/env python3
"""Download a parametric test corpus from NASA SkyView.

Every image comes back as FITS with a correct WCS written by SkyView, so the
ground truth travels with the file - that is what makes this usable as a
capability gate rather than a self-consistency check.

The grid is chosen to exercise the things that actually break a plate solver,
per the index solver sections of README.md:

  * field of view from 0.05 to 10 degrees, crossing the auto-FOV range
    (0.38 to 9.5) and the database-type switch at 6 degrees;
  * survey depth, fetched as two separate images of the same field below one
    degree: a shallow one (2MASS, a few hundred stars/deg^2) and a deep one
    (DSS2 Red, thousands). Depth is not a detail of the cutout, it is the axis
    the index solver lives on - the tier that solves an image has to match its
    detected density, and a cutout at fixed pixel count gets *deeper* as the
    field shrinks (the same field gives 748 stars/deg^2 at 0.5 degrees and 2700
    at 0.1), which is what puts small fields above the default tier ladder;
  * galactic latitude from the plane to the pole, which is the practical knob
    for detected star density - the axis that broke the quad index when the
    index tier and the image density diverged by more than about 2x;
  * a field near the celestial pole, where the tile geometry and the RA
    wraparound are special-cased.

Usage:
    python3 tools/fetch_skyview_corpus.py --out corpus
    python3 tools/fetch_skyview_corpus.py --out corpus --dry-run
    python3 tools/fetch_skyview_corpus.py --out corpus --fov 1 2 --limit 4

Then point the comparison harness at corpus/manifest.csv.

Only the standard library is used. SkyView is a free public service: the script
requests one image at a time with a pause between them, and skips anything it
has already downloaded, so an interrupted run can simply be repeated.
"""

import argparse
import csv
import os
import sys
import time
import urllib.parse
import urllib.request

ENDPOINT = "https://skyview.gsfc.nasa.gov/cgi-bin/images"

# name, RA (deg), Dec (deg), hard?, why it is in the set
#
# "hard" fields are excluded unless --include-hard is given. They are not merely
# difficult: a globular cluster blends its brightest stars into an unresolved
# core, so the quads are geometrically degenerate, and both this port and
# astap_cli fail on a 0.5 degree DSS2 cutout of M13. Keep them for regression
# purposes, not as part of the routine solve-rate figure.
FIELDS = [
    ("mid_lat_a", 53.0, 24.0, False, "ordinary mid galactic latitude field"),
    ("mid_lat_b", 160.0, -20.0, False, "ordinary mid galactic latitude field"),
    ("m31", 10.68, 41.27, False, "extended galaxy, moderate density"),
    ("galactic_pole", 192.86, 27.13, False, "north galactic pole, sparse"),
    ("sparse_south", 40.0, -55.0, False, "high galactic latitude, sparse"),
    ("near_ncp", 90.0, 88.5, False, "near the celestial pole, tile geometry edge case"),
    ("near_scp", 270.0, -88.5, False, "near the south celestial pole"),
    ("ra_zero", 0.3, 15.0, False, "straddles RA 0h, wraparound edge case"),
    ("galactic_plane", 284.0, 0.5, True, "Milky Way plane, very dense"),
    ("m13", 250.42, 36.46, True, "globular cluster, degenerate quads"),
    ("orion", 83.82, -5.39, True, "strong nebulosity"),
]

# Which surveys to fetch at a given field size. Below one degree both are taken,
# because depth is an independent axis there and the two answer different
# questions: 2MASS stands in for a short exposure (400-4000 stars/deg^2 over this
# grid) and DSS2 Red for a deep one (1500-6300). The same field at the same size
# can solve in one and not the other, so a single survey would report a solve
# rate that is really a statement about the cutout.
#
# Above one degree only DSS2 is fetched: 2MASS at 2 degrees and beyond is too
# sparse to be a fair wide-field test, and the wide end is not where depth binds.
# Below 0.1 degrees only 2MASS, since at 0.05 degrees every field is down to
# 3-10 detected stars and neither solver gets near it - those images are kept as
# the negative control that fixes where the floor is, and one survey shows that.
SURVEYS_FOR_FOV = [
    (0.1, ["2MASS-J"]),  # below 0.1 degrees
    (1.0, ["2MASS-J", "DSS2 Red"]),  # 0.1 to 1 degree: both depths
    (99.0, ["DSS2 Red"]),  # 1 degree and above
]
SURVEY_TAGS = {"2MASS-J": "2massj", "DSS2 Red": "dss2r", "DSS1 Red": "dss1r"}

DEFAULT_FOVS = [0.05, 0.1, 0.15, 0.25, 0.35, 0.5, 1.0, 2.0, 5.0, 10.0]


def surveys_for(fov):
    for limit, names in SURVEYS_FOR_FOV:
        if fov < limit:
            return names
    return SURVEYS_FOR_FOV[-1][1]


def build_url(survey, ra, dec, fov, pixels):
    params = {
        "Survey": survey,
        "Position": f"{ra},{dec}",
        "Coordinates": "J2000",
        "Projection": "Tan",
        "Size": f"{fov}",
        "Pixels": f"{pixels},{pixels}",
        "Return": "FITS",
    }
    return ENDPOINT + "?" + urllib.parse.urlencode(params)


def fetch(url, dest, timeout, retries=3):
    for attempt in range(retries):
        try:
            req = urllib.request.Request(
                url, headers={"User-Agent": "astap-cpp-platesolver-corpus/1.0"})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                data = r.read()
            # SkyView answers errors with an HTML page and HTTP 200.
            if not data.startswith(b"SIMPLE"):
                head = data[:200].decode("ascii", "replace").strip().replace("\n", " ")
                return False, f"not FITS ({head[:120]})"
            with open(dest, "wb") as f:
                f.write(data)
            return True, f"{len(data) / 1e6:.1f} MB"
        except Exception as exc:  # network, timeout, HTTP error
            if attempt == retries - 1:
                return False, str(exc)
            time.sleep(3 * (attempt + 1))
    return False, "unreachable"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="corpus", help="output directory")
    ap.add_argument("--fov", nargs="*", type=float, default=DEFAULT_FOVS,
                    help="field sizes in degrees")
    ap.add_argument("--pixels", type=int, default=1024,
                    help="image is pixels x pixels; with --fov this sets the plate scale")
    ap.add_argument("--surveys", nargs="*", default=None,
                    help="override the automatic survey choice (default: depth matched to FOV)")
    ap.add_argument("--fields", nargs="*", default=None, help="subset of the named fields")
    ap.add_argument("--include-hard", action="store_true",
                    help="also fetch the degenerate fields (globular cluster, galactic plane)")
    ap.add_argument("--limit", type=int, default=0, help="stop after this many downloads")
    ap.add_argument("--pause", type=float, default=2.0, help="seconds between requests")
    ap.add_argument("--timeout", type=float, default=180.0)
    ap.add_argument("--dry-run", action="store_true", help="list what would be fetched")
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    manifest_path = os.path.join(args.out, "manifest.csv")

    jobs = []
    for name, ra, dec, hard, note in FIELDS:
        if args.fields is not None and name not in args.fields:
            continue
        if hard and not args.include_hard and args.fields is None:
            continue
        for fov in args.fov:
            for survey in (args.surveys if args.surveys else surveys_for(fov)):
                tag = SURVEY_TAGS.get(survey, survey.replace(" ", "").lower())
                fname = f"{name}_{tag}_{fov:g}deg.fits"
                jobs.append({
                    "file": fname,
                    "field": name,
                    "survey": survey,
                    "ra_deg": ra,
                    "dec_deg": dec,
                    "fov_deg": fov,
                    "pixels": args.pixels,
                    "scale_arcsec_px": fov * 3600.0 / args.pixels,
                    "hard": hard,
                    "note": note,
                })

    if args.limit:
        jobs = jobs[:args.limit]

    print(f"{len(jobs)} images -> {args.out}/")
    if args.dry_run:
        for j in jobs:
            print(f"  {j['file']:44s} {j['scale_arcsec_px']:7.2f} \"/px  {j['note']}")
        print("\n(dry run, nothing downloaded)")
        return 0

    rows, ok, skipped, failed = [], 0, 0, 0
    for i, j in enumerate(jobs, 1):
        dest = os.path.join(args.out, j["file"])
        if os.path.exists(dest) and os.path.getsize(dest) > 0:
            print(f"[{i}/{len(jobs)}] {j['file']:44s} already present")
            rows.append(j)
            skipped += 1
            continue

        url = build_url(j["survey"], j["ra_deg"], j["dec_deg"], j["fov_deg"], j["pixels"])
        print(f"[{i}/{len(jobs)}] {j['file']:44s} ", end="", flush=True)
        good, detail = fetch(url, dest, args.timeout)
        if good:
            print(detail)
            rows.append(j)
            ok += 1
        else:
            print(f"FAILED: {detail}")
            if os.path.exists(dest):
                os.remove(dest)
            failed += 1
        time.sleep(args.pause)

    if rows:
        with open(manifest_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)

    print(f"\ndownloaded {ok}, already present {skipped}, failed {failed}")
    print(f"manifest: {manifest_path}")
    print("\nThe true position and scale of every image are in its own FITS header "
          "(CRVAL1/2, CDELT1/2), written by SkyView; the manifest repeats them for "
          "convenience. Compare a solver's output against those.")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
