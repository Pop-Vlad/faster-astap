# ASTAP plate solver — C++ port

A C++17 port of ASTAP's astrometric plate solving algorithm by Han Kleijn
(www.hnsky.org), translated from the Free Pascal sources in this repository.
The Pascal command line solver (`command-line_version/`) was used as the
reference.

## The method

```
      => Image <=                                    |  => Star database <=
 1) Find background, noise and star level            |
 2) Find stars and their CCD x, y position           | Extract the same amount of stars
    (standard coordinates)                           | (area corrected) from the area of
                                                     | interest, convert the α, δ equatorial
                                                     | coordinates into standard coordinates
 3) Build the smallest irregular tetrahedrons        | Idem
    ("quads") from four close stars, calculate the   |
    six distances and the mean x,y of the quad       |
 4) Sort the six distances, d1 longest, d6 shortest  | Idem
 5) Scale them as (d1, d2/d1 … d6/d1) — the hash code| Idem

                        => matching process <=
 6) Find hash code matches where the five ratios agree within a small tolerance.
 7) Take the median of the d1_found/d1_reference ratios, drop the outliers.
 8) Solve the overdetermined system A·S = X_ref for the six plate constants
    (least squares, Givens rotations), then derive the image centre position,
    the pixel scale and the rotation.
```

The database side of the search walks a squared spiral over the sky around the
start position until a match is found, then re-solves once from the found
position for maximum accuracy.

## Layout

| File                       | Contents                                                                                       |
|----------------------------|------------------------------------------------------------------------------------------------|
| `include/astap/types.h`    | `ImageArray`, `RowList` (the Pascal `Timage_array` / `Tstar_list`), `Header`                   |
| `astro_math.{h,cpp}`       | median, angular separation, gnomonic projection both ways, position angle, sexagesimal parsing |
| `fits.{h,cpp}`             | uncompressed FITS reader, FITS header writer                                                   |
| `star_detection.{h,cpp}`   | steps 1–2: `get_hist`, `get_background`, `hfd`, `find_stars`                                   |
| `quads.{h,cpp}`            | steps 3–5: `find_quads`, `find_many_quads`                                                     |
| `matching.{h,cpp}`         | steps 6–8: `find_fit`, `find_fit_using_hash`, `lsq_fit`, `find_offset_and_rotation`            |
| `star_database.{h,cpp}`    | `.290` / `.1476` / `.001` database reader, tile selection                                      |
| `calc_trans_cubic.{h,cpp}` | cubic fit behind the SIP coefficients                                                          |
| `solver.{h,cpp}`           | `Solver::solve` — spiral search, WCS derivation, `.ini` / `.wcs` output                        |
| `src/main.cpp`             | command line front end, modelled on `astap_command_line.lpr`                                   |

## Build

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build          # or: ./build/solver_tests
```

No dependencies beyond the standard library. Do not enable `-ffast-math`: the
port relies on IEEE-754 semantics and on the default round-half-to-even mode,
which is what FPC's `Round()` uses.

## Usage

```
astap_solve -f image.fits [options]

-f   filename {FITS file}
-r   radius_area_to_search[degrees]        default 180
-fov diameter_field[degrees]               0 = auto
-ra  right_ascension[hours]
-spd south_pole_distance[degrees]          (declination + 90)
-s   max_number_of_stars                   default 500
-t   quad_tolerance                        default 0.007
-m   minimum_star_size["]                  default 1.5
-z   downsample_factor[0,1,2,3,4,…]        0 = auto
-check   apply the check pattern filter (raw OSC images, binning 1x1 only)
-d   path to the star database
-D   star database abbreviation (d80, d50, g05, w08, …)
-o   base path & file name for the output files
-sip     add SIP (Simple Image Polynomial) coefficients
-speed   auto | slow
-wcs     also write a .wcs file in the FITS header format
-log     write the solver log to a .log file
-progress  log every search step
```

The solution is written to `<output>.ini` (and `<output>.wcs` with `-wcs`), in
the same format as `astap_cli`. Exit status: 0 solved, 1 no solution, 2 not
enough stars, 16 image read error, 32 no star database, 33 star database read
error.

A star database is required — the same `.1476`, `.290` or `.001` files the
original uses, downloadable from www.hnsky.org.

## Verification

`solver_tests` covers the numerics (projection round trip, `lsq_fit`, `smedian`,
sexagesimal parsing), star detection on a synthetic field with planted
gaussians, quad matching that has to recover a known rotation/scale/offset, the
`find_many_quads` path for sparse fields, and the database tile numbering.

End to end, the port was compared against `astap_cli` 2026.05.18 on the
6020×4015 sample in this directory (`test_img.png`, converted with
`tools/png_to_fits.py`) using the D80 database:

|                                      | astap_cli                       | this port               |
|--------------------------------------|---------------------------------|-------------------------|
| blind solve (`-fov 0 -r 180`)        | 18:51:14.9 +32°43'26" in 22.2 s | same position in 11.8 s |
| centre agreement                     |                                 | 0.045″ = 0.07 pixel     |
| pixel scale                          | 0.64192″/px                     | 0.64197″/px             |
| hinted solve (`-fov 0.72 -z 3 -r 5`) | 32/32 quads, 0.1 s              | 32/32 quads, 0.1 s      |
| centre agreement                     |                                 | 0.016″ = 0.025 pixel    |

The residual differences come from a slightly different set of detected stars
(127 vs 124 at the same settings), not from the geometry.

## Scope

Ported: the complete solving pipeline — star detection, quad construction and
matching, the least squares fit, the star database readers (`.290`, `.1476`,
`.001`), the spiral sky search, auto FOV, the second maximum-accuracy pass,
SIP distortion coefficients, and the `.ini` / `.wcs` output.

Not ported, because they are input handling rather than the solving algorithm:
Rice compressed (`.fz`) FITS, the TIFF/PNG/JPEG/RAW loaders (use
`tools/png_to_fits.py`), the `-analyse` / `-extract` reporting modes, the
`-update` in-place FITS header rewrite, and the star database record sizes
other than 5 and 6 bytes (astap_cli has the same restriction).
