# Faster ASTAP plate solver - C++ port, and an index solver built on it

A plate solver: give it an astronomical image and it tells you where on the sky the camera was pointed, at what scale
and rotation.

This is a C++17 port of ASTAP's astrometric solving algorithm by Han Kleijn (www.hnsky.org), translated from the Pascal
sources in this repository, plus an optimized solver built from the same pieces. Two binaries come out of the build:

|                         |                                                                                                                                                                                                                                   |
|-------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **`astap_solve`**       | The faithful port. Reproduces `astap_cli`'s answer **bit for bit**, 8x faster. Use it when you need the reference result.                                                                                                         |
| **`astap_index_solve`** | The index solver. Replaces the sky search with a pre-built whole-sky quad index, so a solve is milliseconds and does not depend on knowing where the telescope was pointed. Use it for throughput, blind solves, and wide fields. |

`astap_solve` is the one to reach for when you need the reference answer digit for digit; `astap_index_solve` does
everything else faster.

Both read the same star databases as the original, take the same options, and write the same `.ini` and `.wcs` files, so
either can stand in for `astap_cli`.

Both are also a [Python package](#python-package), `pip install faster-astap`, which is the way to use them from a
script or a notebook - and the only way to solve an image already in memory, without it going through a file first.

[DESIGN.md](DESIGN.md) has the rest: how the matching works, the depth ladder the index solver rests on, where each
solver stops working and why, and how the port was made fast.

## Quick start

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j
ctest --test-dir build
```

On Windows, see [Building on Windows](#building-on-windows): the same CMake, in an MSYS2 UCRT64 shell, plus one
`cmake --install` that produces binaries which run on any Windows 10 machine with nothing installed.

No dependencies beyond the standard library. You also need a star database - the same `.1476`, `.290` or `.001` files
the original uses, from
[www.hnsky.org](https://www.hnsky.org/astap.htm). D80 (1.33 GB) is recommended and suits most images; D50 or W08 are
smaller and shallower.

```sh
# the port: blind solve, searching the whole sky
./build/astap_solve -f image.fits -fov 0 -r 180 -d /path/to/database

# the index solver: same thing, no start position needed at all
./build/astap_index_solve -f image.fits -d /path/to/database

# many images at once - the index is loaded once and reused
./build/astap_index_solve -d /path/to/database *.fits
```

The first `astap_index_solve` run builds the index (about 5 s) and caches it under `~/.cache/faster-astap/`
(`%LOCALAPPDATA%\faster-astap\cache\` on Windows), one file per star database. Later runs read it back. The solution
lands in `image.ini`, in the same format `astap_cli` writes; add `-wcs` for a FITS-header `.wcs` file. Exit status is 0
solved, 1 no solution, 2 not enough stars, 16 image read error, 32 no star database, 33 star database read error.

### Image formats

The same extensions `astap_cli` accepts, read in process - no conversion step, no temporary files:

|                                      |                                                                |
|--------------------------------------|----------------------------------------------------------------|
| `.fit` `.fits` `.fts` `.new`         | uncompressed FITS, BITPIX 8/16/32/-32/-64                      |
| `.fz`                                | compressed FITS: Rice (`RICE_1`) and, unlike `astap_cli`, GZIP |
| `.ppm` `.pgm` `.pfm`                 | binary Netpbm (`P5`/`P6`) and Portable Float Map               |
| `.bmp`                               | 1/4/8/16/24/32 bit, uncompressed                               |
| `.png` `.jpg` `.jpeg` `.tif` `.tiff` | through libpng, libjpeg and libtiff                            |
| `.cr2` `.cr3` `.nef` `.arw` `.dng` … | raw camera files through LibRaw                                |

The four libraries are optional: CMake links whichever it finds and the binary reports the rest as unsupported
(`astap_solve -h` lists what your build reads). FITS, Netpbm and BMP never need anything beyond the compiler.

Raw files are handed to the solver as the undemosaiced sensor frame, margins included, which is what the ASTAP GUI gets
from `dcraw`/`unprocessed_raw` - byte for byte the same pixels, without the detour through a PGM.

## Performance

Measured on a 13600K (14 cores / 20 threads), solving a 6020×4015 frame against the D80 database. Best of three runs
with a warm page cache.

"Blind" starts the search **90° from the true position** - a deliberately average blind solve, since a start position
that happens to be close flatters the spiral search and one that is antipodal punishes it. "Hinted" gives the field size
and a position within 5°, which is the normal case when a mount reports where it is pointing.

|                                           | blind (90° off) | vs `astap_cli` | hinted      | vs `astap_cli` |
|-------------------------------------------|-----------------|----------------|-------------|----------------|
| `astap_cli` 2026.05.18 (reference)        | 24.39 s         | 1.0x           | 0.341 s     | 1.0x           |
| `astap_solve`, 1 thread                   | 9.86 s          | 2.5x           | 0.180 s     | 1.9x           |
| `astap_solve`, 20 threads                 | **3.07 s**      | **7.9x**       | **0.158 s** | **2.2x**       |
| `astap_index_solve`, per run              | 0.165 s         | 148x           | 0.153 s     | 2.2x           |
| `astap_index_solve`, per image in a batch | **0.005 s**     | **4900x**      | **0.002 s** | **170x**       |

Reading the two index solver rows together matters. The per-run figure is a whole invocation - map the index, read the
145 MB frame, detect stars, solve; the 5 ms is the solve stage alone, which is what a second image in the same
invocation adds. The 2.7 GB index is memory-mapped rather than read, so with a warm page cache starting a run costs
nothing worth measuring and eight images solve in 1.15 s total. The exception is a genuinely cold cache: the first run
after a boot pays about 4.9 s faulting the index in from disk, once.

Things worth noting:

* **The index solver does not care where the telescope was pointed.** Its blind and hinted times are the same, because
  it never searches the sky - it looks the image's quads up in a whole-sky table. The spiral search, by contrast, costs
  whatever the distance from the start position happens to be: 0.16 s when the hint is good, 3.1 s at 90°, and worse
  further out.
* **The port is bit-identical to `astap_cli`.** Every optimisation behind the 7.9x was checked to leave the solution
  unchanged to the last digit, at every thread count.

### Capability, not just speed

Speed is only half the question - a solver that is fast because it quietly gives up on hard images is not faster. The
corpus is 120 images from NASA SkyView (`tools/fetch_skyview_corpus.py`): eight fields from the galactic plane to the
galactic pole, both celestial poles and the RA 0h wraparound, cut at ten sizes from 0.05° to 10°. Below one degree every
field is taken twice, shallow (2MASS-J, 400-4000 detected stars/deg² here) and deep (DSS2 Red, 1500-6300), because depth
decides which index tier can solve an image and one survey would report a solve rate that is really a statement about
the cutout. Ground truth is the WCS SkyView wrote into each file. Everything below is blind - no position, no field
size.

|             | port   | index, default ladder | ceiling 3600 |
|-------------|--------|-----------------------|--------------|
| solved      | 69/120 | 80/120                | **92/120**   |
| index cache | -      | 2.7 GB                | 10.0 GB      |
| index build | -      | 5 s                   | 24 s         |

By field size, and this is the whole story in one table:

| field | images | port | index, default | index, ceiling 3600 |
|-------|--------|------|----------------|---------------------|
| 0.05° | 8      | 0    | 0              | 0                   |
| 0.1°  | 16     | 3    | 0              | 5                   |
| 0.15° | 16     | 2    | 3              | 10                  |
| 0.25° | 16     | 13   | 13             | 13                  |
| 0.35° | 16     | 15   | 16             | 16                  |
| 0.5°  | 16     | 16   | 16             | 16                  |
| 1°    | 8      | 7    | 8              | 8                   |
| 2°    | 8      | 8    | 8              | 8                   |
| 5°    | 8      | 4    | 8              | 8                   |
| 10°   | 8      | 1    | 8              | 8                   |

The two index columns are identical at every size from 0.25° up, which is the point
of [density matching](DESIGN.md#small-fields-and-where-the-ladder-ends): the default 2.7 GB cache reaches as far as the
10 GB one until the fields get small enough that thinning the image would leave too few stars. Buying the deep rungs is
worth it below 0.25° and pointless above.

The two solvers fail at opposite ends, for opposite reasons. At 5° and 10° the spiral's auto-FOV sweep gives up and the
index takes all sixteen. Below half a degree the port's blind sweep is guessing the field size off a geometric grid that
stops at 0.37°, so it degrades from 15/16 at 0.35° to 2/16 at 0.15°, while the index - which never needs the field
size - holds 16/16 down to 0.35° on the stock cache and is ahead of the port at every size below 0.5°. Accuracy is the
same either way: median 0.23 px for the index and 0.27 px for the port against the sky, and the worst cases are the DSS
plates' own astrometry rather than the fit.

Position error is quoted in pixels throughout, which for the small fields is the only honest unit: 60 arcsec is a sixth
of a 0.1° frame.

`tools/corpus_harness.cpp` runs both solvers and is the acceptance gate for any change to the index solver;
see [the capability gate](DESIGN.md#the-capability-gate-for-the-index-solver) for what it currently reports and what it
holds the solver to.

## Usage

### `astap_solve` - the port

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
-threads N   worker threads; omitted or 0 = one per available hardware thread
```

The solution is written to `<output>.ini` (and `<output>.wcs` with `-wcs`), in the same format as `astap_cli`. Exit
status: 0 solved, 1 no solution, 2 not enough stars, 16 image read error, 32 no star database, 33 star database read
error.

### `astap_index_solve` - the index solver

Takes the same options, plus the index ones. `-ra` / `-spd` / `-r` are absent because it does not search from a start
position, and `-sip` because it produces a linear WCS only.

```
astap_index_solve -f image.fits [options] [more.fits ...]

-f   filename {FITS file. May be repeated, or list files after the options}
-d   path to the star database (needed to build an index)
-D   star database abbreviation (d80, d50, g05, w08, …)
-fov diameter_field[degrees]               orders the depth sweep; does not restrict it
-s   max_number_of_stars                   default 500
-t   quad_tolerance                        default 0.007; a cache belongs to one tolerance
-m   minimum_star_size["]                  default 1.5, applied only with -fov
-z   downsample_factor[0,1,2,3,4,…]        0 = auto
-o   base path & file name for the output files
-sip      add SIP distortion coefficients (written to the .wcs file)
-norefine skip the second pass, leaving the raw index solution
-wcs / -log / -progress / -threads N       as above

-i        index cache file, overriding the default location
-tiers    d1,d2,… depth ladder in stars/deg² for a newly built index
-maxtier  d  raise the default ladder's ceiling to d stars/deg², for fields
             below 0.25°. Adds a rung of 2.5 GB at 1800 and 4.8 GB at 3600
-rebuild  rebuild the index even when a usable cache exists
-nocache  build in memory, neither read nor write a cache
-cacheinfo  report which rungs are cached, then exit
```

Passing several images to one invocation is much faster than one invocation per image: the index is read once.

An imaging application cannot do that - it solves one frame at a time, in a fresh process, minutes apart, so it pays the
index read every single solve and never sees the 5 ms. Keeping the index in memory between those invocations is what [
`nina-plugin/`](nina-plugin/README.md) is for. It is a module of its own and changes nothing here;
`-DASTAP_NINA_INTEGRATION=OFF` leaves it out of the build.

#### The index cache

The index is built once per star database and cached under
`~/.cache/faster-astap/` (`$XDG_CACHE_HOME` is honoured;
`%LOCALAPPDATA%\faster-astap\cache\` on Windows), one file per depth rung, about 2.7 GB in total for D80. Rungs are
independent files, so ladders compose: adding a deeper one leaves the rest alone, `-cacheinfo` lists what exists, and
deleting a rung's file reclaims its space.

Two options change what gets built, and both cost disk rather than accuracy:

```sh
# smaller: drop the rungs your images will never use. 1.28 GB instead of 2.7,
# covering roughly 30-1000 stars/deg^2, which suits a typical amateur frame
astap_index_solve -f image.fits -d <db> -tiers 60,125,250,500

# deeper: for fields below 0.25°, where the default ladder runs out of depth.
# Adds 7.3 GB, built and cached once
astap_index_solve -f image.fits -d <db> -maxtier 3600
```

Neither is needed for ordinary work - the default ladder covers 0.25° upwards.
[DESIGN.md](DESIGN.md#the-index-solvers-depth-ladder) explains what a rung is, why the ceiling sets the smallest field,
and what the measurements say.

## Python package

```sh
pip install faster-astap
```

You still need a star database - the package cannot bundle one, and looks in the same places the command line tools do:

```python
import faster_astap as fa

fa.find_databases()  # where it will look, in order
```

### Solving

The point of the package, rather than shelling out to `astap_index_solve`, is that the index stays in memory. `load()`
is the expensive call and happens once; after that a solve is milliseconds.

```python
import faster_astap as fa

solver = fa.IndexSolver(database_path=r"D:\astap\d80")
solver.load()  # reads the index; seconds the first time
print(solver)  # <IndexSolver d80 12 tiers 61.6M quads 2.9 GiB>

for frame in frames:  # a numpy array, or a path
    sol = solver.solve(frame, fov=1.2)
    if sol:
        print(sol.ra_deg, sol.dec_deg, sol.scale_arcsec, sol.rotation_deg)
```

`fov` is the field diameter in degrees when you know it. It only orders the depth sweep, so a wrong value costs time
rather than a solution, and leaving it out is fine - the index solver needs no position and no scale.

A `Solution` is falsey when nothing was found, so `if sol:` is the check. It carries the WCS in the units people
actually write down:

|                                      |                                                                   |
|--------------------------------------|-------------------------------------------------------------------|
| `ra_deg` `dec_deg` `ra_hours`        | image centre; `ra` and `dec` are the same in radians              |
| `scale_arcsec`                       | pixel scale, unsigned                                             |
| `rotation_deg`                       | rotation at the centre                                            |
| `fov_deg`                            | `(width, height)` in degrees, from the solved scale               |
| `stars` `nr_inliers` `solve_seconds` | what the solve did                                                |
| `messages`                           | the lines the command line tools would have printed               |
| `to_astropy_wcs()`                   | the solution as an `astropy.wcs.WCS`; astropy is not a dependency |

### Solving an array

Any 2-D array, or 3-D as `(colours, height, width)` or `(height, width, colours)`. Integer camera frames need no
conversion:

```python
from astropy.io import fits

frame = fits.getdata("light_0001.fits")  # uint16 straight off the sensor
sol = solver.solve(frame, fov=1.2)
```

### The port, for a reference answer

Same lifecycle, same `Solution`, so swapping one for the other is a one line change. It searches outward from a start
position, so it is much quicker told roughly where the telescope was pointing and much slower when not.

```python
import math

port = fa.SpiralSolver(database_path=r"D:\astap\d80").load()
sol = port.solve(frame,
                 ra=math.radians(53.0),  # radians, as the solver works in
                 dec=math.radians(24.0),
                 radius=10)  # degrees to search around it
```

Reach for the port when you need the answer `astap_cli` would have given, digit for digit. For everything else
`IndexSolver` is faster and needs less from you.

## Build

The plain build is in [Quick start](#quick-start) and needs nothing beyond a C++17 compiler.

CMake targets: `astap_solve`, `astap_index_solve`, `corpus_harness`,
`solver_tests`, `image_io_tests`, `quad_batch_tests`, `quad_index_bench`.

Two libraries are built: `astap_solver` is the solving core, `astap_image`
everything that turns a file into pixels. Only the image module looks for optional dependencies, and it reports what it
found:

```
-- Image formats: FITS, FITS.fz (Rice), Netpbm, BMP, FITS.fz (GZIP), PNG, JPEG, TIFF, raw camera files
```

zlib, libpng, libjpeg, libtiff and LibRaw are each picked up when present and skipped when not; `-DASTAP_PNG=OFF` and
friends leave one out on purpose. A build without any of them still reads FITS, `.fz` Rice, Netpbm and BMP.

### Building on Windows

The toolchain is [MSYS2](https://www.msys2.org/)'s **UCRT64** environment: the same GCC as the Linux build, and the
environment where `pkg-config` finds all five optional libraries, LibRaw included. Install MSYS2, then from the **MSYS2
UCRT64** shell:

```sh
pacman -S --needed \
  mingw-w64-ucrt-x86_64-gcc mingw-w64-ucrt-x86_64-cmake \
  mingw-w64-ucrt-x86_64-ninja mingw-w64-ucrt-x86_64-pkgconf \
  mingw-w64-ucrt-x86_64-zlib mingw-w64-ucrt-x86_64-libpng \
  mingw-w64-ucrt-x86_64-libjpeg-turbo mingw-w64-ucrt-x86_64-libtiff \
  mingw-w64-ucrt-x86_64-libraw
```

Then build, and install the result into `dist`:

```sh
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
cmake --install build --prefix dist
```

You will find the 2 built executables in the **`dist` directory:

```
astap_solve.exe
astap_index_solve.exe
```

#### Running the tests

The test binaries are built in `build` and are not part of the bundle, so this is the one job that uses the build tree
directly:

```sh
ctest --test-dir build
```

`image_io_tests` needs the reference files in `tests/data`, which are not in the repository, and reports
`Error, accessing the file!` for each one that is missing. Write them once and the suite is green:

```sh
python tools/make_test_images.py tests/data
```

That needs astropy and Pillow, which MSYS2 does not package for UCRT64. The simplest source is a normal Windows Python,
from outside the MSYS2 shell:

```pwsh
py -m venv .venv
.venv\Scripts\python.exe -m pip install astropy pillow
.venv\Scripts\python.exe tools\make_test_images.py tests\data
```

`quad_batch_tests` compares the batched quad construction against
`find_quads()` value by value with `memcmp`, across every group size the solver uses, so the batching cannot change a
solution.

Its reference files are written by `tools/make_test_images.py`; see
[Building on Windows](#building-on-windows) for the one time setup, which is the same on any platform.

`image_io_tests` checks every loader against reference files decoded by other implementations (astropy, i.e. CFITSIO,
for the compressed FITS variants, and Pillow for PNG/TIFF/JPEG), pixel by pixel and including the row order.

