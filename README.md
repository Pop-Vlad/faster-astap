# Faster ASTAP plate solver — C++ port, and an index solver built on it

A plate solver: give it an astronomical image and it tells you where on the sky
the camera was pointed, at what scale and rotation.

This is a C++17 port of ASTAP's astrometric solving algorithm by Han Kleijn
(www.hnsky.org), translated from the Pascal sources in this repository,
plus an optimized solver built from the same pieces. Two binaries come out of the
build:

| | |
| --- | --- |
| **`astap_solve`** | The faithful port. Reproduces `astap_cli`'s answer **bit for bit**, 8x faster. Use it when you need the reference result. |
| **`astap_index_solve`** | The index solver. Replaces the sky search with a pre-built whole-sky quad index, so a solve is milliseconds and does not depend on knowing where the telescope was pointed. Use it for throughput, blind solves, and wide fields. |

`astap_solve` is the one to reach for when you need the reference answer
digit for digit; `astap_index_solve` does everything else faster.

Both read the same star databases as the original, take the same options, and
write the same `.ini` and `.wcs` files, so either can stand in for `astap_cli`.

## Quick start

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j
ctest --test-dir build
```

On Windows, see [Building on Windows](#building-on-windows): the same CMake, in
an MSYS2 UCRT64 shell, plus one `cmake --install` that produces binaries which
run on any Windows 10 machine with nothing installed.

No dependencies beyond the standard library. You also need a star database —
the same `.1476`, `.290` or `.001` files the original uses, from
[www.hnsky.org](https://www.hnsky.org/astap.htm). D80 (1.33 GB) is recommended and suits most
images; D50 or W08 are smaller and shallower.

```sh
# the port: blind solve, searching the whole sky
./build/astap_solve -f image.fits -fov 0 -r 180 -d /path/to/database

# the index solver: same thing, no start position needed at all
./build/astap_index_solve -f image.fits -d /path/to/database

# many images at once — the index is loaded once and reused
./build/astap_index_solve -d /path/to/database *.fits
```

The first `astap_index_solve` run builds the index (about 5 s) and caches it
under `~/.cache/faster-astap/` (`%LOCALAPPDATA%\faster-astap\cache\` on
Windows), one file per star database. Later runs read it back. The solution lands
in `image.ini`, in the same format `astap_cli` writes;
add `-wcs` for a FITS-header `.wcs` file. Exit status is 0 solved, 1 no
solution, 2 not enough stars, 16 image read error, 32 no star database, 33 star
database read error.

### Image formats

The same extensions `astap_cli` accepts, read in process — no conversion step,
no temporary files:

| | |
| --- | --- |
| `.fit` `.fits` `.fts` `.new` | uncompressed FITS, BITPIX 8/16/32/-32/-64 |
| `.fz` | compressed FITS: Rice (`RICE_1`) and, unlike `astap_cli`, GZIP |
| `.ppm` `.pgm` `.pfm` | binary Netpbm (`P5`/`P6`) and Portable Float Map |
| `.bmp` | 1/4/8/16/24/32 bit, uncompressed |
| `.png` `.jpg` `.jpeg` `.tif` `.tiff` | through libpng, libjpeg and libtiff |
| `.cr2` `.cr3` `.nef` `.arw` `.dng` … | raw camera files through LibRaw |

The four libraries are optional: CMake links whichever it finds and the binary
reports the rest as unsupported (`astap_solve -h` lists what your build reads).
FITS, Netpbm and BMP never need anything beyond the compiler.

Raw files are handed to the solver as the undemosaiced sensor frame, margins
included, which is what the ASTAP GUI gets from `dcraw`/`unprocessed_raw` —
byte for byte the same pixels, without the detour through a PGM.

## Performance

Measured on a 13600K (14 cores / 20 threads), solving a
6020×4015 frame against the D80 database. Best of three runs with a warm page cache.

"Blind" starts the search **90° from the true position** — a deliberately
average blind solve, since a start position that happens to be close flatters
the spiral search and one that is antipodal punishes it. "Hinted" gives the
field size and a position within 5°, which is the normal case when a mount
reports where it is pointing.

| | blind (90° off) | vs `astap_cli` | hinted | vs `astap_cli` |
| --- | --- | --- | --- | --- |
| `astap_cli` 2026.05.18 (reference) | 24.39 s | 1.0x | 0.341 s | 1.0x |
| `astap_solve`, 1 thread | 9.86 s | 2.5x | 0.180 s | 1.9x |
| `astap_solve`, 20 threads | **3.07 s** | **7.9x** | **0.158 s** | **2.2x** |
| `astap_index_solve`, per run | 1.65 s | 14.8x | 1.63 s | 0.2x |
| `astap_index_solve`, per image in a batch | **0.005 s** | **4900x** | **0.002 s** | **170x** |

Reading the two index solver rows together matters, and the second one is the
honest headline only if you are solving more than one image. The solve itself is
5 ms; the 1.6 s is almost entirely the 2.7 GB index being read from disk, paid
once per invocation. Hand it eight images and it solves all eight in 2.08 s
total. Hand it one, and for a *hinted* solve it is slower than the reference —
the index buys you nothing when the answer was nearly known already.

Things worth noting:

* **The index solver does not care where the telescope was pointed.** Its blind
  and hinted times are the same, because it never searches the sky — it looks
  the image's quads up in a whole-sky table. The spiral search, by contrast,
  costs whatever the distance from the start position happens to be: 0.16 s when
  the hint is good, 3.1 s at 90°, and worse further out.
* **The port is bit-identical to `astap_cli`.** Every optimisation behind the
  7.9x was checked to leave the solution unchanged to the last digit, at every
  thread count.

### Capability, not just speed

Speed is only half the question — a solver that is fast because it quietly gives
up on hard images is not faster. Over a 40 image corpus fetched from NASA
SkyView (`tools/fetch_skyview_corpus.py`), spanning 0.5° to 10° fields, the
galactic plane to the galactic pole, both celestial poles and the RA 0h
wraparound, with ground truth taken from each file's own WCS:

| | |
| --- | --- |
| `astap_solve` solved | 28/40 |
| `astap_index_solve` solved | **40/40** |
| images the port solves that the index misses | **0** |
| position error against the sky | median 0.15 px, worst 2.2 px |

The twelve extra are the 5° and 10° fields, where the spiral's auto-FOV sweep
gives up. `tools/corpus_harness.cpp` runs both solvers and exits non-zero if the
index solver ever loses an image the port could solve; that is the acceptance
gate for any change to it.

## Usage

### `astap_solve` — the port

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

The solution is written to `<output>.ini` (and `<output>.wcs` with `-wcs`), in
the same format as `astap_cli`. Exit status: 0 solved, 1 no solution, 2 not
enough stars, 16 image read error, 32 no star database, 33 star database read
error.

### `astap_index_solve` — the index solver

Takes the same options, plus the index ones. `-ra` / `-spd` / `-r` are absent
because it does not search from a start position, and `-sip` because it produces
a linear WCS only.

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
-rebuild  rebuild the index even when a usable cache exists
-nocache  build in memory, neither read nor write a cache
-cacheinfo  report the cache that would be used, then exit
```

Passing several images to one invocation is much faster than one invocation per
image: the index is read once.

#### Depth tiers, and making the index smaller

The index is a ladder of depth tiers, and which tier solves an image depends on
how many stars the image yielded per square degree. That is not a tuning detail,
it is the central constraint: a database quad is only findable when all four of
its stars were bright enough to be *detected in the image*, so an index built too
deep holds quads whose stars the image never saw, and one built too shallow has
too few. A tier reaches about a factor of two in image density either side of
itself, which is why the default ladder steps by 2.5x.

Whole sky, D80, on a 13600K — 12 tiers, 61.6 M quads, 3.1 GB resident (2.7 GB on
disk), 5.1 s to build:

| tier (stars/deg²) | quads | size | | tier | quads | size |
| --- | --- | --- | --- | --- | --- | --- |
| 0.5 | 0.14 M | 19 MB | | 32 | 0.98 M | 59 MB |
| 1 | 0.33 M | 28 MB | | 60 | 1.9 M | 101 MB |
| 2 | 0.73 M | 47 MB | | 125 | 3.9 M | 198 MB |
| 4 | 1.5 M | 85 MB | | 250 | 7.8 M | 385 MB |
| 8 | 0.24 M | 24 MB | | 500 | 15.6 M | 759 MB |
| 16 | 0.49 M | 35 MB | | 900 | 28.0 M | 1358 MB |

The four sparsest rungs are the ones built with the larger quad groups described
below, which is why they do not fall off smoothly. The two deepest rungs hold
two thirds of the whole ladder, so if you know your images are not deep, `-tiers`
cuts the index down sharply for free. A typical amateur frame yields 100-500
stars/deg²:

```sh
# 1.28 GB instead of 2.7 GB, covering roughly 30-1000 stars/deg²
astap_index_solve -f image.fits -d <db> -tiers 60,125,250,500
```

That halves the cache and, because the per-run cost is dominated by reading it,
takes a single solve from 1.65 s to 0.91 s. Dropping the 900 rung is where
almost all of the saving is — trimming only the sparse end is not worth doing,
since those rungs are tiny and are exactly what lets it solve 5-10° fields.
Changing `-tiers` invalidates the cache, so the next run rebuilds it (2.3 s
here).

#### The second pass against the database, and SIP

Once the field is known, the solver reads the star database *once* at that
position — at a depth matched to the image's own star count, over a square large
enough to take in the whole frame — and redoes the match there. This is the same
work the spiral does at each of its many positions, done once at the right one.

End to end, including rebuilding the image quads, the match and the cubic fits.
It costs around 0.5ms. Against a 5 ms solve that is ~18%. `-progress` reports it per image.

It buys two things. Accuracy, because the index's own fit rests on quad centres
alone: over the corpus the median error drops from 0.214 px to 0.145 px, and the
worst individual gains are large (a 10° field went from 0.67 px to 0.07 px). And
SIP, because a cubic distortion fit needs at least twenty matched quads, which
the raw index consensus reaches on 13 of 40 corpus images and the second pass
lifts to 19.

`-norefine` turns it off. `-sip` adds the coefficients, which go to the `.wcs`
file (the `.ini` format has no place for them) and set `CTYPE1/2` to
`RA---TAN-SIP`. When there are too few quads the solver says so and writes a linear WCS.

#### Why it solves wide fields the port cannot

Quads are normally built from each star's three nearest neighbours, one per star.
At a few stars per square degree that is fragile: a quad matches only when the
image and the catalogue picked the same four stars, so a single detection the
catalogue subset lacks changes which three neighbours a star has and replaces its
quad outright. A 10° field at 1 star/deg² yields only 76 quads, with no
redundancy to absorb that.

So when the ordinary pass finds nothing, the solver rebuilds the quads from every
combination of each star's six nearest — 15 per star, 1024 quads on that field
instead of 76, which took its true matches from 3 to 15. It is a strict superset,
so it can only add matches; it is held back as a retry purely because it costs
fifteen times the index queries. Three of the 40 corpus images need it, all at
1-8 stars/deg² and 5-10°.

The cache lives at `~/.cache/faster-astap/<database>_<type>_t<tolerance>.qix`
(`$XDG_CACHE_HOME` is honoured), one file per star database, so several
databases can coexist. It is ~2.7 GB for the default 12 tier ladder over D80. A
cache written by a different version, byte order, tolerance or ladder is
rejected and rebuilt rather than partially trusted, so a stale one fails loudly
instead of matching against the wrong thing.

A star database is required — the same `.1476`, `.290` or `.001` files the
original uses, downloadable from www.hnsky.org.

## Build

The plain build is in [Quick start](#quick-start) and needs nothing beyond a C++17 compiler.

CMake targets: `astap_solve`, `astap_index_solve`, `corpus_harness`,
`solver_tests`, `image_io_tests`, `quad_batch_tests`, `quad_index_bench`.

Two libraries are built: `astap_solver` is the solving core, `astap_image`
everything that turns a file into pixels. Only the image module looks for
optional dependencies, and it reports what it found:

```
-- Image formats: FITS, FITS.fz (Rice), Netpbm, BMP, FITS.fz (GZIP), PNG, JPEG, TIFF, raw camera files
```

zlib, libpng, libjpeg, libtiff and LibRaw are each picked up when present and
skipped when not; `-DASTAP_PNG=OFF` and friends leave one out on purpose. A
build without any of them still reads FITS, `.fz` Rice, Netpbm and BMP.

### Building on Windows

The toolchain is [MSYS2](https://www.msys2.org/)'s **UCRT64** environment: the
same GCC as the Linux build, and the environment where `pkg-config` finds all
five optional libraries, LibRaw included. Install MSYS2, then from the
**MSYS2 UCRT64** shell:

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

The test binaries are built in `build` and are not part of the bundle, so this
is the one job that uses the build tree directly:

```sh
ctest --test-dir build
```

`image_io_tests` needs the reference files in `tests/data`, which are not in the
repository, and reports `Error, accessing the file!` for each one that is
missing. Write them once and the suite is green:

```sh
python tools/make_test_images.py tests/data
```

That needs astropy and Pillow, which MSYS2 does not package for UCRT64. The
simplest source is a normal Windows Python, from outside the MSYS2 shell:

```pwsh
py -m venv .venv
.venv\Scripts\python.exe -m pip install astropy pillow
.venv\Scripts\python.exe tools\make_test_images.py tests\data
```

`quad_batch_tests` compares the batched quad construction against
`find_quads()` value by value with `memcmp`, across every group size the solver
uses, so the batching cannot change a solution.

Its reference files are written by `tools/make_test_images.py`; see
[Building on Windows](#building-on-windows) for the one time setup, which is
the same on any platform.

`image_io_tests` checks every loader against reference files decoded by other
implementations (astropy, i.e. CFITSIO, for the compressed FITS variants, and
Pillow for PNG/TIFF/JPEG), pixel by pixel and including the row order.

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

Steps 1-5 are the same for both solvers. They differ in step 6, in how they find
the database quads to compare against.

**`astap_solve` walks a squared spiral** over the sky from the start position,
reading the database and rebuilding its quads at every position until a match is
found, then re-solves once from the found position for maximum accuracy. This is
the original algorithm. Its cost is the distance from the start position, and a
blind solve on this machine issues about 273 000 database queries.

**`astap_index_solve` looks the quads up.** The hashes of the whole sky are built
once into an index, so a solve is a single pass over the image's ~90 quads
against it. The catch is that a quad is only findable when all four of its stars
were bright enough to be detected in the image, and one index commits to one
depth — so the index is a *ladder* of depth tiers from 0.5 to 900 stars/deg², and
the solver sweeps it. Candidate matches then vote jointly on plate scale and sky
position, and the winning cluster has to survive a RANSAC consensus before it
reaches the same least squares fit and the same acceptance checks the port uses.

Hash ratios are stored as fp32 in the index: over two million quads the worst
deviation from the fp64 value is 3.2e-7, against a matching tolerance of 0.007 —
a margin of more than 20 000x — and the final fit is fp64 throughout.

## Code Layout

| File                       | Contents                                                                                       |
|----------------------------|------------------------------------------------------------------------------------------------|
| `include/astap/types.h`    | `ImageArray`, `RowList` (the Pascal `Timage_array` / `Tstar_list`), `Header`                   |
| `astro_math.{h,cpp}`       | median, angular separation, gnomonic projection both ways, position angle, sexagesimal parsing |
| `image/image_io.{h,cpp}`   | `load_image`, the dispatch on the extension and the synthetic FITS header                      |
| `image/fits.{h,cpp}`       | FITS reader, plain and tile compressed, and the FITS header writer                             |
| `image/rice.{h,cpp}`       | Rice and GZIP tile decoding, from CFITSIO's `ricecomp.c` by way of the Pascal                  |
| `image/image_pnm.cpp`      | Netpbm and Portable Float Map, no dependencies                                                 |
| `image/image_bmp.cpp`      | Windows bitmap, no dependencies                                                                |
| `image/image_{png,jpeg,tiff,raw}.cpp` | the optional decoders, one system library each                                      |
| `star_detection.{h,cpp}`   | steps 1–2: `get_hist`, `get_background`, `hfd`, `find_stars`                                   |
| `quads.{h,cpp}`            | steps 3–5: `find_quads`, `find_many_quads`                                                     |
| `matching.{h,cpp}`         | steps 6–8: `find_fit`, `find_fit_using_hash`, `lsq_fit`, `find_offset_and_rotation`            |
| `star_database.{h,cpp}`    | `.290` / `.1476` / `.001` database reader, tile selection                                      |
| `calc_trans_cubic.{h,cpp}` | cubic fit behind the SIP coefficients                                                          |
| `solver.{h,cpp}`           | `Solver::solve` — spiral search, WCS derivation, `.ini` / `.wcs` output                        |
| `parallel.{h,cpp}`         | thread pool and range splitting used by the parallel stages                                    |
| `src/main.cpp`             | `astap_solve`, modelled on `astap_command_line.lpr`                                            |
| `quad_index.{h,cpp}`       | the whole-sky quad index: tier ladder build, binned query, on-disk cache                       |
| `index_solver.{h,cpp}`     | joint scale/position vote, RANSAC consensus, fit                                               |
| `src/index_main.cpp`       | `astap_index_solve`                                                                            |
| `quad_batch.{h,cpp}`       | quad construction for a whole batch of search positions, over the thread pool                  |
| `tools/corpus_harness.cpp` | runs both solvers over a corpus and reports the capability gate                                |
| `tools/fetch_skyview_corpus.py` | downloads the test corpus from NASA SkyView, with ground truth                            |
| `tools/png_to_fits.py`     | PNG/TIFF/JPEG to 16 bit FITS, for a build without the image libraries                          |
| `tools/make_test_images.py` | writes the PNG/TIFF/JPEG reference files `image_io_tests` reads from `tests/data`             |

## Verification

`solver_tests` covers the numerics (projection round trip, `lsq_fit`, `smedian`,
sexagesimal parsing), star detection on a synthetic field with planted
gaussians, quad matching that has to recover a known rotation/scale/offset, the
`find_many_quads` path for sparse fields, and the database tile numbering.

The port's own output is the tighter check: its blind and hinted `.ini` must stay
**bit identical** across refactors, and every optimisation in this README was
accepted only after that held. Save an `.ini` before a change and `diff` it
after; only the `CMDLINE=` line may differ.

### The capability gate for the index solver

Bit-identity is not available for the index solver — it is a different algorithm
and will not agree to the last digit. Its gate is statistical instead, and it is
a subset relation rather than a solve rate: **it does not have to solve images
the port cannot, but it must not lose any the port can.**

```sh
python3 tools/fetch_skyview_corpus.py --out corpus   # ~25 min, resumable
./build/corpus_harness corpus /path/to/database --csv results.csv
```

Every image carries the WCS SkyView wrote when it made the cutout, so each solve
is checked against the sky rather than against the other solver. The harness
exits non-zero if the index solver misses anything the port solved. It currently
reports 28/40 for the port, 40/40 for the index solver, and 0 gate failures.

## How the port was made fast

Every optimisation below was checked to produce a **bit-identical** solution to
the unoptimised port, at every thread count. The starting point was a direct
translation of the Pascal that solved the sample blind in 11.9 s; the changes
below and the parallel search took that to about 3 s.

A blind solve and a hinted solve stress completely different code. `perf` on the
original port:

* blind — `find_many_quads` 47%, `find_quads` 16%, `StarDatabase::read_star` 15%:
  the spiral search rebuilds the database quads at every sky position.
* hinted — `load_fits` 49%, `bin_mono_and_crop` 16%: the image never gets read
  more than once, so ingest dominates.

What was changed:

1. **`find_many_quads` distance table.** All C(7,4) = 35 quads of a group are
   built from the same 7 stars, so their six side lengths are a subset of the
   group's 21 pairwise distances. Computing those once per star instead of six
   per quad cuts 210 square roots per star down to 21. 172 → 56 µs per call at
   60 database stars.
2. **Duplicate quad rejection in O(1).** The original scans every quad found so
   far. A uniform grid with cell size 1.0 gives the *same* answer — two centres
   with |dx| < 1 and |dy| < 1 are at most one cell apart, so probing the 3×3
   neighbourhood is exactly equivalent. The grid is two flat arrays (bucket
   heads plus a next-index chain) reused across calls; a first attempt using a
   vector per bucket was *slower* than the linear scan it replaced, because this
   runs once per spiral position.
3. **Parallel spiral search**, the main win for blind solves. Positions are
   evaluated in batches; the lowest index that solves wins. That is exactly what
   the sequential search returns, since it stops at the first solving position
   and batches are generated in spiral order. Each worker has its own database
   file handle, tile cache and match state. The batch starts at one position and
   grows, so a solve that succeeds immediately never pays for a parallel batch.
4. **`load_fits`**: the per-pixel `switch (BITPIX)` is hoisted out of the pixel
   loop and a whole plane is read at once, which lets the byte swap vectorise;
   the conversion then runs in parallel. `measured_max` is a max-reduction, so
   it is order independent. The one order-dependent case — a NaN taking the
   running maximum — is detected and redone sequentially.
5. **`bin_mono_and_crop`**: restructured so each source row is walked once,
   sequentially, with a row accumulator. Every output pixel still sums its
   inputs in the original order (colour, row, column), so the result is
   unchanged. Parallel over output rows.
6. **`get_hist`**: per-thread partial histograms. The pixel sum is accumulated
   as `int64` rather than `double`, which is both exact and independent of the
   thread count.
7. **Allocation removal in the inner loops**: `hfd()` sorted its background
   annulus through a heap-allocated copy on every star candidate (now in place
   on the stack, as the Pascal does), and the sigma-clipped mean reallocated a
   65500-bin histogram per image tile.

Thread scaling, blind solve from 90° off, best of three:

| threads | 1    | 2    | 4    | 8    | 14   | 20   |
|---------|------|------|------|------|------|------|
| seconds | 9.86 | 7.13 | 5.01 | 3.26 | 3.13 | 3.47 |

Past about 8 threads this flattens and stops being monotonic — a repeat of the
sweep can put 20 threads ahead of 14. That is not measurement noise alone. The
spiral evaluates positions in batches sized from the thread count, so the count
decides which batch straddles the answer and how much speculative work that
batch does before one of its positions solves. On the throughput-bound part of
the curve (1 to 8 threads) scaling is clean.

The default is one thread per thread the process may actually run on,
SMT siblings included — the spiral search gives every sibling independent work,
so the last step from 14 to 20 threads is still worth about 14%. The count comes
from the CPU affinity mask rather than from the size of the machine, so it stays
correct under `taskset`, inside a container, or on a batch scheduler that pins
jobs to a subset of the CPUs. `-threads N` overrides it.
