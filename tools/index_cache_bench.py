#!/usr/bin/env python3
"""Measure what an index solve costs when the index is not in the page cache.

The index solver maps its tier files instead of reading them (src/mapped_file.cpp)
and asks for MADV_RANDOM (src/quad_index.cpp), so a tier page is fetched only
when a query reaches it. Warm, that is the whole point: a run touches a few per
cent of a multi-gigabyte cache and finishes in milliseconds. Cold, every touched
page is its own synchronous fault, and the run time becomes the fault count
times the drive's queue-depth-1 latency - which is why the interesting number
here is not seconds but *how much of the index a run touches*.

The failure path is what this tool exists for. A solve stops at the depth tier
that matches; a solve that finds nothing sweeps the whole ladder, so it touches
an order of magnitude more of the index than a success does. Whether that is
acceptable is a question about the cold case only, and the cold case cannot be
measured by rerunning the solver - the first run leaves the pages behind.

So each image is measured twice: best of N warm runs, then one run with the
index evicted from the page cache with posix_fadvise(DONTNEED). Major faults
come from the child's rusage, and the bytes actually pulled in come from
mincore() over the tier files, read after the solver has exited. The tool also
times a plain sequential read of the whole cache, which is the break-even any
demand-faulting scheme is measured against: once a run touches more of the index
than that, reading the file up front would have been cheaper.

Usage:
    python3 tools/index_cache_bench.py -d /path/to/star_database corpus/*.fits
    python3 tools/index_cache_bench.py -d db --warm-runs 5 --no-cold img.fits
    python3 tools/index_cache_bench.py -d db img.fits -- -fov 0.15 -norefine

Linux only: it needs posix_fadvise(DONTNEED) and mincore(). Nothing else may
hold the index mapped while it runs — stop astap_index_server first, or the
eviction silently does nothing and every run reports as warm.

Only the standard library is used.
"""

import argparse
import ctypes
import glob
import mmap
import os
import resource
import subprocess
import sys
import tempfile
import time

DEFAULT_CACHE = os.path.expanduser("~/.cache/faster-astap")
DEFAULT_SOLVER = "build/astap_index_solve"
PAGE = os.sysconf("SC_PAGE_SIZE")


def tier_files(cache):
    files = sorted(glob.glob(os.path.join(cache, "*.qix")))
    if not files:
        sys.exit("no index cache in %s — build one with a first solve, or pass --cache" % cache)
    return files


def evict(files):
    """Drop the index from the page cache. Clean, unmapped pages only."""
    for path in files:
        fd = os.open(path, os.O_RDONLY)
        try:
            os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED)
        finally:
            os.close(fd)


def resident_bytes(files, libc):
    """Bytes of the index in the page cache, via mincore() over a fresh mapping.

    The mapping is never touched, only asked about, so measuring does not itself
    fault anything in.
    """
    total = 0
    for path in files:
        size = os.path.getsize(path)
        fd = os.open(path, os.O_RDONLY)
        try:
            mm = mmap.mmap(fd, size, flags=mmap.MAP_PRIVATE)
        except ValueError:  # empty file
            os.close(fd)
            continue
        addr = ctypes.addressof(ctypes.c_char.from_buffer(mm))
        vec = (ctypes.c_ubyte * ((size + PAGE - 1) // PAGE))()
        if libc.mincore(ctypes.c_void_p(addr), ctypes.c_size_t(size), vec) != 0:
            raise OSError(ctypes.get_errno(), "mincore")
        total += sum(v & 1 for v in vec) * PAGE
        del addr  # the mapping cannot be closed while a pointer into it is live
        mm.close()
        os.close(fd)
    return total


def sequential_read(files):
    """Time streaming the whole cache: the cost demand faulting has to beat."""
    buf = bytearray(1 << 20)
    start = time.perf_counter()
    for path in files:
        with open(path, "rb", buffering=0) as f:
            while f.readinto(buf):
                pass
    return time.perf_counter() - start


def solve(cmd):
    """Run the solver once, returning wall time, major faults and the outcome."""
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    start = time.perf_counter()
    proc = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    elapsed = time.perf_counter() - start
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    return elapsed, after.ru_majflt - before.ru_majflt, proc.returncode


OUTCOMES = {0: "solved", 1: "no solution", 2: "too few stars",
            16: "image unreadable", 32: "no database", 33: "database unreadable"}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("images", nargs="+", help="image files to solve")
    ap.add_argument("-d", "--database", help="star database path, passed to the solver")
    ap.add_argument("--solver", default=DEFAULT_SOLVER, help="solver binary")
    ap.add_argument("--cache", default=DEFAULT_CACHE, help="index cache directory")
    ap.add_argument("--warm-runs", type=int, default=3,
                    help="warm runs to take the best of, after one priming run")
    ap.add_argument("--no-cold", action="store_true", help="skip the eviction pass")
    ap.add_argument("--no-baseline", action="store_true",
                    help="skip the sequential read of the whole cache")

    # Everything after a bare -- goes to the solver untouched, so its flags
    # cannot collide with this script's.
    argv = sys.argv[1:]
    extra = []
    if "--" in argv:
        cut = argv.index("--")
        argv, extra = argv[:cut], argv[cut + 1:]
    args = ap.parse_args(argv)

    if not hasattr(os, "posix_fadvise"):
        sys.exit("needs posix_fadvise(): Linux only")

    libc = ctypes.CDLL("libc.so.6", use_errno=True)
    files = tier_files(args.cache)
    total = sum(os.path.getsize(f) for f in files)

    print("index: %d tiers, %.2f GB in %s" % (len(files), total / 1e9, args.cache))
    if not args.no_baseline:
        evict(files)
        dt = sequential_read(files)
        print("sequential read of the whole index: %.2f s (%.2f GB/s) — the break-even"
              % (dt, total / dt / 1e9))
    print()

    header = "%-34s %-13s %9s %9s %9s %14s" % (
        "image", "outcome", "warm", "cold", "faults", "index touched")
    print(header)
    print("-" * len(header))

    with tempfile.TemporaryDirectory(prefix="index_cache_bench.") as outdir:
        for image in args.images:
            cmd = [args.solver, "-f", image, "-o", os.path.join(outdir, "out")]
            if args.database:
                cmd += ["-d", args.database]
            cmd += extra

            # The priming run pulls in the image and whatever tiers this solve
            # needs, so the timed warm runs measure the solve and not the I/O.
            _, _, code = solve(cmd)
            warm = min(solve(cmd)[0] for _ in range(args.warm_runs))

            cold = faults = touched = None
            if not args.no_cold:
                evict(files)
                left = resident_bytes(files, libc)
                if left > total * 0.01:
                    print("warning: %.2f GB stayed resident after eviction — is something "
                          "else holding the index mapped?" % (left / 1e9), file=sys.stderr)
                cold, faults, _ = solve(cmd)
                touched = resident_bytes(files, libc)

            print("%-34s %-13s %8.3fs %9s %9s %14s" % (
                os.path.basename(image)[:34],
                OUTCOMES.get(code, "exit %d" % code),
                warm,
                "%.2f s" % cold if cold is not None else "-",
                "%d" % faults if faults is not None else "-",
                "%.0f MB (%.1f%%)" % (touched / 1e6, 100 * touched / total)
                if touched is not None else "-"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
