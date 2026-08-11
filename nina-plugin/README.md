# N.I.N.A. plugin

A plugin for [N.I.N.A.](https://nighttime-imaging.eu/) that makes it plate solve against the whole-sky quad index instead
of ASTAP.

This directory is a separate module holding the whole integration: the solver N.I.N.A. launches, and the C# plugin that
puts it in ASTAP's place. Nothing under `../src` or `../include` refers to anything here — the module consumes the solver
libraries and adds to them. `astap_solve` and `astap_index_solve` are what they were, and

    cmake -S . -B build -DASTAP_NINA_INTEGRATION=OFF

leaves this out of the build entirely.

## Why a plugin is needed at all, and what it does not do

N.I.N.A. 3.x has no extension point for plate solvers. `PlateSolverEnum` is a closed set — ASTAP, ASTROMETRY_NET,
PLATESOLVE2, PINPOINT, TSX_IMAGELINK — and
`PlateSolverFactory.GetSolver` switches on it. Plugins can supply star detection, star annotation and autofocus
behaviours; they cannot supply a solver. So this plugin cannot appear in the plate solver dropdown, and the dropdown
will keep saying "ASTAP" while this is in use.

What N.I.N.A. *does* offer is that its ASTAP solver is a CLI solver: it launches whatever executable
`PlateSolveSettings.ASTAPLocation` points at, with

    -f "<image>" -z <bin> -fov <deg> -r <radius> -ra <hours> -spd <dec+90> -s <stars>

and reads the `.ini` written next to the image. `astap_nina_solve` takes that option set and writes that `.ini`. The
plugin therefore:

- ships `astap_nina_solve.exe` in its own plugin folder,
- writes `faster-astap.ini` next to it, holding the star database directory and the depth ladder — the settings ASTAP's
  option set has nowhere to carry,
- points `ASTAPLocation` at the shipped executable, keeping the previous value so that disabling the plugin puts real
  ASTAP back,
- shows what is cached, and what the last solves did.

Everything that goes through N.I.N.A.'s plate solving path is covered by this, including the Framing Assistant's "solve
this loaded image" prompt — which is the case that benefits most, since an image loaded from disk has no pointing hint
and is therefore a blind solve.

Nothing of the plugin runs while a solve does. N.I.N.A. starts the executable, the executable writes the `.ini` and
exits.

## `astap_nina_solve`

The executable N.I.N.A. actually launches. It is `astap_index_solve` with the two differences that come of being
launched by a program rather than typed by a person: it accepts and ignores the options ASTAP takes that an index solve
has no use for (`-r`, `-ra`, `-spd`), and it reads the settings that option set cannot carry from `faster-astap.ini`
beside it. The `.ini` it writes is byte for byte what `astap_index_solve` writes for the same image.

    astap_nina_solve -f image.fits     solve, taking ASTAP's options
    astap_nina_solve -prepare          build or read the index and exit
    astap_nina_solve -h                the whole option set

Measured here on a 1024×1024 corpus frame against D80 with the full twelve-tier ladder (3.10 GB) cached:

|                                            | wall clock |
|--------------------------------------------|------------|
| `astap_index_solve`, one invocation        | 0.61 s     |
| `astap_nina_solve`, one invocation         | 0.60 s     |
| `astap_nina_solve -prepare` (no image)     | 0.03 s     |

The two are the same because they run the same code. The third line is the point: mapping the whole 3.10 GB ladder and
reporting it ready costs 30 ms in a process that has just started, so a solve pays essentially nothing for having the
index — the rest of that 0.6 s is reading the image, detecting stars and sweeping the tiers.

### Why there is no longer a server

Earlier versions shipped `astap_index_server`: a process that held the index in memory and answered solve requests over
a named pipe, because reading a 2.7 GB ladder into a buffer took 1.65 seconds and the solve after it took 5 ms — a cost
worth paying once for a batch and worth avoiding entirely for an imaging application, which solves one frame at a time,
minutes apart, in a fresh process every time.

The index cache is now memory mapped rather than read into a buffer (see `MappedFile`, and the reasons in
`include/astap/mapped_file.h`), and that removes the cost the server existed to amortise. Mapping a file is not reading
it: the pages a solve touches are faulted in from the operating system's page cache, they stay there when the process
exits, and the next process maps the same file and finds them resident. Nothing has to hold them — not holding them is
what a page cache is for. What is left for a server to save is the process start, and that is the 30 ms above.

So the server is gone, and with it the named pipe, the wire format, the parent-process watchdog, the idle timeout, and
every failure mode that came of a second process holding gigabytes: a stranded server after a crash, a solve arriving
while the server was still warming up, a client that could not reach one. A plate solve is now one program launch that
either writes an `.ini` or does not.

Two costs are real and were worth naming rather than hiding:

- **The first solve after the machine boots** faults its pages in off the disk, and is slower than the ones after it by
  whatever that disk charges. Any later solve finds them cached.
- **A ladder that has never been built** costs minutes, once per star database. That used to happen while N.I.N.A.
  started up; now it would happen on the first frame that needs it, which is a worse place for it. **Build the index
  now** on the options page runs `-prepare` and pays it at a chosen moment instead.

Settings that ASTAP's option set cannot carry are read from `faster-astap.ini` next to the executable — the star
database directory, the depth ladder, tolerance, cache location, thread count, and the log file. The plugin writes that
file; a person can edit it, and a rebuild of the plugin does not overwrite it.

## Layout

    nina-plugin/
      README.md                        this file
      CHANGELOG.md
      CMakeLists.txt                   builds astap_nina_solve
      install.ps1                      puts the built plugin where N.I.N.A. looks
      install.cmd                      the same, for double-clicking
      uninstall.cmd                    install.ps1 -Uninstall, for double-clicking
      solver/nina_solve_main.cpp       the ASTAP-compatible front end
      solver/version.rc                the version N.I.N.A. checks before it will call us
      FasterAstap.sln                  the plugin solution
      plugin/FasterAstapPlugin.cs      the manifest: settings, lifecycle, ASTAP path swap
      plugin/Solver.cs                 the settings file, the log, and the cache on disk
      plugin/Options.xaml              the options page

The two halves are built by different tools and share nothing but the contents of `faster-astap.ini`.

## Building

The solver comes from the parent CMake build, into `build/nina-plugin/`. Build it first — the plugin build copies the
executable into its own output, and warns rather than fails if it is not there yet.

The plugin needs the .NET 8 SDK (a newer SDK is fine, the target framework is what matters: N.I.N.A. 3.x runs on .NET 8
and will not load a plugin built for a later runtime):

    cmake --build build
    dotnet build nina-plugin/FasterAstap.sln -c Release

Output lands in `plugin/bin/x64/Release/`: the plugin assembly and the solver, and nothing else — the N.I.N.A.
assemblies are referenced for compilation only, since the application supplies them at runtime.

## Installing

Run the install script that copies the faster-astap files:

    powershell -ExecutionPolicy Bypass -File nina-plugin\install.ps1

or double-click `install.cmd`, which is the same thing with the execution policy already dealt with — as is
`uninstall.cmd` for the other direction. It takes the build output from `plugin\bin\x64\Release` and puts it in

    %localappdata%\NINA\Plugins\3.0.0\FasterAstap\

`3.0.0` is not the application version. It is the plugin generation, and it stays
`3.0.0` across all of N.I.N.A. 3.x. Anything placed directly in `Plugins\` is treated as a leftover from before that
folder existed and is migrated *into*
`3.0.0` on the next start, so a folder named for the application version ends up one level too deep and is never seen.

Run it again after a rebuild and it updates in place. Only the four files are replaced; `faster-astap.ini` and
`faster-astap.log`, which the plugin writes into that same folder, are left as they are, so settings survive an update.
Installing over a version that shipped `astap_index_server.exe` stops a copy of it that is still running and removes the
file, since nothing will launch it again.

Nothing is written to `Program Files` and no administrator rights are involved. N.I.N.A. holds the assembly open while
it runs, so it has to be closed first — the script checks, and says so, rather than failing halfway through a copy.

What it takes:

|                |                                                                    |
|----------------|--------------------------------------------------------------------|
| `-Source`      | where the build output is, if not `plugin\bin\x64\Release`          |
| `-Destination` | the plugin folder, if not the one above                            |
| `-Uninstall`   | remove the plugin folder instead                                   |
| `-RemoveCache` | with `-Uninstall`, delete the index cache as well                  |
| `-Force`       | kill an old `astap_index_server` that ignored the request to stop  |
| `-WhatIf`      | say what would be copied or deleted, and do neither                |

The star database is not copied there. It is gigabytes, and there is usually already one beside the existing ASTAP
installation, which is where the plugin looks first.

## Using it

Nothing, if the defaults suit. On the first start the plugin writes its settings and points `ASTAPLocation` at its own
solver.

Options ▸ Plugins ▸ Faster ASTAP shows which path is in force, what is cached, and the last few lines the solver wrote —
which is the only account of a solve there is, since N.I.N.A. launches it with no console anybody can see. Worth
checking there:

- the **database directory**, filled in from the ASTAP installation N.I.N.A. was already configured for, and only wrong
  if the `.1476` or `.290` files live somewhere else;
- **depth tiers**, if the full ladder is more disk than you want to spend;
- **Build the index now**, before the first session on a machine that has never built one — that is the one operation
  here that takes minutes.

Unticking **Use for plate solving**, or removing the plugin, puts the old path back. The solver dropdown will still say
ASTAP — N.I.N.A. has no way to list a plugin there.

## Removing it

N.I.N.A. offers a plugin no way to know it is being uninstalled. It moves the plugin folder aside and deletes it on the
next start, so by the time the removal is complete there is nothing of this plugin left running to tidy up after itself.
Two things follow, and they are worth separating because one is guaranteed and the other is not.

**The ASTAP path is always restored.** Not on uninstall — on *every* exit, before anything is checked or decided. The
path points here only while the plugin is loaded, and it is re-adopted on the next start. So a plugin folder that
disappears between one run and the next cannot leave
`ASTAPLocation` standing against a missing executable, however it disappeared. This is the failure that would actually
cost you a night, and it does not depend on detecting anything.

**The index cache is removed on a best-effort basis.** The cache is the one thing that outlives the plugin folder:
several gigabytes under
`%localappdata%\faster-astap\cache`, put there deliberately so that it survives plugin updates. Deleting it needs the
plugin to know it is being uninstalled, and what it actually knows is narrower — at teardown it looks for
`astap_nina_solve.exe` beside its own assembly, and takes its absence to mean the folder has already been moved to
`PluginDeletion` while the session ran. That holds when you uninstall from a running N.I.N.A. and then close it. It does
not hold if you delete the plugin folder while N.I.N.A. is closed, because then none of this code ever runs again. In
that case the cache is simply left where it is.

Leaving a cache behind costs disk space and nothing else, so best effort is an acceptable answer here in a way it would
not be for the solver path.

**Restore ASTAP and remove everything** on the options page is the reliable route:
it restores the path and deletes the cache and the generated files there and then, with no inference involved. It cannot
delete the plugin itself — N.I.N.A. does that afterwards.

With N.I.N.A. already closed, the reliable route is the script, which has the same nothing-inferred property:

    powershell -ExecutionPolicy Bypass -File nina-plugin\install.ps1 -Uninstall -RemoveCache

or `uninstall.cmd`, which is that without the `-RemoveCache`: double-clicking it takes the plugin folder and leaves the
cache, and says where it left it. Pass
`-RemoveCache` to take that too.

Before leaving *delete the index cache when the plugin is uninstalled* ticked, note that the cache is keyed to the star
database and the quad tolerance, not to this plugin: a command line `astap_index_solve` on the same machine is sharing
those very files. Deleting costs no data, only the minutes to rebuild.
