# N.I.N.A. plugin

A plugin for [N.I.N.A.](https://nighttime-imaging.eu/) that makes it plate solve
against the resident quad index instead of ASTAP.

This directory is a separate module holding the whole integration: the resident
solver process, the local channel it is spoken to over, and the C# plugin.
Nothing under `../src` or `../include` refers to anything here — the module
consumes the solver libraries and adds to them. `astap_solve` and
`astap_index_solve` are what they were, and

    cmake -S . -B build -DASTAP_NINA_INTEGRATION=OFF

leaves this out of the build entirely.

## Why a plugin is needed at all, and what it does not do

N.I.N.A. 3.x has no extension point for plate solvers. `PlateSolverEnum` is a
closed set — ASTAP, ASTROMETRY_NET, PLATESOLVE2, PINPOINT, TSX_IMAGELINK — and
`PlateSolverFactory.GetSolver` switches on it. Plugins can supply star detection,
star annotation and autofocus behaviours; they cannot supply a solver. So this
plugin cannot appear in the plate solver dropdown, and the dropdown will keep
saying "ASTAP" while this is in use.

What N.I.N.A. *does* offer is that its ASTAP solver is a CLI solver: it launches
whatever executable `PlateSolveSettings.ASTAPLocation` points at, with

    -f "<image>" -z <bin> -fov <deg> -r <radius> -ra <hours> -spd <dec+90> -s <stars>

and reads the `.ini` written next to the image. `astap_index_server` takes that
option set and writes that `.ini`. The plugin therefore:

- ships `astap_index_server.exe` in its own plugin folder,
- writes `faster-astap.ini` next to it, holding the star database directory and
  the depth ladder — the settings ASTAP's option set has nowhere to carry,
- starts the server when N.I.N.A. starts, so the index is already in memory
  before the first solve of the night,
- points `ASTAPLocation` at the shipped executable, keeping the previous value so
  that disabling the plugin puts real ASTAP back,
- shows what is resident, and how long recent solves took.

Everything that goes through N.I.N.A.'s plate solving path is covered by this,
including the Framing Assistant's "solve this loaded image" prompt — which is
the case that benefits most, since an image loaded from disk has no pointing
hint and is therefore a blind solve.

## `astap_index_server`

The executable N.I.N.A. actually launches. One binary in two roles:

    astap_index_server -serve        hold the index and answer requests
    astap_index_server -f image.fits solve through it, taking ASTAP's options
    astap_index_server -status       what is resident, and recent solve times
    astap_index_server -stop         ask a running server to exit

Measured on a 6000×4000 DSLR frame against D80, all twelve tiers resident
(3.10 GB):

| | wall clock |
| --- | --- |
| `astap_index_solve`, one invocation | 4.43 s |
| `astap_index_server`, index not yet loaded | 3.46 s |
| `astap_index_server`, index resident | **0.13 s** |

The solve itself is 7 ms of that; the rest is starting a process and reading a
48 MB file. The `.ini` is byte for byte what `astap_index_solve` writes for the
same image.

Settings that ASTAP's option set cannot carry are read from `faster-astap.ini`
next to the executable — the star database directory, the depth ladder,
tolerance, cache location, endpoint, thread count, and whether a client may start
a server itself. When no server is running and none can be started, the client
solves in its own process rather than failing: slow, but a slow solve beats a
lost frame in the middle of a night.

## Layout

    nina-plugin/
      README.md                        this file
      CHANGELOG.md
      CMakeLists.txt                   builds astap_index_server
      server/index_server_main.cpp     the server, and the ASTAP-compatible client
      server/ipc.cpp, server/ipc.h     named pipe on Windows, Unix socket elsewhere
      server/version.rc                the version N.I.N.A. checks before it will call us
      FasterAstap.sln                  the plugin solution
      plugin/FasterAstapPlugin.cs      the manifest: settings, lifecycle, ASTAP path swap
      plugin/IndexServer.cs            starts and stops the server, and asks its status
      plugin/Options.xaml              the options page
      dist/FasterAstap.zip             the installable package

The two halves are built by different tools and share nothing but the wire
format in `IndexServer.cs` and `ipc.cpp`.

## Building

The server comes from the parent CMake build, into `build/nina-plugin/`. Build
it first — the plugin build copies the executable into its own output, and warns
rather than fails if it is not there yet.

The plugin needs the .NET 8 SDK (a newer SDK is fine, the target framework is
what matters: N.I.N.A. 3.x runs on .NET 8 and will not load a plugin built for a
later runtime):

    cmake --build build
    dotnet build nina-plugin/FasterAstap.sln -c Release

Output lands in `plugin/bin/x64/Release/`: the plugin assembly and the solver,
and nothing else — the N.I.N.A. assemblies are referenced for compilation only,
since the application supplies them at runtime.

## Installing

Four files: `FasterAstap.dll`, its `.pdb` and `.deps.json`, and
`astap_index_server.exe`. They belong together in one folder, because the plugin
points N.I.N.A. at the executable sitting beside it.

The reliable route is N.I.N.A.'s own installer: **Plugins ▸ install from
archive**, pointed at `dist/FasterAstap.zip`. That puts the folder wherever the
running version expects it, which since 3.0 is a version subfolder under
`%localappdata%\NINA\Plugins\` whose name follows the application version.
Copying the four files there by hand works too, if you would rather pick the
folder yourself.

Nothing is written to `Program Files` and no administrator rights are involved.

## Using it

Options ▸ Plugins ▸ Faster ASTAP:

1. Check the **database directory**. It is filled in from the ASTAP installation
   N.I.N.A. is already configured for, and only needs changing if the `.1476` or
   `.290` files live somewhere else.
2. Press **Start server**. The first run on a machine that has never built the
   index will take a while and then cache it; later runs read it back in
   seconds.
3. Tick **Use for plate solving**. That is the moment `ASTAPLocation` changes.
   Unticking it, or removing the plugin, puts the old path back.

The solver dropdown will still say ASTAP — N.I.N.A. has no way to list a plugin
there. The options page reports what is actually serving.

## Installing

Plugins live under `%localappdata%\NINA\Plugins\<version>\<name>\`, which is a
per-user directory: no administrator rights are involved, nothing is written to
`Program Files`, and removing the plugin removes the solver with it. The star
database is *not* copied there — it is gigabytes, and there is usually already
one next to the existing ASTAP installation, which is where the plugin looks
first.
