# N.I.N.A. plugin

A plugin for [N.I.N.A.](https://nighttime-imaging.eu/) that makes it plate solve against the resident quad index instead
of ASTAP.

This directory is a separate module holding the whole integration: the resident solver process, the local channel it is
spoken to over, and the C# plugin. Nothing under `../src` or `../include` refers to anything here — the module consumes
the solver libraries and adds to them. `astap_solve` and
`astap_index_solve` are what they were, and

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

and reads the `.ini` written next to the image. `astap_index_server` takes that option set and writes that `.ini`. The
plugin therefore:

- ships `astap_index_server.exe` in its own plugin folder,
- writes `faster-astap.ini` next to it, holding the star database directory and the depth ladder — the settings ASTAP's
  option set has nowhere to carry,
- starts the server when N.I.N.A. starts, so the index is already in memory before the first solve of the night,
- points `ASTAPLocation` at the shipped executable, keeping the previous value so that disabling the plugin puts real
  ASTAP back,
- shows what is resident, and how long recent solves took.

Everything that goes through N.I.N.A.'s plate solving path is covered by this, including the Framing Assistant's "solve
this loaded image" prompt — which is the case that benefits most, since an image loaded from disk has no pointing hint
and is therefore a blind solve.

## `astap_index_server`

The executable N.I.N.A. actually launches. One binary in two roles:

    astap_index_server -serve        hold the index and answer requests
    astap_index_server -f image.fits solve through it, taking ASTAP's options
    astap_index_server -status       what is resident, and recent solve times
    astap_index_server -stop         ask a running server to exit

Measured on a 6000×4000 DSLR frame against D80, all twelve tiers resident (3.10 GB):

|                                            | wall clock |
|--------------------------------------------|------------|
| `astap_index_solve`, one invocation        | 4.43 s     |
| `astap_index_server`, index not yet loaded | 3.46 s     |
| `astap_index_server`, index resident       | **0.13 s** |

The solve itself is 7 ms of that; the rest is starting a process and reading a 48 MB file. The `.ini` is byte for byte
what `astap_index_solve` writes for the same image.

Settings that ASTAP's option set cannot carry are read from `faster-astap.ini`
next to the executable — the star database directory, the depth ladder, tolerance, cache location, endpoint, thread
count, and whether a client may start a server itself. When no server is running and none can be started, the client
solves in its own process rather than failing: slow, but a slow solve beats a lost frame in the middle of a night.

## Layout

    nina-plugin/
      README.md                        this file
      CHANGELOG.md
      CMakeLists.txt                   builds astap_index_server
      install.ps1                      puts the built plugin where N.I.N.A. looks
      install.cmd                      the same, for double-clicking
      uninstall.cmd                    install.ps1 -Uninstall, for double-clicking
      server/index_server_main.cpp     the server, and the ASTAP-compatible client
      server/ipc.cpp, server/ipc.h     named pipe on Windows, Unix socket elsewhere
      server/version.rc                the version N.I.N.A. checks before it will call us
      FasterAstap.sln                  the plugin solution
      plugin/FasterAstapPlugin.cs      the manifest: settings, lifecycle, ASTAP path swap
      plugin/IndexServer.cs            starts and stops the server, and asks its status
      plugin/Options.xaml              the options page

The two halves are built by different tools and share nothing but the wire format in `IndexServer.cs` and `ipc.cpp`.

## Building

The server comes from the parent CMake build, into `build/nina-plugin/`. Build it first — the plugin build copies the
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

Nothing is written to `Program Files` and no administrator rights are involved. N.I.N.A. holds the assembly open while
it runs, so it has to be closed first — the script checks, and says so, rather than failing halfway through a copy. It
also asks a leftover `astap_index_server` to stop before overwriting it, which is the state a crashed N.I.N.A. leaves
behind.

What it takes:

|                |                                                               |
|----------------|---------------------------------------------------------------|
| `-Source`      | where the build output is, if not `plugin\bin\x64\Release`    |
| `-Destination` | the plugin folder, if not the one above                       |
| `-Uninstall`   | remove the plugin folder instead                              |
| `-RemoveCache` | with `-Uninstall`, delete the index cache as well             |
| `-Force`       | kill an `astap_index_server` that ignored the request to stop |
| `-WhatIf`      | say what would be copied or deleted, and do neither           |

The star database is not copied there. It is gigabytes, and there is usually already one beside the existing ASTAP
installation, which is where the plugin looks first.

## Using it

Nothing, if the defaults suit. On the first start the plugin writes its settings, reads the index into memory, and then
points `ASTAPLocation` at its own solver.

The order there is the whole design. The path is taken over **after** the index is resident, never before: warming up
takes seconds when the ladder is cached and minutes the first time it has to be built, and a plate solve arriving in
that window would otherwise reach a solver that cannot answer yet. Until then plate solving stays with whatever was
configured, and the options page says so.

Warming up during N.I.N.A.'s own startup is the point: it costs nothing there, where paying the same seconds on the
first solve of the night makes something wait. That is what *read the index into memory when N.I.N.A. starts* does, and
it is separate from taking the ASTAP path over — a resident index is useful even when something else is doing the
solving, since any client finds a running server whoever launched it.

Options ▸ Plugins ▸ Faster ASTAP shows what is resident, what recent solves cost, and which path is in force. Worth
checking there:

- the **database directory**, filled in from the ASTAP installation N.I.N.A. was already configured for, and only wrong
  if the `.1476` or `.290` files live somewhere else;
- **depth tiers**, if the full ladder is more memory than you want to spend.

Unticking **Use for plate solving**, or removing the plugin, puts the old path back. The solver dropdown will still say
ASTAP — N.I.N.A. has no way to list a plugin there.

## Shutting down, and shutting down badly

The server is a second process holding gigabytes, so it must never outlive the application it was started for. Two
mechanisms, because one of them is not enough:

- N.I.N.A. closing normally runs the plugin's teardown, which asks the server to stop. It finishes the solve it is on,
  releases the index, and exits.
- N.I.N.A. crashing, or being killed, never reaches that. So the server is also given N.I.N.A.'s process id with
  `-parent` and watches it: when that process ends, by whatever means, the server logs it and exits. Verified by killing
  the parent outright — the memory comes back within a second or two.

`-idle-exit` is a third line, off by default: a server that has been asked to solve nothing for N minutes stops on its
own. It exists for a machine where something started a server outside of N.I.N.A. entirely.

## Removing it

N.I.N.A. offers a plugin no way to know it is being uninstalled. It moves the plugin folder aside and deletes it on the
next start, so by the time the removal is complete there is nothing of this plugin left running to tidy up after itself.
Two things follow, and they are worth separating because one is guaranteed and the other is not.

**The ASTAP path is always restored.** Not on uninstall — on *every* exit, before anything is checked or decided. The
path points here only while the plugin is loaded and the index is resident, and it is re-adopted on the next start. So a
plugin folder that disappears between one run and the next cannot leave
`ASTAPLocation` standing against a missing executable, however it disappeared. This is the failure that would actually
cost you a night, and it does not depend on detecting anything.

**The index cache is removed on a best-effort basis.** The cache is the one thing that outlives the plugin folder:
several gigabytes under
`%localappdata%\faster-astap\cache`, put there deliberately so that it survives plugin updates. Deleting it needs the
plugin to know it is being uninstalled, and what it actually knows is narrower — at teardown it looks for
`astap_index_server.exe` beside its own assembly, and takes its absence to mean the folder has already been moved to
`PluginDeletion` while the session ran. That holds when you uninstall from a running N.I.N.A. and then close it. It does
not hold if you delete the plugin folder while N.I.N.A. is closed, because then none of this code ever runs again. In
that case the cache is simply left where it is.

Leaving a cache behind costs disk space and nothing else, so best effort is an acceptable answer here in a way it would
not be for the solver path.

**Restore ASTAP and remove everything** on the options page is the reliable route:
it restores the path, stops the server, and deletes the cache and the generated files there and then, with no inference
involved. It cannot delete the plugin itself — N.I.N.A. does that afterwards.

With N.I.N.A. already closed, the reliable route is the script, which has the same nothing-inferred property:

    powershell -ExecutionPolicy Bypass -File nina-plugin\install.ps1 -Uninstall -RemoveCache

or `uninstall.cmd`, which is that without the `-RemoveCache`: double-clicking it takes the plugin folder and leaves the
cache, and says where it left it. Pass
`-RemoveCache` to take that too.

Before leaving *delete the index cache when the plugin is uninstalled* ticked, note that the cache is keyed to the star
database and the quad tolerance, not to this plugin: a command line `astap_index_solve` on the same machine is sharing
those very files. Deleting costs no data, only the minutes to rebuild.
