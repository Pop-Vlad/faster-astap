# Changelog

## 1.1.0.0

- The resident index server is gone. The index cache is memory mapped, so the rungs a solve touches stay in the
  operating system's page cache and the next solve finds them there — a process holding them saved 30 ms and cost a
  second process holding gigabytes. N.I.N.A. now launches `astap_nina_solve.exe` per solve and nothing else runs.
- **Build the index now** on the options page builds a ladder that has never been built, which is the one operation
  here that takes minutes, at a moment of your choosing rather than on the first frame of a session.
- The options page reports what the solver last wrote, in place of what a server was holding.
- Installing over the previous version stops a leftover `astap_index_server.exe`, removes it, and takes an ASTAP path
  that was pointed at it back to the real ASTAP.

## 1.0.0.0

First release.

- Keeps the whole-sky quad index in memory in a background process, so a plate solve costs milliseconds instead of the
  1.6 seconds a cold solver spends reading the index back from disk.
- Points N.I.N.A.'s ASTAP path at the resident solver, remembering the previous value so that turning the plugin off
  restores it.
- Finds the star database beside the ASTAP installation already configured.
- Options page reporting what is resident, and recent solve times.
