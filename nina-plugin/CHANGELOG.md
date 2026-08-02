# Changelog

## 1.0.0.0

First release.

- Keeps the whole-sky quad index in memory in a background process, so a plate
  solve costs milliseconds instead of the 1.6 seconds a cold solver spends
  reading the index back from disk.
- Points N.I.N.A.'s ASTAP path at the resident solver, remembering the previous
  value so that turning the plugin off restores it.
- Finds the star database beside the ASTAP installation already configured.
- Options page reporting what is resident, and recent solve times.
