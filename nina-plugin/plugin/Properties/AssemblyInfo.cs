using System.Reflection;
using System.Runtime.InteropServices;

// The unique identifier of this plugin. N.I.N.A. also uses it to key the
// plugin's profile-specific settings, so it must not change between releases.
[assembly: Guid("3138b765-ce6e-46f3-9036-68d229bcec15")]

[assembly: AssemblyVersion("1.0.0.0")]
[assembly: AssemblyFileVersion("1.0.0.0")]

// The name. The options page finds its user interface by a DataTemplate keyed
// "<AssemblyTitle>_Options", so this string and the key in Options.xaml have to
// agree exactly or the options tab comes up empty.
[assembly: AssemblyTitle("Faster ASTAP")]
[assembly: AssemblyDescription("Plate solves against a whole-sky quad index held in memory, in place of ASTAP.")]

[assembly: AssemblyCompany("Vlad Pop")]
[assembly: AssemblyProduct("Faster ASTAP")]
[assembly: AssemblyCopyright("Copyright © 2026 Vlad Pop. Original ASTAP algorithm © 2018-2026 Han Kleijn.")]

// The guideline is to state the version of the NINA.Plugin package this was
// built against, and claiming less than that is not a favour to anyone: it lets
// an older N.I.N.A. load a plugin compiled against interfaces it does not have,
// which fails at the point of use rather than at the point of loading.
[assembly: AssemblyMetadata("MinimumApplicationVersion", "3.2.0.9001")]

[assembly: AssemblyMetadata("License", "MPL-2.0")]
[assembly: AssemblyMetadata("LicenseURL", "https://www.mozilla.org/en-US/MPL/2.0/")]
[assembly: AssemblyMetadata("Repository", "https://github.com/Pop-Vlad/faster-astap")]

[assembly: AssemblyMetadata("Homepage", "https://github.com/Pop-Vlad/faster-astap")]
[assembly: AssemblyMetadata("Tags", "Plate Solving,ASTAP,Astrometry")]
[assembly: AssemblyMetadata("ChangelogURL", "https://github.com/Pop-Vlad/faster-astap/blob/main/nina-plugin/CHANGELOG.md")]
[assembly: AssemblyMetadata("FeaturedImageURL", "")]
[assembly: AssemblyMetadata("ScreenshotURL", "")]
[assembly: AssemblyMetadata("AltScreenshotURL", "")]
[assembly: AssemblyMetadata("LongDescription", @"Plate solving without the wait for the index.

The index solver looks an image's star quads up in a pre-built whole-sky table instead of searching the sky, which makes a solve a few milliseconds and makes it indifferent to where the telescope was pointed. The catch is that the table is gigabytes, and a solver launched afresh for every frame spends over a second reading it back before it can start.

This plugin keeps that table in memory in a small background process, and points N.I.N.A.'s ASTAP setting at a program that solves through it. A warm solve, end to end and including reading the image off disk, takes about a tenth of a second.

Everything that plate solves in N.I.N.A. goes through it, including the Framing Assistant's prompt when you load an image from a file — which is a blind solve, and where the difference is largest.

Requires a star database (the .1476 or .290 files ASTAP uses); an existing ASTAP installation already has one, and the plugin will find it.

N.I.N.A. has no extension point for plate solvers, so this appears in the solver dropdown as ASTAP. The plugin's options page shows what is actually serving.")]

[assembly: ComVisible(false)]
