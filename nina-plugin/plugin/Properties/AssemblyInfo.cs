using System.Reflection;
using System.Runtime.InteropServices;

// The unique identifier of this plugin. N.I.N.A. also uses it to key the
// plugin's profile-specific settings, so it must not change between releases.
[assembly: Guid("3138b765-ce6e-46f3-9036-68d229bcec15")]

[assembly: AssemblyVersion("1.1.0.0")]
[assembly: AssemblyFileVersion("1.1.0.0")]

// The name. The options page finds its user interface by a DataTemplate keyed
// "<AssemblyTitle>_Options", so this string and the key in Options.xaml have to
// agree exactly or the options tab comes up empty.
[assembly: AssemblyTitle("Faster ASTAP")]
[assembly: AssemblyDescription("Plate solves against a pre-built whole-sky quad index, in place of ASTAP.")]

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
[assembly: AssemblyMetadata("LongDescription", @"Plate solving that does not search the sky.

The index solver looks an image's star quads up in a pre-built whole-sky table instead of walking a spiral outwards from where the telescope thinks it is pointed. That makes a solve indifferent to the pointing hint, and makes a blind solve cost the same as any other.

This plugin points N.I.N.A.'s ASTAP setting at that solver, so everything that plate solves in N.I.N.A. goes through it — including the Framing Assistant's prompt when you load an image from a file, which is a blind solve and where the difference is largest.

The table is gigabytes, and it is read through a memory mapping rather than copied: the parts one solve touches stay in the operating system's cache for the solves after it, so there is no background process here and nothing to warm up. It is built once per star database, which takes minutes, and the options page has a button for doing that before a session rather than during one.

Requires a star database (the .1476 or .290 files ASTAP uses); an existing ASTAP installation already has one, and the plugin will find it.

N.I.N.A. has no extension point for plate solvers, so this appears in the solver dropdown as ASTAP. The plugin's options page shows which solver is actually in use.")]

[assembly: ComVisible(false)]
