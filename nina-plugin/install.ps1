<#
.SYNOPSIS
    Installs, updates or removes the Faster ASTAP plugin in N.I.N.A.'s plugin folder.

.DESCRIPTION
    N.I.N.A. 3.x has no "install from archive" for plugins: a plugin is a folder
    under %localappdata%\NINA\Plugins\3.0.0\ and nothing else. This script puts
    the build output there, or takes it away again.

    Four files travel together — FasterAstap.dll, its .pdb and .deps.json, and
    astap_nina_solve.exe — because the plugin points N.I.N.A.'s ASTAP setting
    at the executable sitting beside its own assembly.

    Nothing is written outside %localappdata% and no administrator rights are
    involved. N.I.N.A. holds the plugin assembly open while it runs, so it has to
    be closed first; the script checks and says so rather than failing halfway
    through a copy.

    faster-astap.ini and faster-astap.log, which the plugin writes into its own
    folder, are left alone by an update — only the four shipped files are
    replaced.

.PARAMETER Source
    Folder holding the build output. Defaults to plugin\bin\x64\Release beside
    this script.

.PARAMETER Destination
    Plugin folder to install into. Defaults to
    %localappdata%\NINA\Plugins\3.0.0\FasterAstap.

    3.0.0 is not the application version. It is the plugin generation, and it
    stays 3.0.0 across all of N.I.N.A. 3.x.

.PARAMETER Uninstall
    Remove the plugin folder instead of installing into it.

.PARAMETER RemoveCache
    With -Uninstall, also delete the index cache under
    %localappdata%\faster-astap\cache. The cache is keyed to the star database
    and the quad tolerance rather than to this plugin, so a command line
    astap_index_solve on the same machine is sharing those files; deleting costs
    no data, only the minutes to rebuild.

.PARAMETER Force
    Kill a still-running astap_index_server — the resident server earlier
    versions shipped, which this installs over — if it ignores the request to
    stop. Never applied to N.I.N.A. itself.

.EXAMPLE
    .\install.ps1
    Install or update from the default build output.

.EXAMPLE
    .\install.ps1 -Uninstall -RemoveCache
    Remove the plugin folder and the index cache.
#>

# Below the help block, not above it: a #Requires ahead of comment-based help
# stops Windows PowerShell 5.1 associating the two, and -? then shows nothing
# but the syntax line.
#Requires -Version 5.1

[CmdletBinding(SupportsShouldProcess = $true)]
param(
    [string] $Source,
    [string] $Destination,
    [switch] $Uninstall,
    [switch] $RemoveCache,
    [switch] $Force
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

# The plugin assembly, the runtime it needs to be resolved against, and the
# solver it points N.I.N.A. at. The .pdb is optional: a Release build without
# symbols is still a working plugin, only a worse crash report.
$RequiredFiles = @('FasterAstap.dll', 'FasterAstap.deps.json', 'astap_nina_solve.exe')
$OptionalFiles = @('FasterAstap.pdb')

# Shipped by the versions that kept the index in a resident server, and removed
# by an update rather than left behind: nothing launches it any more, and a copy
# of it left running holds gigabytes.
$LegacyFiles = @('astap_index_server.exe')

function Resolve-SourceFolder {
    if ($Source) {
        if (-not (Test-Path -LiteralPath $Source -PathType Container)) {
            throw "Source folder not found: $Source"
        }
        return (Resolve-Path -LiteralPath $Source).Path
    }
    return (Join-Path $PSScriptRoot 'plugin\bin\x64\Release')
}

function Resolve-DestinationFolder {
    if ($Destination) { return $Destination }
    if (-not $env:LOCALAPPDATA) {
        throw 'LOCALAPPDATA is not set; pass -Destination with the plugin folder to use.'
    }
    return (Join-Path $env:LOCALAPPDATA 'NINA\Plugins\3.0.0\FasterAstap')
}

function Get-CacheFolder {
    if (-not $env:LOCALAPPDATA) { return $null }
    return (Join-Path $env:LOCALAPPDATA 'faster-astap\cache')
}

# N.I.N.A. loads the plugin into its own process and holds the assembly open for
# as long as it runs. Copying over it fails partway, which is a worse outcome
# than not starting, so this is a hard stop and -Force does not override it.
function Assert-NinaClosed {
    $nina = @(Get-Process -Name 'NINA' -ErrorAction SilentlyContinue)
    if ($nina.Count -gt 0) {
        $pids = ($nina | ForEach-Object { $_.Id }) -join ', '
        throw "N.I.N.A. is running (pid $pids) and holds the plugin assembly open. Close it and run this again."
    }
}

# The resident server of an earlier version, still running: from this session's
# N.I.N.A., or from one that crashed before its teardown ran. It keeps
# astap_index_server.exe locked and gigabytes resident, and nothing installed
# here will ever talk to it again. Ask it to stop the way that plugin would.
function Stop-LegacyServer {
    param([string[]] $ExeCandidates)

    if (@(Get-Process -Name 'astap_index_server' -ErrorAction SilentlyContinue).Count -eq 0) {
        return
    }

    Write-Host '  astap_index_server is running; asking it to stop'
    foreach ($exe in $ExeCandidates) {
        if ($exe -and (Test-Path -LiteralPath $exe -PathType Leaf)) {
            try { & $exe -stop 2>$null | Out-Null } catch { }
            break
        }
    }

    # It finishes the solve it is on before releasing the index, so give it a
    # few seconds rather than assuming an immediate exit.
    for ($i = 0; $i -lt 40; $i++) {
        if (@(Get-Process -Name 'astap_index_server' -ErrorAction SilentlyContinue).Count -eq 0) {
            Write-Host '  stopped'
            return
        }
        Start-Sleep -Milliseconds 250
    }

    $left = @(Get-Process -Name 'astap_index_server' -ErrorAction SilentlyContinue)
    if ($left.Count -eq 0) { return }

    if (-not $Force) {
        $pids = ($left | ForEach-Object { $_.Id }) -join ', '
        throw "astap_index_server (pid $pids) did not stop and holds its own executable open. Re-run with -Force to kill it."
    }

    Write-Host '  did not stop; killing it'
    $left | Stop-Process -Force
    Start-Sleep -Milliseconds 500
}

function Get-FileVersionString {
    param([string] $Path)
    try {
        $v = (Get-Item -LiteralPath $Path).VersionInfo.FileVersion
        if ($v) { return $v }
    } catch { }
    return ''
}

function Invoke-Install {
    $src = Resolve-SourceFolder
    $dst = Resolve-DestinationFolder

    $missing = @($RequiredFiles | Where-Object { -not (Test-Path -LiteralPath (Join-Path $src $_) -PathType Leaf) })
    if ($missing.Count -gt 0) {
        throw @"
Missing from ${src}:
  $($missing -join "`r`n  ")

Build both halves first, from the repository root:

  cmake --build build
  dotnet build nina-plugin/FasterAstap.sln -c Release

The solver comes from the CMake build and the plugin build copies it into its
own output, so the CMake build has to come first.
"@
    }

    Write-Host "Source      $src"
    Write-Host "Destination $dst"

    $existing = Test-Path -LiteralPath $dst -PathType Container
    if ($existing) {
        Write-Host 'Updating an existing installation'
    }

    Assert-NinaClosed
    Stop-LegacyServer -ExeCandidates @((Join-Path $dst 'astap_index_server.exe'))

    if (-not $existing) {
        if ($PSCmdlet.ShouldProcess($dst, 'Create plugin folder')) {
            New-Item -ItemType Directory -Path $dst -Force | Out-Null
        }
    }

    $copied = 0
    foreach ($name in ($RequiredFiles + $OptionalFiles)) {
        $from = Join-Path $src $name
        if (-not (Test-Path -LiteralPath $from -PathType Leaf)) { continue }
        if ($PSCmdlet.ShouldProcess((Join-Path $dst $name), 'Copy')) {
            Copy-Item -LiteralPath $from -Destination $dst -Force
            $copied++
        }
        Write-Host (("  {0,-28} {1}" -f $name, (Get-FileVersionString $from)).TrimEnd())
    }

    foreach ($name in $LegacyFiles) {
        $stale = Join-Path $dst $name
        if (-not (Test-Path -LiteralPath $stale -PathType Leaf)) { continue }
        if ($PSCmdlet.ShouldProcess($stale, 'Remove file from the previous design')) {
            Remove-Item -LiteralPath $stale -Force
            Write-Host ("  {0,-28} removed, no longer part of the plugin" -f $name)
        }
    }

    Write-Host ''
    if ($WhatIfPreference) {
        Write-Host 'Nothing was copied: -WhatIf.'
        return
    }
    Write-Host "$copied file(s) in place."

    $kept = @(@('faster-astap.ini', 'faster-astap.log') |
        Where-Object { Test-Path -LiteralPath (Join-Path $dst $_) -PathType Leaf })
    if ($kept.Count -gt 0) {
        Write-Host "Left as they were: $($kept -join ', ')"
    }

    Write-Host 'Start N.I.N.A. and look under Options > Plugins > Faster ASTAP.'
}

function Invoke-Uninstall {
    $dst = Resolve-DestinationFolder

    if (-not (Test-Path -LiteralPath $dst -PathType Container)) {
        Write-Host "Nothing installed at $dst"
    } else {
        Write-Host "Removing $dst"
        Assert-NinaClosed
        Stop-LegacyServer -ExeCandidates @((Join-Path $dst 'astap_index_server.exe'))
        if ($PSCmdlet.ShouldProcess($dst, 'Remove plugin folder')) {
            Remove-Item -LiteralPath $dst -Recurse -Force
            Write-Host '  removed'
            # The plugin restores N.I.N.A.'s ASTAP path on every exit, before
            # anything is checked or decided, so a N.I.N.A. that has already
            # been closed is not pointing at the executable this just deleted.
            Write-Host '  N.I.N.A.''s plate solver setting was restored when it last exited.'
        }
    }

    if ($RemoveCache) {
        $cache = Get-CacheFolder
        if ($cache -and (Test-Path -LiteralPath $cache -PathType Container)) {
            Write-Host "Removing $cache"
            if ($PSCmdlet.ShouldProcess($cache, 'Remove index cache')) {
                Remove-Item -LiteralPath $cache -Recurse -Force
                Write-Host '  removed'
            }
        } else {
            Write-Host 'No index cache to remove.'
        }
    } else {
        $cache = Get-CacheFolder
        if ($cache -and (Test-Path -LiteralPath $cache -PathType Container)) {
            Write-Host "Index cache left at $cache (pass -RemoveCache to delete it)."
        }
    }
}

try {
    if ($Uninstall) { Invoke-Uninstall } else { Invoke-Install }
} catch {
    Write-Host ''
    Write-Host $_.Exception.Message -ForegroundColor Red
    exit 1
}
