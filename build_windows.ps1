$ErrorActionPreference = "Stop"

$repo = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$buildEnv = Join-Path $PSScriptRoot ".build-venv"
$python = Join-Path $buildEnv "Scripts\python.exe"
$requirements = Join-Path $PSScriptRoot "requirements-build.txt"
$buildDir = Join-Path $PSScriptRoot "build"
$distDir = Join-Path $repo "dist"
$stagingDist = Join-Path $buildDir "staging-dist"
$appDist = Join-Path $distDir "BIOME"
$cutadaptDist = Join-Path $buildDir "cutadapt-dist"
$cutadaptExe = Join-Path $cutadaptDist "cutadapt.exe"

$buildPythonVersion = ""
if (Test-Path $python) {
    $buildPythonVersion = & $python -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')"
}

if ($buildPythonVersion.Trim() -ne "3.13") {
    $basePython = Join-Path $env:LOCALAPPDATA "Programs\Python\Python313\python.exe"
    if (-not (Test-Path $basePython)) {
        $basePython = ""
        $pythonCommand = Get-Command python.exe -ErrorAction SilentlyContinue
        if ($pythonCommand -and $pythonCommand.Source -notlike "*\WindowsApps\*") {
            $basePython = $pythonCommand.Source
        }
    }
    if (-not $basePython) {
        throw "Python 3.13 is required on the build machine because Cutadapt publishes a Windows wheel for that runtime. Python is bundled in the output app and is not required by end users."
    }
    $pythonVersion = & $basePython -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')"
    if ($LASTEXITCODE -ne 0 -or $pythonVersion.Trim() -ne "3.13") {
        throw "The packaging build requires Python 3.13; found $($pythonVersion.Trim())."
    }
    if (Test-Path $buildEnv) {
        Remove-Item -LiteralPath $buildEnv -Recurse -Force
    }
    & $basePython -m venv $buildEnv
    if ($LASTEXITCODE -ne 0) { throw "Could not create the packaging virtual environment." }
}

& $python -m pip install -r $requirements
if ($LASTEXITCODE -ne 0) { throw "Could not install pinned packaging dependencies." }

New-Item -ItemType Directory -Force -Path $buildDir, $cutadaptDist, $distDir | Out-Null
Push-Location $repo
try {
    & $python -m PyInstaller --noconfirm --clean --onefile --console `
        --name cutadapt `
        --distpath $cutadaptDist `
        --workpath (Join-Path $buildDir "cutadapt-work") `
        --collect-all cutadapt `
        --collect-all dnaio `
        --collect-all xopen `
        --collect-all isal `
        (Join-Path $PSScriptRoot "cutadapt_entry.py")
    if ($LASTEXITCODE -ne 0) { throw "Could not freeze the Cutadapt executable." }

    & $python -m PyInstaller --noconfirm --clean --onedir --windowed `
        --name BIOME `
        --distpath $stagingDist `
        --workpath (Join-Path $buildDir "app-work") `
        --icon=assets/biome.ico `
        --add-data=run_pipeline.R:. `
        --add-data=pipeline:pipeline `
        --add-data=packaging/build/cutadapt-dist/cutadapt.exe:tools `
        --add-data=packaging/LICENSES.md:licenses `
        --add-data=assets/biome.ico:assets `
        --add-data=assets/biome-rounded.png:assets `
        (Join-Path $repo "run_pipeline.py")
    if ($LASTEXITCODE -ne 0) { throw "Could not freeze the desktop application." }
}
finally {
    Pop-Location
}

$stagedAppDist = Join-Path $stagingDist "BIOME"

$sitePackages = Join-Path $buildEnv "Lib\site-packages"
$thirdPartyLicenses = Join-Path $stagedAppDist "_internal\licenses\third-party"
New-Item -ItemType Directory -Force -Path $thirdPartyLicenses | Out-Null
Get-ChildItem -Path $sitePackages -Directory -Filter "*.dist-info" | ForEach-Object {
    $licenseDir = Join-Path $_.FullName "licenses"
    $packageLicenseDir = Join-Path $thirdPartyLicenses $_.Name
    if (Test-Path $licenseDir) {
        Get-ChildItem -Path $licenseDir -File -Recurse | ForEach-Object {
            $relativePath = $_.FullName.Substring($licenseDir.Length).TrimStart("\")
            $destination = Join-Path $packageLicenseDir $relativePath
            New-Item -ItemType Directory -Force -Path (Split-Path $destination) | Out-Null
            Copy-Item -LiteralPath $_.FullName -Destination $destination -Force
        }
    }
    Get-ChildItem -Path $_.FullName -File | Where-Object {
        $_.Name -match "^(LICENSE|COPYING|NOTICE)"
    } | ForEach-Object {
        New-Item -ItemType Directory -Force -Path $packageLicenseDir | Out-Null
        Copy-Item -LiteralPath $_.FullName -Destination $packageLicenseDir -Force
    }
}

$null = & robocopy $stagedAppDist $appDist /MIR /R:2 /W:1
if ($LASTEXITCODE -ge 8) {
    throw "Could not publish the rebuilt application to $appDist (robocopy exit code $LASTEXITCODE). Close any running BIOME instance and retry."
}

$zipPath = Join-Path $distDir "BIOME-portable.zip"
try {
    if (Test-Path $zipPath) { Remove-Item -LiteralPath $zipPath -Force }
    Compress-Archive -Path (Join-Path $appDist "*") -DestinationPath $zipPath -CompressionLevel Optimal -ErrorAction Stop
    Write-Output "Portable package created: $zipPath"
}
catch {
    Write-Warning "Could not create the optional portable ZIP: $($_.Exception.Message)"
}

$iscc = Get-Command ISCC.exe -ErrorAction SilentlyContinue
$isccPath = if ($iscc) { $iscc.Source } else { "" }
if (-not $isccPath) {
    foreach ($candidate in @(
        (Join-Path ${env:ProgramFiles(x86)} "Inno Setup 7\ISCC.exe"),
        (Join-Path $env:ProgramFiles "Inno Setup 7\ISCC.exe"),
        (Join-Path ${env:ProgramFiles(x86)} "Inno Setup 6\ISCC.exe"),
        (Join-Path $env:ProgramFiles "Inno Setup 6\ISCC.exe")
    )) {
        if (Test-Path $candidate) {
            $isccPath = $candidate
            break
        }
    }
}
if ($isccPath) {
    & $isccPath (Join-Path $PSScriptRoot "BIOME.iss")
    if ($LASTEXITCODE -ne 0) { throw "Inno Setup failed to create Setup.exe." }
    Write-Output "Installer created in $distDir"
}
else {
    Write-Warning "Inno Setup was not found. The portable ZIP is ready; install Inno Setup and rerun this script to create BIOME-Setup.exe."
}
