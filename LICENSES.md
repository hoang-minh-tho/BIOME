# Third-party components

The Windows build freezes Cutadapt and its Python dependencies into
`tools/cutadapt.exe`. The build environment is pinned in
`requirements-build.txt`.

Before distributing a build, review and retain the license and notice files for
Cutadapt, its transitive dependencies, PyInstaller, Python, and any other
bundled components. This project does not include R or R/Bioconductor packages;
those are installed separately into the user's R library from CRAN and
Bioconductor and remain subject to their respective licenses.
