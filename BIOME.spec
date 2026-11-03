# -*- mode: python ; coding: utf-8 -*-


a = Analysis(
    ['C:/Users/Grace_Hoang/Desktop/Tool_project/databio_16S_python_R/run_pipeline.py'],
    pathex=[],
    binaries=[],
    datas=[('run_pipeline.R', '.'), ('pipeline', 'pipeline'), ('packaging/build/cutadapt-dist/cutadapt.exe', 'tools'), ('packaging/LICENSES.md', 'licenses'), ('assets/biome.ico', 'assets'), ('assets/biome-rounded.png', 'assets')],
    hiddenimports=[],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name='BIOME',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon=['assets/biome.ico'],
)
coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name='BIOME',
)
