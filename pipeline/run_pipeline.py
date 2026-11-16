"""Select paired-end amplicon inputs and run the supplied R pipeline.

Python 3.9+ with Tkinter. No pip packages are required.
The analysis runs in one Rscript process, with stdout/stderr saved to disk.
"""

import csv
import ctypes
import errno
import json
import math
import os
from pathlib import Path
import queue
import re
import shutil
import subprocess
import sys
import tempfile
import threading
import traceback
from datetime import datetime
from functools import partial

LOCAL_APP_DATA = Path(os.environ.get("LOCALAPPDATA") or Path.home() / "AppData/Local")
SETTINGS_DIR = LOCAL_APP_DATA / "BIOME"
LEGACY_APP_SETTINGS = LOCAL_APP_DATA / "Databio16S" / "settings.json"


def _write_startup_log(details):
    log_path = SETTINGS_DIR / "startup.log"
    try:
        SETTINGS_DIR.mkdir(parents=True, exist_ok=True)
        with log_path.open("a", encoding="utf-8") as log:
            log.write(f"\n--- {datetime.now().isoformat(timespec='seconds')} ---\n")
            log.write(details)
    except OSError as exc:
        return None, exc
    return log_path, None


def _show_native_startup_error(message):
    if os.name != "nt":
        print(message, file=sys.stderr)
        return
    message_box = ctypes.windll.user32.MessageBoxW
    message_box.argtypes = (
        ctypes.c_void_p,
        ctypes.c_wchar_p,
        ctypes.c_wchar_p,
        ctypes.c_uint,
    )
    message_box.restype = ctypes.c_int
    message_box(None, message, "BIOME could not start", 0x10)


def _uncaught_startup_exception(exc_type, exc_value, exc_traceback):
    details = "".join(traceback.format_exception(exc_type, exc_value, exc_traceback))
    log_path, log_error = _write_startup_log(details)
    message = (
        "BIOME encountered an unexpected error and could not start.\n\n"
        f"{exc_type.__name__}: {exc_value}"
    )
    if log_path is not None:
        message += f"\n\nDetails were saved to:\n{log_path}"
    elif log_error is not None:
        message += f"\n\nCould not write startup log: {log_error}\n\n{details}"
    try:
        _show_native_startup_error(message)
    except (AttributeError, OSError, ctypes.ArgumentError) as display_error:
        print(f"{message}\nCould not display startup error: {display_error}", file=sys.stderr)


def _uncaught_thread_exception(args):
    _uncaught_startup_exception(args.exc_type, args.exc_value, args.exc_traceback)


sys.excepthook = _uncaught_startup_exception
threading.excepthook = _uncaught_thread_exception

import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext, ttk
import zipfile

from app_setup import SetupWizard, find_cutadapt, find_rscript

BASE = Path(getattr(sys, "_MEIPASS", Path(__file__).resolve().parent))
MARKER_PROFILES_PATH = BASE / "pipeline" / "marker_profiles.json"
SETTINGS = SETTINGS_DIR / "settings.json"
APP_TEMP_ROOT = SETTINGS_DIR / "temp"
LEGACY_SETTINGS = Path(__file__).resolve().parent / "last_selection.json"
PAIR_PATTERN = re.compile(r"(.+)(?:_R([12])|_([12]))(?:_001)?\.(?:fastq|fq)(?:\.gz)?$")
ILLUMINA_LANE_SUFFIX = re.compile(r"(?:_S\d+)?_L\d+$")
FASTQ_SUFFIXES = (".fastq", ".fastq.gz", ".fq", ".fq.gz")
MESSAGES = {
        "file_missing": "File does not exist or is empty: {path}",
        "fastq_name": "FASTQ names must end with _R1/_R2 or _1/_2: {name}",
        "duplicate_sample": "Duplicate sample {sample}, R{direction}. Select one pair per sample.",
        "select_fastq": "Select paired-end FASTQ files to analyze.",
        "missing_mates": "Missing R1/R2 for samples: {samples}",
        "pair_order": "FASTQ pair order differs between folders. Keep R1/R2 in the same folder.",
        "metadata_csv": "Metadata must be a CSV file with a SampleID column, a header, and data.",
        "duplicate_columns": "Metadata contains duplicate column names.",
        "metadata_row": "A metadata row has a different number of fields than the CSV header.",
        "sample_ids": "The SampleID column must contain unique, non-empty values.",
        "window_title": "BIOME — Bioinformatic Interface for Optimized Multi-marker Ecology",
        "heading": "Analyze your amplicon sequencing data",
        "clear_temp_cache": "🧹 Clear Temporary Cache",
        "temp_cache_description": "Cleans up temporary files left behind by completed or cancelled runs. Does not delete your analysis results, FASTQ files, or any .rds files.",
        "temp_cache_cleared": "Temporary cache cleared! (Freed {size})",
        "temp_cache_error": "Could not clear the temporary cache: {error}",
        "temp_cache_running": "Stop the pipeline before clearing temporary files.",
        "settings_error": "Could not save settings: {error}",
        "fastq_input": "Input FASTQ",
        "no_files": "No files selected",
        "choose_files": "Choose files…",
        "choose_folder": "Choose folder…",
        "choose": "Browse…",
        "metadata_file": "Metadata CSV",
        "reference_db_path": "Reference database",
        "primer_fwd": "Forward primer",
        "primer_rev": "Reverse primer",
        "trim_length_fwd": "Forward Trim Length (bp)",
        "trim_length_rev": "Reverse Trim Length (bp)",
        "rarefaction_label": "Rarefaction",
        "output_folder": "Output folder",
        "options": "Parameters — check primers and lengths for your experiment",
        "group_column": "Group column",
        "reference_group": "Reference group",
        "asv_min": "ASV min (bp)",
        "asv_max": "ASV max (bp)",
        "quality_threshold": "Quality threshold Q",
        "custom_depth": "Depth for custom mode",
        "run_ancombc": "Run ANCOM-BC (step 16)",
        "control_options": "Control sample handling",
        "metadata_hint": "Metadata: first column is SampleID matching FASTQ names; requires ≥2 groups, ≥2 samples/group.",
        "check": "Check environment / data",
        "analyze": "Run analysis",
        "stop": "Stop",
        "open_results": "Open results",
        "ready": "Ready. Check the environment before running.",
        "browse_output": "Select output folder",
        "browse_key": "Select {label}",
        "browse_fastq": "Select all R1 and R2 files (Ctrl/Shift)",
        "all_files": "All files",
        "dialog_error": "Could not open the file picker: {error}",
        "selected_files": "{count} FASTQ files selected",
        "metadata_title": "Metadata",
        "invalid_file": "Select a valid file for {key}.",
        "select_output": "Select an output folder.",
        "invalid_group": "The group column must exist and use a name such as Group/Treatment (letters, numbers, underscore).",
        "missing_metadata": "FASTQ samples missing from metadata: {samples}",
        "invalid_group_value": "Group is empty, NA, or has surrounding whitespace for sample {sample}.",
        "minimum_groups": "Diversity analysis requires at least 2 groups with at least 2 samples per group.",
        "invalid_reference": "The reference group must be present among the selected FASTQ samples.",
        "invalid_decontam_controls": "decontam requires at least one selected Negative Control and one non-negative sample.",
        "integer": "{key} must be an integer.",
        "positive": "{key} must be greater than zero.",
        "asv_order": "ASV min must not exceed ASV max.",
        "rarefaction": "Choose rarefaction auto_rarefy, min, or custom.",
        "primer": "{key} must be a DNA/IUPAC sequence with no spaces.",
        "invalid_primer_method": "Choose a valid primer removal method.",
        "nonnegative": "{key} must be zero or greater.",
        "scripts": "The pipeline folder beside this Python file must contain all 17 scripts from 00 to 16.",
        "cannot_run": "Cannot run",
        "checking": "Checking…",
        "analyzing": "Analyzing…",
        "run_folder": "Run folder: {path}",
        "launch_error": "Launch error: {error}",
        "stopped_status": "Stopped. Current results are incomplete.",
        "check_success": "Check completed successfully.",
        "analysis_success": "Analysis completed.",
        "complete": "Complete",
        "failed": "Failed (code {code}). See the log below.",
        "pipeline_stopped": "Pipeline stopped",
        "log_at": "See logs/pipeline.log in:\n{path}",
        "stopping": "Stopping…",
        "results": "Results",
        "running": "Running",
        "stop_before_close": "Stop the pipeline before closing the window.",
}

try:
    MARKER_PROFILES = json.loads(MARKER_PROFILES_PATH.read_text(encoding="utf-8"))
except (OSError, ValueError) as exc:
    raise RuntimeError(f"Could not load marker profiles from {MARKER_PROFILES_PATH}: {exc}") from exc


def get_marker_profile(marker):
    marker = str(marker).upper()
    try:
        profile = MARKER_PROFILES[marker]
    except KeyError:
        choices = ", ".join(MARKER_PROFILES)
        raise ValueError(
            f"Unsupported marker {marker!r}. Choose a configured marker: {choices}."
        ) from None
    if not isinstance(profile, dict):
        raise ValueError(f"Marker profile {marker!r} must be an object.")
    return profile


def saved_asv_length_filter_enabled(settings):
    if "asv_length_filter_enabled" in settings:
        return _parse_bool(settings["asv_length_filter_enabled"])
    return bool(settings) and (
        "asv_min_len" in settings or "asv_max_len" in settings
    )


def message(message_id, **values):
    return MESSAGES[message_id].format(**values)


def save_settings(settings):
    SETTINGS_DIR.mkdir(parents=True, exist_ok=True)
    SETTINGS.write_text(json.dumps(settings, ensure_ascii=False, indent=2), encoding="utf-8")


def load_settings():
    for path in (SETTINGS, LEGACY_APP_SETTINGS, LEGACY_SETTINGS):
        try:
            saved = json.loads(path.read_text(encoding="utf-8"))
            if isinstance(saved, dict):
                return saved
        except (OSError, ValueError, TypeError):
            continue
    return {}


def write_configuration_report(config, path):
    """Write an ordinary UTF-8 text report, separate from R's runtime JSON."""
    sections = {
        "Project": ("marker", "fastq_dir", "metadata_file", "reference_db_path",
                    "output_dir", "group_column", "reference_group", "steps"),
        "Read processing": ("primer_method", "primer_fwd", "primer_rev", "trim_left_fwd",
                            "trim_left_rev", "automatic_truncation", "trunc_qmin", "manual_read_length", "trunc_len",
                            "max_ee", "max_n", "trunc_q", "min_len", "remove_phix", "pooling",
                            "chimera_method", "mergepairs_strategy", "mergepairs_min_overlap",
                            "mergepairs_max_mismatch", "asv_length_filter_enabled",
                            "asv_min_len", "asv_max_len"),
        "Experimental controls": ("mock_control_column", "mock_control_level",
                                  "negative_control_column", "negative_control_level",
                                  "run_decontam", "exclude_mock_controls", "exclude_controls"),
        "Downstream analysis": ("rarefaction_mode", "custom_rarefaction_depth", "max_sample_loss",
                                "plateau_gain", "plateau_samples", "read_step", "beta_distances",
                                "alpha_measures", "run_alpha_statistics",
                                "ordination", "run_permanova", "run_betadisper",
                                "permanova_permutations", "differential_method"),
        "Plot output": ("figure_format", "abundance_plot_type", "heatmap_rank", "heatmap_top_n",
                        "min_rank1_abundance_percent", "min_rank2_abundance_percent"),
        "Execution": ("threads", "run_dir", "rscript", "cutadapt"),
    }
    labels = {
        "fastq_dir": "FASTQ folder", "metadata_file": "Metadata CSV",
        "reference_db_path": "Reference database", "group_column": "Comparison group",
        "primer_fwd": "Forward primer", "primer_rev": "Reverse primer",
        "trim_left_fwd": "Forward leading bases removed", "trim_left_rev": "Reverse leading bases removed",
        "trunc_qmin": "Automatic quality threshold",
        "manual_read_length": "Requested read lengths after primer removal (forward, reverse)",
        "trunc_len": "Internal DADA2 truncation positions (forward, reverse)",
        "max_ee": "Maximum expected errors (forward, reverse)", "max_n": "Maximum ambiguous bases",
        "trunc_q": "Per-read quality cutoff", "min_len": "Minimum read length (bp)",
        "step_1": "Read processing", "step_2": "Taxonomy", "step_3": "Downstream analysis",
        "min_rank1_abundance_percent": "First rank collapse threshold",
        "min_rank2_abundance_percent": "Second rank collapse threshold",
    }
    lines = ["BIOME - Run configuration", "=" * 40,
             f"Report created: {datetime.now().astimezone().isoformat(timespec='seconds')}", ""]

    def append_value(key, value, indent=""):
        label = labels.get(key, key.replace("_", " ").capitalize())
        if isinstance(value, dict):
            lines.append(f"{indent}{label}:")
            for child, item in value.items():
                append_value(child, item, indent + "  ")
            return
        if isinstance(value, bool):
            value = "Yes" if value else "No"
        elif value is None or value == "":
            value = "Not set"
        elif key.endswith("_percent"):
            value = f"{float(value):g}%"
        elif key in ("max_sample_loss", "plateau_gain", "plateau_samples"):
            value = f"{float(value) * 100:g}%"
        elif isinstance(value, (list, tuple)):
            value = ", ".join(str(item) for item in value) or "None"
        lines.append(f"{indent}{label}: {value}")

    reported = {"config_path", "marker_profile", "fastq_files", "sample_id_map"}
    for heading, keys in sections.items():
        lines.extend([heading, "-" * len(heading)])
        for key in keys:
            if key in config:
                append_value(key, config[key])
                reported.add(key)
        lines.append("")
    remaining = sorted(set(config) - reported)
    if remaining:
        lines.extend(["Additional parameters", "---------------------"])
        for key in remaining:
            append_value(key, config[key])
        lines.append("")
    lines.extend(["Selected input files", "--------------------",
                  "Sample-to-file mapping: tables/input_path.csv"])
    lines.extend(str(item) for item in config.get("fastq_files", []))
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def pair_files(files):
    """Reject missing mates, duplicate SampleIDs, and unsupported names."""
    pairs = {}
    for filename in sorted(set(files)):
        path = Path(filename).resolve()
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(message("file_missing", path=path))
        match = PAIR_PATTERN.fullmatch(path.name)
        if not match:
            raise ValueError(message("fastq_name", name=path.name))
        sample, direction = match.group(1), match.group(2) or match.group(3)
        if sample in pairs and direction in pairs[sample]:
            raise ValueError(message("duplicate_sample", sample=sample, direction=direction))
        pairs.setdefault(sample, {})[direction] = path.as_posix()
    if not pairs:
        raise ValueError(message("select_fastq"))
    incomplete = [sample for sample, mates in pairs.items() if set(mates) != {"1", "2"}]
    if incomplete:
        raise ValueError(message("missing_mates", samples=", ".join(incomplete[:15])))
    # Step 02 sorts full paths before comparing names. Keep that contract.
    forward = sorted(mates["1"] for mates in pairs.values())
    reverse = sorted(mates["2"] for mates in pairs.values())
    forward_ids = [PAIR_PATTERN.fullmatch(Path(p).name).group(1) for p in forward]
    reverse_ids = [PAIR_PATTERN.fullmatch(Path(p).name).group(1) for p in reverse]
    if forward_ids != reverse_ids:
        raise ValueError(message("pair_order"))
    return pairs


def is_fastq_path(path):
    return Path(path).name.lower().endswith(FASTQ_SUFFIXES)


def _raise_walk_error(error):
    raise error


def count_fastq_input_files(input_path):
    """Count FASTQ files in a folder or ZIP without blocking the caller's UI."""
    input_path = Path(input_path)
    if input_path.is_file() and input_path.suffix.lower() == ".zip":
        archives = [input_path]
        count = 0
    elif input_path.is_dir():
        archives = []
        count = 0
        for directory, _, names in os.walk(input_path, onerror=_raise_walk_error):
            for name in names:
                path = Path(directory) / name
                if is_fastq_path(path):
                    count += 1
                elif path.suffix.lower() == ".zip":
                    archives.append(path)
    else:
        archives = []
        count = int(input_path.is_file() and is_fastq_path(input_path))
    return count + sum(len(_zip_fastq_members(archive)) for archive in archives)


def _windows_fastq_input_dialog(parent, title, initial_dir):
    from ctypes import wintypes

    class BrowseInfo(ctypes.Structure):
        _fields_ = [
            ("hwndOwner", wintypes.HWND),
            ("pidlRoot", ctypes.c_void_p),
            ("pszDisplayName", wintypes.LPWSTR),
            ("lpszTitle", wintypes.LPCWSTR),
            ("ulFlags", wintypes.UINT),
            ("lpfn", ctypes.c_void_p),
            ("lParam", ctypes.c_ssize_t),
            ("iImage", ctypes.c_int),
        ]

    bffm_initialized = 1
    bffm_set_selection_w = 0x400 + 103
    initial_path = ctypes.create_unicode_buffer(initial_dir)
    user32 = ctypes.windll.user32
    user32.SendMessageW.argtypes = [
        wintypes.HWND, wintypes.UINT, wintypes.WPARAM, wintypes.LPARAM,
    ]
    user32.SendMessageW.restype = wintypes.LPARAM

    callback_type = ctypes.WINFUNCTYPE(
        wintypes.INT, wintypes.HWND, wintypes.UINT, wintypes.LPARAM, wintypes.LPARAM
    )

    def browse_callback(hwnd, message_id, _wparam, _lparam):
        if message_id == bffm_initialized:
            user32.SendMessageW(
                hwnd,
                bffm_set_selection_w,
                1,
                ctypes.cast(initial_path, ctypes.c_void_p).value,
            )
        return 0

    callback = callback_type(browse_callback)
    display_name = ctypes.create_unicode_buffer(32768)
    info = BrowseInfo(
        hwndOwner=parent.winfo_id(),
        pidlRoot=None,
        pszDisplayName=ctypes.cast(display_name, wintypes.LPWSTR),
        lpszTitle=title,
        ulFlags=0x10 | 0x40 | 0x4000,
        lpfn=ctypes.cast(callback, ctypes.c_void_p),
        lParam=0,
        iImage=0,
    )
    shell32 = ctypes.windll.shell32
    shell32.SHBrowseForFolderW.restype = ctypes.c_void_p
    shell32.SHBrowseForFolderW.argtypes = [ctypes.POINTER(BrowseInfo)]
    ole32 = ctypes.windll.ole32
    ole32.CoInitializeEx.argtypes = [ctypes.c_void_p, wintypes.DWORD]
    ole32.CoInitializeEx.restype = ctypes.c_long
    ole32.CoUninitialize.argtypes = []
    ole32.CoUninitialize.restype = None
    result = ole32.CoInitializeEx(None, 0x2)
    if result < 0:
        raise OSError(result, "Could not initialize the Windows file picker.")
    pidl = None
    try:
        pidl = shell32.SHBrowseForFolderW(ctypes.byref(info))
        if not pidl:
            return ""
        selected_path = ctypes.create_unicode_buffer(32768)
        shell32.SHGetPathFromIDListW.argtypes = [ctypes.c_void_p, wintypes.LPWSTR]
        shell32.SHGetPathFromIDListW.restype = wintypes.BOOL
        if not shell32.SHGetPathFromIDListW(pidl, selected_path):
            raise OSError("Could not retrieve the selected FASTQ folder or ZIP archive.")
        return selected_path.value
    finally:
        if pidl:
            ole32.CoTaskMemFree.argtypes = [ctypes.c_void_p]
            ole32.CoTaskMemFree.restype = None
            ole32.CoTaskMemFree(pidl)
        ole32.CoUninitialize()


def _path_is_within(path, parent):
    try:
        Path(path).resolve().relative_to(Path(parent).resolve())
        return True
    except ValueError:
        return False


def _zip_fastq_members(archive):
    try:
        with zipfile.ZipFile(archive) as zipped:
            return [
                item.filename for item in zipped.infolist()
                if not item.is_dir() and is_fastq_path(item.filename)
            ]
    except (OSError, zipfile.BadZipFile) as exc:
        raise ValueError(f"Could not read FASTQ ZIP archive {archive}: {exc}") from exc


def _extract_zip_safely(archive, destination):
    destination = Path(destination).resolve()
    with zipfile.ZipFile(archive) as zipped:
        for item in zipped.infolist():
            member = Path(item.filename.replace("\\", "/"))
            if member.is_absolute() or ".." in member.parts:
                raise ValueError(f"Unsafe path in ZIP archive {archive}: {item.filename}")
            target = (destination / member).resolve()
            if not _path_is_within(target, destination):
                raise ValueError(f"Unsafe path in ZIP archive {archive}: {item.filename}")
        zipped.extractall(destination)


def resolve_fastq_input(input_path, output_dir):
    """Return the read directory and files, extracting FASTQ ZIPs without expanding gzip files."""
    input_path = Path(input_path).expanduser()
    output_dir = Path(output_dir).expanduser().resolve()
    if input_path.is_file():
        if input_path.suffix.lower() == ".zip":
            archives = [input_path.resolve()]
            fastq_paths = []
        elif is_fastq_path(input_path):
            archives = []
            fastq_paths = [input_path.resolve()]
        else:
            raise ValueError("Select an existing FASTQ folder or ZIP archive.")
    elif input_path.is_dir():
        input_path = input_path.resolve()
        archives = []
        fastq_paths = []
        if output_dir == input_path:
            excluded = [output_dir / "extracted_fastq"]
        elif _path_is_within(output_dir, input_path):
            excluded = [output_dir]
        else:
            excluded = []
        for directory, dirs, names in os.walk(input_path):
            current = Path(directory).resolve()
            dirs[:] = [
                name for name in dirs
                if not any(
                    _path_is_within(current / name, excluded_path)
                    for excluded_path in excluded
                )
            ]
            for name in names:
                path = current / name
                if is_fastq_path(path):
                    fastq_paths.append(path)
                elif path.suffix.lower() == ".zip":
                    archives.append(path)
    else:
        raise ValueError("Select an existing FASTQ folder or ZIP archive.")

    if not fastq_paths:
        readable_archives = [
            (archive, _zip_fastq_members(archive))
            for archive in sorted(archives)
        ]
        readable_archives = [
            (archive, members) for archive, members in readable_archives if members
        ]
    else:
        readable_archives = []
    if not fastq_paths and not readable_archives:
        raise ValueError("The FASTQ folder or ZIP archive contains no FASTQ files.")

    extracted_paths = []
    if readable_archives:
        extraction_root = output_dir / "extracted_fastq"
        extraction_root.mkdir(parents=True, exist_ok=True)
        run_extract_dir = Path(tempfile.mkdtemp(prefix="run_", dir=extraction_root))
        for index, (archive, _members) in enumerate(readable_archives, start=1):
            destination = run_extract_dir / f"{index:03d}_{archive.stem}"
            destination.mkdir()
            _extract_zip_safely(archive, destination)
            extracted_paths.extend(
                path.resolve() for path in destination.rglob("*")
                if path.is_file() and is_fastq_path(path)
            )
        read_paths = fastq_paths + extracted_paths
        read_directory = Path(os.path.commonpath([str(path.parent) for path in extracted_paths]))
    else:
        read_paths = fastq_paths
        read_directory = Path(os.path.commonpath([str(path.parent) for path in read_paths]))

    return read_directory.resolve(), sorted(str(path.resolve()) for path in read_paths)


def match_metadata_sample_ids(fastq_sample_ids, metadata_sample_ids):
    """Match Illumina lane-suffixed names to metadata IDs when needed."""
    metadata_ids = set(metadata_sample_ids)
    sample_id_map = {}
    for sample_id in fastq_sample_ids:
        metadata_id = sample_id
        if metadata_id not in metadata_ids:
            candidate = ILLUMINA_LANE_SUFFIX.sub("", sample_id)
            if candidate != sample_id and candidate in metadata_ids:
                metadata_id = candidate
        if metadata_id in sample_id_map.values():
            raise ValueError(
                f"Multiple FASTQ sample names match metadata SampleID {metadata_id!r}."
            )
        sample_id_map[sample_id] = metadata_id
    return sample_id_map


def parse_control_prefixes(value):
    """Parse comma-separated sample-ID prefixes, ignoring empty entries."""
    return list(dict.fromkeys(prefix.strip() for prefix in value.split(",") if prefix.strip()))


def is_control_sample(sample_id, prefixes):
    sample_id = sample_id.casefold()
    return any(sample_id.startswith(prefix.casefold()) for prefix in prefixes)


def read_metadata(path):
    path = Path(path)
    if path.suffix.casefold() != ".csv":
        raise ValueError("Select a metadata file with a .csv extension.")
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, strict=True)
        columns = reader.fieldnames or []
        try:
            rows = list(reader)
        except csv.Error as exc:
            raise ValueError(f"Could not read metadata CSV: {exc}") from exc
    if len(columns) < 2 or "SampleID" not in columns or not rows:
        raise ValueError(message("metadata_csv"))
    if len(set(columns)) != len(columns):
        raise ValueError(message("duplicate_columns"))
    if any(None in row or any(value is None for value in row.values()) for row in rows):
        raise ValueError(message("metadata_row"))
    sample_ids = [row["SampleID"] for row in rows]
    if any(not value.strip() for value in sample_ids) or len(set(sample_ids)) != len(sample_ids):
        raise ValueError(message("sample_ids"))
    return columns, rows


def _parse_bool(value):
    if isinstance(value, bool):
        return value
    if value is None:
        return False
    return str(value).strip().lower() in {"1", "true", "yes", "y", "on"}


def create_app_temp_dir(root=APP_TEMP_ROOT):
    root = Path(root)
    if root.is_symlink():
        raise ValueError("BIOME temporary folder must not be a symbolic link.")
    root.mkdir(parents=True, exist_ok=True)
    directory = Path(tempfile.mkdtemp(prefix="session_", dir=root))
    (directory / ".owner").write_text(str(os.getpid()), encoding="ascii")
    return directory


def _process_is_running(pid):
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    except OSError as exc:
        return exc.errno != errno.ESRCH
    return True


def _remove_owned_temp_dir(directory, root=APP_TEMP_ROOT, expected_pid=None):
    root = Path(root)
    if root.is_symlink():
        raise ValueError("BIOME temporary folder must not be a symbolic link.")
    root = root.resolve()
    directory = Path(directory)
    if directory.is_symlink() or directory.resolve().parent != root:
        raise ValueError("Temporary directory is outside the BIOME temp folder.")
    if not directory.name.startswith("session_"):
        raise ValueError("Refusing to remove an unrecognized temporary directory.")
    owner_file = directory / ".owner"
    try:
        owner_pid = int(owner_file.read_text(encoding="ascii").strip())
    except (OSError, ValueError) as exc:
        raise ValueError(f"Temporary directory has no valid ownership marker: {directory}") from exc
    if owner_pid != (os.getpid() if expected_pid is None else expected_pid):
        raise ValueError("Temporary directory is not owned by this process.")
    size = 0

    def raise_walk_error(error):
        raise error

    for current, dirs, files in os.walk(
        directory, topdown=False, followlinks=False, onerror=raise_walk_error
    ):
        for name in files:
            path = Path(current) / name
            if name == ".owner" or path.suffix.casefold() == ".rds":
                continue
            if not path.is_symlink():
                size += path.stat().st_size
            path.unlink()
        for name in dirs:
            path = Path(current) / name
            if path.suffix.casefold() == ".rds":
                continue
            if path.is_symlink():
                path.unlink()
                continue
            try:
                path.rmdir()
            except OSError as exc:
                if exc.errno not in {errno.ENOTEMPTY, errno.EEXIST}:
                    raise

    has_rds_files = any(
        path.suffix.casefold() == ".rds"
        for current, _dirs, files in os.walk(
            directory, followlinks=False, onerror=raise_walk_error
        )
        for path in (Path(current) / name for name in files)
    )
    if not has_rds_files:
        owner_file = directory / ".owner"
        size += owner_file.stat().st_size
        owner_file.unlink()
        directory.rmdir()
    return size


def clear_orphaned_temp_dirs(root=APP_TEMP_ROOT):
    root = Path(root)
    if root.is_symlink():
        raise ValueError("BIOME temporary folder must not be a symbolic link.")
    if not root.exists():
        return 0
    freed = 0
    for directory in root.iterdir():
        if (not directory.name.startswith("session_") or not directory.is_dir()
                or directory.is_symlink()):
            continue
        owner_file = directory / ".owner"
        try:
            owner_pid = int(owner_file.read_text(encoding="ascii").strip())
        except (OSError, ValueError):
            continue
        if _process_is_running(owner_pid):
            continue
        freed += _remove_owned_temp_dir(directory, root=root, expected_pid=owner_pid)
    return freed


def format_size(size):
    if size >= 1024 ** 3:
        return f"{size / (1024 ** 3):.2f} GB"
    return f"{size / (1024 ** 2):.2f} MB"


def report_startup_failure(exc_type, exc_value, exc_traceback, parent=None):
    details = "".join(traceback.format_exception(exc_type, exc_value, exc_traceback))
    log_path, log_error = _write_startup_log(details)
    if log_error is not None:
        details += f"\nCould not write startup log: {log_error}\n"

    message = (
        "BIOME encountered an unexpected error and could not start.\n\n"
        f"{exc_type.__name__}: {exc_value}"
    )
    if log_path is not None:
        message += f"\n\nDetails were saved to:\n{log_path}"
    else:
        message += f"\n\n{details}"

    try:
        messagebox.showerror("BIOME could not start", message, parent=parent)
    except tk.TclError:
        _show_native_startup_error(message)


class PipelineApp:
    def __init__(self, root):
        self.root = root
        if os.name == "nt":
            self.root.iconbitmap(str(BASE / "assets" / "biome.ico"))
        # Fit the usable desktop area, including smaller laptop displays.
        available_width = self.root.winfo_screenwidth()
        available_height = self.root.winfo_screenheight()
        if sys.platform == "win32":
            from ctypes import wintypes
            work_area = wintypes.RECT()
            if ctypes.windll.user32.SystemParametersInfoW(0x0030, 0, ctypes.byref(work_area), 0):
                available_width = work_area.right - work_area.left
                available_height = work_area.bottom - work_area.top
        width = max(1, min(1080, available_width - 40))
        height = max(1, min(760, available_height - 80))
        self.root.geometry(f"{width}x{height}")
        self.root.minsize(min(800, width), min(480, height))
        self.events = queue.Queue()
        self.process = None
        self.busy = False
        self.cancelled = threading.Event()
        self.temp_dir = None
        self.fastq_files = []
        self.last_run = None
        self.vars = {}
        self.controls = []
        self.status_key = "ready"
        self.status_args = {}
        self.status = tk.StringVar()
        self.saved = load_settings()
        # Restore runtime setup only; each analysis starts with fresh inputs.
        self.saved = {key: self.saved[key] for key in
                      ("rscript", "cutadapt", "setup_completed") if key in self.saved}
        marker = str(self.saved.get("marker", "16S")).upper()
        marker_profile = get_marker_profile(marker)
        legacy_asv_filter = (
            "asv_length_filter_enabled" not in self.saved
            and saved_asv_length_filter_enabled(self.saved)
        )
        saved_reference = self.saved.get(
            "reference_db_path", self.saved.get("silva_ref_path", "")
        )
        saved_rank1_threshold = self.saved.get(
            "min_rank1_abundance_percent",
            self.saved.get("min_phylum_abundance_percent", "1.0"),
        )
        saved_rank2_threshold = self.saved.get(
            "min_rank2_abundance_percent",
            self.saved.get("min_genus_abundance_percent", "2.0"),
        )
        defaults = {
            "rscript": self.saved.get("rscript") or find_rscript(),
            "cutadapt": self.saved.get("cutadapt") or find_cutadapt(BASE),
            "fastq_dir": "", "metadata_file": "", "reference_db_path": saved_reference,
            "marker": marker,
            "output_dir": "", "group_column": "",
            "reference_group": "",
            "primer_fwd": "",
            "primer_rev": "",
            "primer_mode": marker_profile["default_primer_mode"],
            "trim_left_fwd": "",
            "trim_left_rev": "",
            "truncation_mode": "automatic", "trunc_qmin": "30",
            "trunc_len_fwd": "", "trunc_len_rev": "",
            "max_ee_fwd": "2", "max_ee_rev": "2", "pooling": "false",
            "max_n": "0", "trunc_q": "2",
            "min_len": str(marker_profile["default_min_read_length"]),
            "asv_min_len": "", "asv_max_len": "",
            "chimera_method": "consensus", "mergepairs_strategy": "merge",
            "mergepairs_min_overlap": str(marker_profile["default_merge_overlap"]),
            "mergepairs_max_mismatch": "0",
            "rarefaction_mode": "auto_rarefy", "custom_rarefaction_depth": "",
            "max_sample_loss_percent": "10", "plateau_gain_percent": "2",
            "plateau_samples_percent": "90", "read_step": "500",
            "permutations": "999", "ordination": "PCoA", "differential_method": "ANCOM-BC",
            "min_rank1_abundance_percent": str(saved_rank1_threshold),
            "min_rank2_abundance_percent": str(saved_rank2_threshold),
            "abundance_plot_type": "barplot", "heatmap_rank": "Genus", "heatmap_top_n": "30",
            "figure_format": "PNG",
            "mock_control_mode": "None", "mock_control_column": "", "mock_control_level": "",
            "negative_control_mode": "None", "negative_control_column": "",
            "negative_control_level": "",
        }
        defaults.update({k: str(v) for k, v in self.saved.items() if k in defaults})
        defaults["marker"] = marker
        self.legacy_asv_filter = legacy_asv_filter
        self.asv_length_filter_enabled = tk.BooleanVar(
            value=saved_asv_length_filter_enabled(self.saved)
        )
        old_primer_method = str(self.saved.get("primer_method", ""))
        if "primer_mode" not in self.saved and old_primer_method:
            defaults["primer_mode"] = (
                "length" if "trim_left" in old_primer_method.casefold()
                else "sequences"
            )
        if "truncation_mode" not in self.saved:
            if "disable_truncation" in self.saved:
                defaults["truncation_mode"] = (
                    "automatic" if _parse_bool(self.saved["disable_truncation"]) else "manual"
                )
        if self.saved.get("mergepairs_strategy") == "concatenate":
            defaults["mergepairs_strategy"] = "concatenate"
        if defaults["mergepairs_strategy"] == "merge":
            defaults["mergepairs_strategy"] = "Merge paired reads"
        elif defaults["mergepairs_strategy"] == "concatenate":
            defaults["mergepairs_strategy"] = "Concatenate reads"
        legacy_rarefaction = {
            "auto_rarefy": "Automatic",
            "min": "Minimum depth",
            "custom": "Custom depth",
        }
        defaults["rarefaction_mode"] = legacy_rarefaction.get(
            defaults["rarefaction_mode"], defaults["rarefaction_mode"]
        )
        if self.saved.get("mock_control_column"):
            defaults["mock_control_mode"] = "Metadata column"
        if self.saved.get("negative_control_column"):
            defaults["negative_control_mode"] = "Metadata column"
        chimera_methods = {"consensus": "Consensus", "pooled": "Pooled", "none": "None"}
        defaults["chimera_method"] = chimera_methods.get(
            defaults["chimera_method"].casefold(), defaults["chimera_method"]
        )
        self.fastq_files = self.saved.get("fastq_files", [])
        if not defaults["fastq_dir"] and self.fastq_files:
            defaults["fastq_dir"] = str(Path(self.fastq_files[0]).parent)
        for key in ("fastq_dir", "metadata_file", "reference_db_path", "output_dir", "rscript"):
            if defaults.get(key):
                defaults[key] = str(Path(defaults[key]))
        self.vars = {k: tk.StringVar(value=v) for k, v in defaults.items()}
        self._active_marker = marker
        self.step_vars = {
            f"step_{number}": tk.BooleanVar(value=_parse_bool(self.saved.get(f"step_{number}", True)))
            for number in (1, 2, 3)
        }
        self.run_betadisper = tk.BooleanVar(value=_parse_bool(self.saved.get("run_betadisper", True)))
        self.run_permanova = tk.BooleanVar(value=_parse_bool(self.saved.get("run_permanova", True)))
        self.run_alpha_statistics = tk.BooleanVar(
            value=_parse_bool(self.saved.get("run_alpha_statistics", True))
        )
        self.alpha_vars = {
            metric: tk.BooleanVar(value=metric in self.saved.get(
                "alpha_measures", ["Observed", "Chao1", "Shannon", "Simpson"]
            ))
            for metric in ("Observed", "Chao1", "Shannon", "Simpson", "Pielou")
        }
        self.remove_phix = tk.BooleanVar(value=_parse_bool(self.saved.get("remove_phix", True)))
        self.beta_vars = {
            "bray": tk.BooleanVar(value=_parse_bool(self.saved.get("beta_bray", True))),
            "wunifrac": tk.BooleanVar(value=_parse_bool(self.saved.get("beta_wunifrac", False))),
            "unifrac": tk.BooleanVar(value=_parse_bool(self.saved.get("beta_unifrac", False))),
            "jaccard": tk.BooleanVar(value=_parse_bool(self.saved.get("beta_jaccard", False))),
        }
        self.metadata_columns = []
        self.advanced_open = False
        self.rarefaction_advanced_open = False
        self.build_ui()
        self.refresh_fastq_label()
        if self.vars["metadata_file"].get():
            self.refresh_metadata(show_errors=False)
        self.root.protocol("WM_DELETE_WINDOW", self.close)
        self.root.after(100, self.poll)

    def text(self, message_id, **values):
        return message(message_id, **values)

    def set_status(self, key, **values):
        self.status_key = key
        self.status_args = values
        self.status.set(self.text(key, **values))

    def update_marker_profile(self):
        marker = self.vars["marker"].get().upper()
        profile = get_marker_profile(marker)
        if marker != self._active_marker:
            for key in ("primer_fwd", "primer_rev", "trim_left_fwd", "trim_left_rev",
                        "trunc_len_fwd", "trunc_len_rev", "asv_min_len", "asv_max_len"):
                self.vars[key].set("")
            self.vars["reference_db_path"].set("")
            if profile["default_primer_mode"] == "length":
                self.vars["primer_mode"].set("length")
                if hasattr(self, "sequence_fields"):
                    self.update_primer_fields()
            self._active_marker = marker
        if hasattr(self, "reference_suggestion_label"):
            self.reference_suggestion_label.configure(
                text=f"Suggested reference: {profile['reference_suggestion']}"
            )
        if hasattr(self, "abundance_rank_labels"):
            for label, rank in zip(
                self.abundance_rank_labels,
                profile["relative_abundance_ranks"],
            ):
                label.configure(text=f"Collapse {rank} taxa below")
        if hasattr(self, "heatmap_rank_combo"):
            rank1 = profile["relative_abundance_ranks"][0]
            self.heatmap_rank_combo.configure(values=("Genus", rank1))
            if self.vars["heatmap_rank"].get() in {"Phylum", "Division"}:
                self.vars["heatmap_rank"].set(rank1)

    def update_asv_length_fields(self):
        state = "normal" if self.asv_length_filter_enabled.get() else "disabled"
        self.asv_min_len_entry.configure(state=state)
        self.asv_max_len_entry.configure(state=state)

    def build_ui(self):
        self.controls = []
        self.root.title(self.text("window_title"))
        self.main_frame = ttk.Frame(self.root)
        self.main_frame.pack(fill="both", expand=True)
        self.main_frame.columnconfigure(0, weight=1)
        self.main_frame.rowconfigure(0, weight=1)

        self.notebook = ttk.Notebook(self.main_frame)
        self.notebook.grid(row=0, column=0, sticky="nsew")
        project = self.scrollable_tab("Project")
        processing = self.scrollable_tab("Read processing")
        downstream = self.scrollable_tab("Downstream analysis")

        marker_frame = ttk.LabelFrame(project, text="Marker", padding=9)
        marker_frame.grid(row=0, column=0, sticky="ew", pady=(0, 8))
        for column, marker in enumerate(MARKER_PROFILES):
            radio = ttk.Radiobutton(
                marker_frame,
                text=marker,
                variable=self.vars["marker"],
                value=marker,
                command=self.update_marker_profile,
            )
            radio.grid(row=0, column=column, sticky="w", padx=(0, 18))
            self.controls.append(radio)
        self.reference_suggestion_label = ttk.Label(
            marker_frame, foreground="#536273", wraplength=800, justify="left"
        )
        self.reference_suggestion_label.grid(
            row=1, column=0, columnspan=2, sticky="w", pady=(6, 0)
        )
        self.update_marker_profile()

        inputs = ttk.LabelFrame(project, text="Project files", padding=9)
        inputs.grid(row=1, column=0, sticky="ew", pady=(0, 8))
        inputs.columnconfigure(1, weight=1)
        file_rows = (
            ("FASTQ folder:", "fastq_dir", None),
            ("Metadata file", "metadata_file", [("CSV", "*.csv"), ("All files", "*")]),
            ("Output folder", "output_dir", None),
            ("Reference database", "reference_db_path", [("FASTA", "*.fa *.fasta *.fa.gz *.fasta.gz"), ("All files", "*")]),
        )
        for row, (label, key, types) in enumerate(file_rows):
            self.path_row(inputs, row, label, key, types)
        self.fastq_label = ttk.Label(inputs, text="", font=("Segoe UI", 9))
        self.fastq_label.grid(row=len(file_rows), column=0, columnspan=3, sticky="e")

        project_settings = ttk.LabelFrame(project, text="Analysis settings", padding=9)
        project_settings.grid(row=2, column=0, sticky="ew", pady=(0, 8))
        for col in (1, 3):
            project_settings.columnconfigure(col, weight=1)
        ttk.Label(project_settings, text="Comparison group").grid(row=0, column=0, sticky="w", padx=(0, 8), pady=5)
        self.group_combo = ttk.Combobox(
            project_settings, textvariable=self.vars["group_column"], state="readonly"
        )
        self.group_combo.grid(row=0, column=1, sticky="ew", padx=(0, 18), pady=5)
        self.group_combo.bind("<<ComboboxSelected>>", lambda _event: self.refresh_groups())
        self.controls.append(self.group_combo)
        ttk.Label(project_settings, text="Reference group").grid(row=1, column=0, sticky="w", padx=(0, 8), pady=5)
        self.ref_combo = ttk.Combobox(
            project_settings, textvariable=self.vars["reference_group"], state="readonly"
        )
        self.ref_combo.grid(row=1, column=1, sticky="ew", padx=(0, 18), pady=5)
        self.controls.append(self.ref_combo)

        steps = ttk.LabelFrame(project, text="Pipeline steps", padding=9)
        steps.grid(row=3, column=0, sticky="ew", pady=(0, 8))
        self.step_checkboxes = {}
        for index, (key, label, detail) in enumerate((
            ("step_1", "Read processing", "Prepare high-quality sequence reads."),
            ("step_2", "Taxonomy", "Classify reads using the selected reference database."),
            ("step_3", "Downstream analysis", "Calculate diversity and compare groups."),
        )):
            row = index * 3
            checkbox = ttk.Checkbutton(
                steps, text=label, variable=self.step_vars[key],
                command=self.update_pipeline_steps,
            )
            checkbox.grid(row=row, column=0, sticky="w", pady=(2, 0))
            self.step_checkboxes[key] = checkbox
            ttk.Label(steps, text=detail).grid(
                row=row + 1, column=0, sticky="w", padx=(24, 0), pady=(0, 2)
            )
            if index < 2:
                ttk.Label(steps, text="↓", foreground="#536273").grid(
                    row=row + 2, column=0, sticky="w", padx=(12, 0), pady=(0, 2)
                )
            self.controls.append(checkbox)

        primer = ttk.LabelFrame(processing, text="Primer removal", padding=9)
        primer.grid(row=0, column=0, sticky="ew", pady=(0, 8))
        primer.columnconfigure(0, weight=1)
        ttk.Radiobutton(
            primer, text="I know the primer sequences", variable=self.vars["primer_mode"],
            value="sequences", command=self.update_primer_fields,
        ).grid(row=0, column=0, sticky="w", pady=(0, 7))
        self.sequence_fields = ttk.Frame(primer)
        self.sequence_fields.grid(row=1, column=0, sticky="ew", padx=(24, 0), pady=(0, 9))
        self.sequence_fields.columnconfigure(1, weight=1)
        self.sequence_fields.columnconfigure(3, weight=1)
        for col, label, key in (
            (0, "Forward primer", "primer_fwd"),
            (2, "Reverse primer", "primer_rev"),
        ):
            ttk.Label(self.sequence_fields, text=label).grid(row=0, column=col, sticky="w", padx=(0, 8))
            entry = ttk.Entry(self.sequence_fields, textvariable=self.vars[key])
            entry.grid(row=0, column=col + 1, sticky="ew", padx=(0, 16))
            self.controls.append(entry)
        ttk.Radiobutton(
            primer, text="I only know the primer length", variable=self.vars["primer_mode"],
            value="length", command=self.update_primer_fields,
        ).grid(row=2, column=0, sticky="w", pady=(0, 7))
        self.length_fields = ttk.Frame(primer)
        self.length_fields.grid(row=3, column=0, sticky="ew", padx=(24, 0))
        self.length_fields.columnconfigure(1, weight=1)
        self.length_fields.columnconfigure(3, weight=1)
        for col, label, key in (
            (0, "Forward primer length (bp)", "trim_left_fwd"),
            (2, "Reverse primer length (bp)", "trim_left_rev"),
        ):
            ttk.Label(self.length_fields, text=label).grid(row=0, column=col, sticky="w", padx=(0, 8))
            entry = ttk.Entry(self.length_fields, textvariable=self.vars[key], width=10)
            entry.grid(row=0, column=col + 1, sticky="w", padx=(0, 16))
            self.controls.append(entry)
        self.primer_suggestion_label = ttk.Label(
            primer, text="Enter primers or primer lengths matching your assay; these fields are required.",
            foreground="#536273"
        )
        self.primer_suggestion_label.grid(
            row=4, column=0, sticky="w", padx=(24, 0), pady=(8, 0)
        )

        quality = ttk.LabelFrame(processing, text="Quality truncation", padding=9)
        quality.grid(row=1, column=0, sticky="ew", pady=(0, 8))
        quality.columnconfigure(0, weight=1)
        ttk.Radiobutton(
            quality, text="Automatically determine truncation length", variable=self.vars["truncation_mode"],
            value="automatic", command=self.update_truncation_fields,
        ).grid(row=0, column=0, sticky="w", pady=(0, 6))
        automatic_fields = ttk.Frame(quality)
        automatic_fields.grid(row=1, column=0, sticky="w", padx=(24, 0), pady=(0, 9))
        ttk.Label(automatic_fields, text="Minimum quality score").pack(side="left", padx=(0, 10))
        self.trunc_qmin_entry = ttk.Entry(automatic_fields, textvariable=self.vars["trunc_qmin"], width=8)
        self.trunc_qmin_entry.pack(side="left")
        self.controls.append(self.trunc_qmin_entry)
        ttk.Label(
            automatic_fields,
            text="The software chooses read lengths from your quality profiles.",
        ).pack(side="left", padx=(12, 0))
        ttk.Radiobutton(
            quality, text="Read length after primer removal (bp)", variable=self.vars["truncation_mode"],
            value="manual", command=self.update_truncation_fields,
        ).grid(row=2, column=0, sticky="w", pady=(0, 6))
        self.manual_truncation_fields = ttk.Frame(quality)
        self.manual_truncation_fields.grid(row=3, column=0, sticky="ew", padx=(24, 0))
        self.manual_truncation_fields.columnconfigure(1, weight=1)
        self.manual_truncation_fields.columnconfigure(3, weight=1)
        ttk.Label(
            self.manual_truncation_fields,
            text="Excludes primer bases. BIOME converts these lengths automatically for either primer method.\n"
                 "Quality filtering may shorten or discard reads.",
            justify="left", wraplength=720,
        ).grid(row=1, column=0, columnspan=4, sticky="w", pady=(6, 0))
        for col, label, key in (
            (0, "Forward (bp)", "trunc_len_fwd"),
            (2, "Reverse (bp)", "trunc_len_rev"),
        ):
            ttk.Label(self.manual_truncation_fields, text=label).grid(row=0, column=col, sticky="w", padx=(0, 8))
            entry = ttk.Entry(self.manual_truncation_fields, textvariable=self.vars[key], width=10)
            entry.grid(row=0, column=col + 1, sticky="w", padx=(0, 16))
            self.controls.append(entry)

        self.advanced_toggle = ttk.Button(
            processing, text="Advanced filtering  +", command=self.toggle_advanced_filtering
        )
        self.advanced_toggle.grid(row=2, column=0, sticky="w", pady=(0, 8))
        self.controls.append(self.advanced_toggle)
        self.advanced_filtering = ttk.LabelFrame(processing, text="Advanced filtering", padding=9)
        self.advanced_filtering.grid(row=3, column=0, sticky="ew", pady=(0, 8))
        self.advanced_filtering.grid_remove()
        for col in (1, 3):
            self.advanced_filtering.columnconfigure(col, weight=1)
        ttk.Label(self.advanced_filtering, text="Chimera removal method").grid(row=0, column=0, sticky="w", padx=(0, 8), pady=4)
        self.chimera_combo = ttk.Combobox(
            self.advanced_filtering, textvariable=self.vars["chimera_method"],
            values=("Consensus", "Pooled", "None"), state="readonly", width=16,
        )
        self.chimera_combo.grid(row=0, column=1, sticky="ew", padx=(0, 18), pady=4)
        self.controls.append(self.chimera_combo)
        for col, label, key in (
            (0, "Maximum expected errors ? forward", "max_ee_fwd"),
            (2, "Maximum expected errors ? reverse", "max_ee_rev"),
            (0, "Maximum ambiguous bases", "max_n"),
            (2, "Minimum quality score for truncation", "trunc_q"),
        ):
            row = 1 if key.startswith("max_ee") else 2
            ttk.Label(self.advanced_filtering, text=label).grid(row=row, column=col, sticky="w", padx=(0, 8), pady=4)
            entry = ttk.Entry(self.advanced_filtering, textvariable=self.vars[key], width=10)
            entry.grid(row=row, column=col + 1, sticky="ew", padx=(0, 18), pady=4)
            self.controls.append(entry)
        phix = ttk.Checkbutton(
            self.advanced_filtering, text="Remove PhiX reads", variable=self.remove_phix
        )
        phix.grid(row=3, column=0, columnspan=2, sticky="w", pady=(8, 0))
        self.controls.append(phix)

        asv_length_filter = ttk.LabelFrame(
            processing, text="ASV length filter", padding=9
        )
        asv_length_filter.grid(row=4, column=0, sticky="ew", pady=(0, 8))
        ttk.Radiobutton(
            asv_length_filter,
            text="No ASV length filtering",
            variable=self.asv_length_filter_enabled,
            value=False,
            command=self.update_asv_length_fields,
        ).grid(row=0, column=0, columnspan=4, sticky="w", pady=(0, 5))
        ttk.Radiobutton(
            asv_length_filter,
            text="Filter ASVs by length",
            variable=self.asv_length_filter_enabled,
            value=True,
            command=self.update_asv_length_fields,
        ).grid(row=1, column=0, columnspan=4, sticky="w", pady=(0, 5))
        ttk.Label(asv_length_filter, text="Minimum ASV length (bp):").grid(
            row=2, column=0, sticky="w", padx=(24, 8), pady=4
        )
        self.asv_min_len_entry = ttk.Entry(
            asv_length_filter, textvariable=self.vars["asv_min_len"], width=12
        )
        self.asv_min_len_entry.grid(row=2, column=1, sticky="w", pady=4)
        ttk.Label(asv_length_filter, text="Maximum ASV length (bp):").grid(
            row=2, column=2, sticky="w", padx=(24, 8), pady=4
        )
        self.asv_max_len_entry = ttk.Entry(
            asv_length_filter, textvariable=self.vars["asv_max_len"], width=12
        )
        self.asv_max_len_entry.grid(row=2, column=3, sticky="w", pady=4)
        self.controls.extend((self.asv_min_len_entry, self.asv_max_len_entry))
        self.update_asv_length_fields()

        merging = ttk.LabelFrame(processing, text="Read merging", padding=9)
        merging.grid(row=5, column=0, sticky="ew", pady=(0, 8))
        merging.columnconfigure(1, weight=1)
        ttk.Label(merging, text="Merge strategy").grid(row=0, column=0, sticky="w", padx=(0, 10), pady=4)
        self.merge_combo = ttk.Combobox(
            merging, textvariable=self.vars["mergepairs_strategy"],
            values=("Merge paired reads", "Concatenate reads"), state="readonly", width=24,
        )
        self.merge_combo.grid(row=0, column=1, sticky="w", pady=4)
        self.merge_combo.bind("<<ComboboxSelected>>", lambda _event: self.update_merging_fields())
        self.controls.append(self.merge_combo)
        self.merge_fields = ttk.Frame(merging)
        self.merge_fields.grid(row=1, column=0, columnspan=2, sticky="ew", pady=(7, 0))
        self.merge_fields.columnconfigure(1, weight=1)
        self.merge_fields.columnconfigure(3, weight=1)
        for col, label, key in (
            (0, "Minimum overlap (bp)", "mergepairs_min_overlap"),
            (2, "Maximum mismatches", "mergepairs_max_mismatch"),
        ):
            ttk.Label(self.merge_fields, text=label).grid(row=0, column=col, sticky="w", padx=(0, 8), pady=4)
            entry = ttk.Entry(self.merge_fields, textvariable=self.vars[key], width=10)
            entry.grid(row=0, column=col + 1, sticky="w", padx=(0, 18), pady=4)
            self.controls.append(entry)

        controls = ttk.LabelFrame(processing, text="Experimental controls", padding=9)
        controls.grid(row=6, column=0, sticky="ew", pady=(0, 8))
        controls.columnconfigure(1, weight=1)
        controls.columnconfigure(3, weight=1)
        ttk.Label(controls, text="Mock community").grid(row=0, column=0, sticky="w", padx=(0, 8), pady=5)
        ttk.Label(controls, text="Negative controls").grid(row=1, column=0, sticky="w", padx=(0, 8), pady=5)
        self.control_widgets = {}
        for row, key in ((0, "mock_control"), (1, "negative_control")):
            mode = ttk.Combobox(
                controls, textvariable=self.vars[f"{key}_mode"],
                values=("None", "Metadata column"), state="readonly", width=18,
            )
            mode.grid(row=row, column=1, sticky="ew", padx=(0, 12), pady=5)
            column = ttk.Combobox(
                controls, textvariable=self.vars[f"{key}_column"], state="readonly"
            )
            column.grid(row=row, column=2, sticky="ew", padx=(0, 8), pady=5)
            level = ttk.Combobox(
                controls, textvariable=self.vars[f"{key}_level"], state="readonly"
            )
            level.grid(row=row, column=3, sticky="ew", pady=5)
            column.bind("<<ComboboxSelected>>", lambda _event, control=key: self.refresh_control_levels(control))
            mode.bind("<<ComboboxSelected>>", lambda _event, control=key: self.update_control_fields(control))
            self.control_widgets[key] = (mode, column, level)
            self.controls.extend((mode, column, level))

        taxonomy = ttk.LabelFrame(downstream, text="Taxonomy", padding=9)
        taxonomy.grid(row=0, column=0, sticky="ew", pady=(0, 8))
        ttk.Label(
            taxonomy,
            text="Taxonomy is assigned using the selected DADA2-compatible reference database and marker profile.",
            wraplength=800, justify="left",
        ).grid(row=0, column=0, sticky="w")

        rarefaction = ttk.LabelFrame(downstream, text="Rarefaction", padding=9)
        rarefaction.grid(row=1, column=0, sticky="ew", pady=(0, 8))
        rarefaction.columnconfigure(0, weight=1)
        self.rarefaction_var = self.vars["rarefaction_mode"]
        for row, (label, value) in enumerate((
            ("Automatic", "Automatic"),
            ("Minimum depth", "Minimum depth"),
            ("Custom depth", "Custom depth"),
        )):
            radio = ttk.Radiobutton(
                rarefaction, text=label, variable=self.rarefaction_var, value=value,
                command=self.update_rarefaction_fields,
            )
            radio.grid(row=row, column=0, sticky="w", pady=(2, 0))
            self.controls.append(radio)
        self.minimum_depth_description = ttk.Label(
            rarefaction,
            text="All samples will be rarefied to the smallest sequencing depth.",
            foreground="#536273",
        )
        self.minimum_depth_description.grid(
            row=3, column=0, sticky="w", padx=(24, 0), pady=(0, 6)
        )
        self.custom_depth_row = ttk.Frame(rarefaction)
        self.custom_depth_row.grid(row=4, column=0, sticky="w", padx=(24, 0), pady=(0, 8))
        ttk.Label(self.custom_depth_row, text="Rarefaction depth (reads)").pack(side="left", padx=(0, 8))
        self.custom_depth_entry = ttk.Entry(
            self.custom_depth_row, textvariable=self.vars["custom_rarefaction_depth"], width=12
        )
        self.custom_depth_entry.pack(side="left")
        self.controls.append(self.custom_depth_entry)
        self.rarefaction_advanced_toggle = ttk.Button(
            rarefaction,
            text="Advanced automatic settings  ▼",
            command=self.toggle_rarefaction_advanced,
        )
        self.rarefaction_advanced_toggle.grid(row=5, column=0, sticky="w", pady=(4, 0))
        self.controls.append(self.rarefaction_advanced_toggle)
        self.rarefaction_advanced = ttk.LabelFrame(rarefaction, padding=10)
        self.rarefaction_advanced.grid(row=6, column=0, sticky="ew", pady=(8, 0))
        self.rarefaction_advanced.grid_remove()
        for col in (1, 3):
            self.rarefaction_advanced.columnconfigure(col, weight=1)
        self.rarefaction_advanced_entries = []
        for col, label, key in (
            (0, "Maximum sample loss (%)", "max_sample_loss_percent"),
            (2, "Plateau gain threshold (%)", "plateau_gain_percent"),
            (0, "Minimum samples at plateau (%)", "plateau_samples_percent"),
            (2, "Evaluation step size (reads)", "read_step"),
        ):
            row = 0 if key in {"max_sample_loss_percent", "plateau_gain_percent"} else 1
            ttk.Label(self.rarefaction_advanced, text=label).grid(
                row=row, column=col, sticky="w", padx=(0, 8), pady=4
            )
            entry = ttk.Entry(self.rarefaction_advanced, textvariable=self.vars[key], width=10)
            entry.grid(row=row, column=col + 1, sticky="ew", padx=(0, 18), pady=4)
            self.rarefaction_advanced_entries.append(entry)
            self.controls.append(entry)

        alpha = ttk.LabelFrame(downstream, text="Alpha diversity", padding=9)
        alpha.grid(row=2, column=0, sticky="ew", pady=(0, 8))
        ttk.Label(alpha, text="Metrics", font=("Segoe UI", 10, "bold")).grid(
            row=0, column=0, sticky="w", pady=(0, 5)
        )
        for row, (metric, variable) in enumerate(self.alpha_vars.items(), start=1):
            check = ttk.Checkbutton(
                alpha, text="Pielou's evenness" if metric == "Pielou" else metric, variable=variable
            )
            check.grid(row=row, column=0, sticky="w", pady=2)
            self.controls.append(check)
        alpha_statistics = ttk.Checkbutton(
            alpha, text="Run statistical comparisons", variable=self.run_alpha_statistics
        )
        alpha_statistics.grid(row=0, column=1, sticky="w", padx=(36, 0))
        self.controls.append(alpha_statistics)
        ttk.Label(
            alpha,
            text="Automatic ANOVA + Tukey or Kruskal–Wallis + Dunn (BH)\n"
                 "Pairwise brackets for 2–3 groups; overall p-value for ≥4\n"
                 "Detailed statistics exported to tables",
            justify="left", wraplength=560,
        ).grid(row=1, column=1, rowspan=5, sticky="nw", padx=(36, 0), pady=(6, 0))

        beta = ttk.LabelFrame(downstream, text="Beta diversity", padding=9)
        beta.grid(row=3, column=0, sticky="ew", pady=(0, 8))
        ttk.Label(beta, text="Metrics", font=("Segoe UI", 10, "bold")).grid(
            row=0, column=0, sticky="w", pady=(0, 5)
        )
        for row, (label, key) in enumerate((
            ("Bray-Curtis", "bray"), ("Weighted UniFrac", "wunifrac"),
            ("Unweighted UniFrac", "unifrac"), ("Jaccard", "jaccard"),
        ), start=1):
            check = ttk.Checkbutton(beta, text=label, variable=self.beta_vars[key])
            check.grid(row=row, column=0, sticky="w", pady=2)
            self.controls.append(check)
        ttk.Label(beta, text="Ordination", font=("Segoe UI", 10, "bold")).grid(
            row=0, column=1, sticky="w", padx=(36, 0), pady=(0, 5)
        )
        self.ordination_var = tk.StringVar(value=self.vars["ordination"].get())
        for row, method in enumerate(("PCoA", "NMDS", "PCA"), start=1):
            radio = ttk.Radiobutton(beta, text=method, variable=self.ordination_var, value=method)
            radio.grid(row=row, column=1, sticky="w", padx=(36, 0), pady=2)
            self.controls.append(radio)

        ttk.Label(beta, text="Statistical comparisons", font=("Segoe UI", 10, "bold")).grid(
            row=5, column=0, columnspan=2, sticky="w", pady=(12, 5)
        )
        statistics = ttk.Frame(beta)
        statistics.grid(row=6, column=0, columnspan=2, sticky="ew")
        permanova = ttk.Checkbutton(
            statistics, text="Run PERMANOVA", variable=self.run_permanova
        )
        permanova.grid(row=0, column=0, sticky="w", pady=(0, 6))
        dispersion = ttk.Checkbutton(
            statistics, text="Run betadisper before PERMANOVA", variable=self.run_betadisper
        )
        dispersion.grid(row=1, column=0, sticky="w")
        ttk.Label(statistics, text="Permutations").grid(row=2, column=0, sticky="w", pady=(8, 0))
        permutations = ttk.Entry(statistics, textvariable=self.vars["permutations"], width=12)
        permutations.grid(row=2, column=1, sticky="w", padx=(8, 0), pady=(8, 0))
        self.controls.extend((permanova, dispersion, permutations))

        differential = ttk.LabelFrame(downstream, text="Differential abundance method", padding=9)
        differential.grid(row=4, column=0, sticky="ew", pady=(0, 8))
        self.differential_var = tk.StringVar(value=self.vars["differential_method"].get())
        descriptions = {
            "ANCOM-BC": "",
            "DESeq2": "Uses a negative-binomial model to compare counts between groups.",
            "ALDEx2": "Uses compositional analysis with Monte Carlo sampling.",
        }
        for row, method in enumerate(("ANCOM-BC", "DESeq2", "ALDEx2")):
            radio = ttk.Radiobutton(
                differential, text=method, variable=self.differential_var, value=method,
                command=self.update_differential_description,
            )
            radio.grid(row=row, column=0, sticky="w", pady=2)
            self.controls.append(radio)
        self.differential_description = ttk.Label(
            differential, text="", foreground="#536273", wraplength=760, justify="left"
        )
        self.differential_description.grid(row=3, column=0, sticky="w", pady=(7, 0))
        self.differential_descriptions = descriptions
        abundance = ttk.LabelFrame(downstream, text="Relative abundance visualization", padding=9)
        abundance.grid(row=5, column=0, sticky="ew", pady=(0, 8))
        self.abundance_vars = {
            key: tk.BooleanVar(value=self.vars["abundance_plot_type"].get() in {key, "both"})
            for key in ("barplot", "heatmap")
        }
        for column, (label, value) in enumerate((
            ("Stacked bar plots", "barplot"), ("Abundance heatmap", "heatmap"),
        )):
            check = ttk.Checkbutton(
                abundance, text=label, variable=self.abundance_vars[value],
                command=self.update_abundance_selection,
            )
            check.grid(row=0, column=column, sticky="w", padx=(0, 24), pady=(0, 6))
            self.controls.append(check)
        self.barplot_fields = ttk.Frame(abundance)
        self.barplot_fields.grid(row=1, column=0, columnspan=2, sticky="ew")
        self.heatmap_fields = ttk.Frame(abundance)
        self.heatmap_fields.grid(row=2, column=0, columnspan=2, sticky="ew", pady=(6, 0))
        ttk.Label(self.heatmap_fields, text="Taxonomic rank").grid(row=0, column=0, sticky="w", pady=4)
        self.heatmap_rank_combo = ttk.Combobox(
            self.heatmap_fields, textvariable=self.vars["heatmap_rank"], state="readonly", width=16
        )
        self.heatmap_rank_combo.grid(row=0, column=1, sticky="w", padx=(8, 0), pady=4)
        ttk.Label(self.heatmap_fields, text="Top N taxa").grid(row=1, column=0, sticky="w", pady=4)
        self.heatmap_top_n_spinbox = ttk.Spinbox(
            self.heatmap_fields, from_=10, to=100, increment=1,
            textvariable=self.vars["heatmap_top_n"], width=10,
        )
        self.heatmap_top_n_spinbox.grid(row=1, column=1, sticky="w", padx=(8, 0), pady=4)
        self.controls.extend((self.heatmap_rank_combo, self.heatmap_top_n_spinbox))
        self.abundance_rank_labels = []
        for row, key in (
            (0, "min_rank1_abundance_percent"),
            (1, "min_rank2_abundance_percent"),
        ):
            label = ttk.Label(self.barplot_fields)
            label.grid(row=row, column=0, sticky="w", pady=4)
            self.abundance_rank_labels.append(label)
            entry = ttk.Entry(self.barplot_fields, textvariable=self.vars[key], width=10)
            entry.grid(row=row, column=1, sticky="w", padx=(8, 4), pady=4)
            ttk.Label(self.barplot_fields, text='% into "Other"').grid(row=row, column=2, sticky="w", pady=4)
            self.controls.append(entry)
        ttk.Label(self.barplot_fields, text='Low-abundance taxa are grouped into "Other".').grid(
            row=2, column=0, columnspan=3, sticky="w", pady=(4, 0)
        )
        self.update_marker_profile()

        self.abundance_warning_label = ttk.Label(
            self.barplot_fields, text="", foreground="#965300", wraplength=760, justify="left"
        )
        self.abundance_warning_label.grid(row=3, column=0, columnspan=3, sticky="w", pady=(6, 0))
        for key in ("min_rank1_abundance_percent", "min_rank2_abundance_percent"):
            self.vars[key].trace_add("write", self.update_abundance_warning)
        self.update_abundance_warning()
        self.update_abundance_fields()

        plot_output = ttk.LabelFrame(downstream, text="Plot output", padding=9)
        plot_output.grid(row=6, column=0, sticky="ew", pady=(0, 8))
        ttk.Label(plot_output, text="Figure format").grid(
            row=0, column=0, sticky="w", pady=4
        )
        self.figure_format_var = self.vars["figure_format"]
        self.figure_format_var.set(self.figure_format_var.get().upper())
        for column, figure_format in enumerate(("PNG", "PDF"), start=1):
            radio = ttk.Radiobutton(
                plot_output,
                text=figure_format,
                variable=self.figure_format_var,
                value=figure_format,
            )
            radio.grid(row=0, column=column, sticky="w", padx=(8, 4), pady=4)
            self.controls.append(radio)

        self.console_frame = ttk.LabelFrame(self.main_frame, text="Progress", padding=(8, 4))
        self.console_frame.grid(row=2, column=0, sticky="ew", pady=(8, 0))
        self.console_frame.columnconfigure(0, weight=1)
        self.log = scrolledtext.ScrolledText(
            self.console_frame, height=3, wrap="word", state="disabled"
        )
        self.log.grid(row=0, column=0, sticky="ew")

        action_bar = ttk.Frame(self.main_frame, padding=(0, 6, 0, 0))
        action_bar.grid(row=3, column=0, sticky="ew")
        action_bar.columnconfigure(1, weight=1)
        self.run_button = self.button(action_bar, "Run selected steps", self.start)
        self.run_button.grid(row=0, column=0, sticky="w")
        self.status.set("Ready")
        ttk.Label(action_bar, textvariable=self.status).grid(row=0, column=1, sticky="w", padx=14)
        self.stop_button = ttk.Button(action_bar, text="Stop", command=self.stop, state="disabled")
        self.stop_button.grid(row=0, column=2, sticky="e")
        self.controls.append(self.stop_button)

        self.update_primer_fields()
        self.update_truncation_fields()
        self.update_pipeline_steps()
        self.update_rarefaction_fields()
        self.update_merging_fields()
        self.update_differential_description()
        for key in self.control_widgets:
            self.update_control_fields(key)

    def update_abundance_selection(self):
        selected = [key for key, variable in self.abundance_vars.items() if variable.get()]
        self.vars["abundance_plot_type"].set("both" if len(selected) == 2 else (
            selected[0] if selected else "none"
        ))
        self.update_abundance_fields()

    def update_abundance_fields(self):
        plot_type = self.vars["abundance_plot_type"].get()
        if plot_type in {"heatmap", "both"}:
            self.heatmap_fields.grid()
        else:
            self.heatmap_fields.grid_remove()
        if plot_type in {"barplot", "both"}:
            self.barplot_fields.grid()
        else:
            self.barplot_fields.grid_remove()

    def update_abundance_warning(self, *_args):
        zero_threshold = False
        for key in ("min_rank1_abundance_percent", "min_rank2_abundance_percent"):
            try:
                zero_threshold = zero_threshold or float(self.vars[key].get()) == 0
            except ValueError:
                pass
        if zero_threshold:
            self.abundance_warning_label.configure(
                text="0% keeps all assigned taxa. Plots can become crowded. "
                     "Increase the collapse percentage to group low-abundance taxa into Other."
            )
            self.abundance_warning_label.grid()
        else:
            self.abundance_warning_label.grid_remove()

    def set_running_controls(self, running):
        if running:
            # Include all interactive widgets, including radio buttons that were
            # not explicitly registered when the interface was constructed.
            controls = list(self.controls)
            pending = list(self.main_frame.winfo_children())
            while pending:
                widget = pending.pop()
                pending.extend(widget.winfo_children())
                if isinstance(widget, (ttk.Button, ttk.Checkbutton, ttk.Radiobutton,
                                       ttk.Entry, ttk.Combobox, ttk.Spinbox, ttk.Scale)):
                    if widget not in controls:
                        controls.append(widget)
            # Snapshot immediately before each run, after user changes.
            self.enabled_states = [(widget, widget.cget("state")) for widget in controls]
            for widget, _state in self.enabled_states:
                widget.configure(state="disabled")
            self.stop_button.configure(state="normal")
        else:
            for widget, state in self.enabled_states:
                widget.configure(state=state)
            self.stop_button.configure(state="disabled")

    def scrollable_tab(self, title):
        tab = ttk.Frame(self.notebook)
        tab.columnconfigure(0, weight=1)
        tab.rowconfigure(0, weight=1)
        canvas = tk.Canvas(tab, highlightthickness=0, borderwidth=0)
        scrollbar = ttk.Scrollbar(tab, orient="vertical", command=canvas.yview)
        canvas.configure(yscrollcommand=scrollbar.set)
        canvas.grid(row=0, column=0, sticky="nsew")
        scrollbar.grid(row=0, column=1, sticky="ns")
        content = ttk.Frame(canvas, padding=10)
        content.columnconfigure(0, weight=1)
        window = canvas.create_window((0, 0), window=content, anchor="nw")
        content.bind(
            "<Configure>", lambda _event: canvas.configure(scrollregion=canvas.bbox("all"))
        )
        canvas.bind(
            "<Configure>", lambda event: canvas.itemconfigure(window, width=event.width)
        )
        self.notebook.add(tab, text=title)
        return content

    def toggle_advanced_filtering(self):
        self.advanced_open = not self.advanced_open
        if self.advanced_open:
            self.advanced_filtering.grid()
            self.advanced_toggle.configure(text="Advanced filtering  -")
        else:
            self.advanced_filtering.grid_remove()
            self.advanced_toggle.configure(text="Advanced filtering  +")

    def toggle_rarefaction_advanced(self):
        self.rarefaction_advanced_open = not self.rarefaction_advanced_open
        self.update_rarefaction_fields()

    def update_pipeline_steps(self):
        if self.step_vars["step_3"].get() and not self.step_vars["step_2"].get():
            self.step_vars["step_2"].set(True)
        if self.step_vars["step_2"].get() and not self.step_vars["step_1"].get():
            self.step_vars["step_1"].set(True)

        if hasattr(self, "step_checkboxes"):
            requires_processing = self.step_vars["step_2"].get() or self.step_vars["step_3"].get()
            requires_taxonomy = self.step_vars["step_3"].get()
            self.step_checkboxes["step_1"].configure(
                state="disabled" if requires_processing else "normal"
            )
            self.step_checkboxes["step_2"].configure(
                state="disabled" if requires_taxonomy else "normal"
            )

    def update_primer_fields(self, _event=None):
        if self.vars["primer_mode"].get() == "sequences":
            self.sequence_fields.grid()
            self.length_fields.grid_remove()
        else:
            self.sequence_fields.grid_remove()
            self.length_fields.grid()

    def update_truncation_fields(self):
        automatic = self.vars["truncation_mode"].get() == "automatic"
        if automatic:
            self.trunc_qmin_entry.configure(state="normal")
            self.manual_truncation_fields.grid_remove()
        else:
            self.trunc_qmin_entry.configure(state="disabled")
            self.manual_truncation_fields.grid()

    def update_rarefaction_fields(self):
        mode = self.rarefaction_var.get()
        automatic = mode == "Automatic"
        if mode == "Custom depth":
            self.custom_depth_row.grid()
        else:
            self.custom_depth_row.grid_remove()
        if mode == "Minimum depth":
            self.minimum_depth_description.grid()
        else:
            self.minimum_depth_description.grid_remove()
        if automatic:
            self.rarefaction_advanced_toggle.grid()
        else:
            self.rarefaction_advanced_toggle.grid_remove()
        for entry in self.rarefaction_advanced_entries:
            entry.configure(state="normal" if automatic else "disabled")
        if automatic and self.rarefaction_advanced_open:
            self.rarefaction_advanced.grid()
            self.rarefaction_advanced_toggle.configure(
                text="Advanced automatic settings  ▲"
            )
        else:
            self.rarefaction_advanced.grid_remove()
            self.rarefaction_advanced_toggle.configure(
                text="Advanced automatic settings  ▼"
            )

    def update_merging_fields(self):
        if self.vars["mergepairs_strategy"].get() == "Merge paired reads":
            self.merge_fields.grid()
        else:
            self.merge_fields.grid_remove()

    def update_control_fields(self, control):
        mode, column, level = self.control_widgets[control]
        enabled = self.vars[f"{control}_mode"].get() == "Metadata column"
        state = "readonly" if enabled else "disabled"
        column.configure(state=state)
        level.configure(state=state)
        self.refresh_control_levels(control)

    def refresh_control_levels(self, control):
        if not hasattr(self, "control_widgets"):
            return
        column_widget = self.control_widgets[control][1]
        column_name = self.vars[f"{control}_column"].get()
        levels = []
        metadata_file = self.vars.get("metadata_file")
        if metadata_file and metadata_file.get() and column_name:
            try:
                _, rows = read_metadata(metadata_file.get())
                levels = sorted({row.get(column_name, "") for row in rows} - {"", "NA"})
            except (OSError, ValueError, UnicodeError):
                levels = []
        column_widget.configure(
            values=[column for column in self.metadata_columns if column != "SampleID"]
        )
        level_widget = self.control_widgets[control][2]
        level_widget.configure(values=levels)
        if self.vars[f"{control}_level"].get() not in levels:
            self.vars[f"{control}_level"].set(levels[0] if levels else "")

    def update_differential_description(self):
        description = self.differential_descriptions[self.differential_var.get()]
        self.differential_description.configure(text=description)
        if description:
            self.differential_description.grid()
        else:
            self.differential_description.grid_remove()

    def button(self, parent, text, command):
        widget = ttk.Button(parent, text=text, command=command)
        self.controls.append(widget)
        return widget

    def open_setup(self):
        if self.busy:
            messagebox.showinfo(self.text("running"), self.text("stop_before_close"), parent=self.root)
            return None
        self.setup_wizard = SetupWizard(
            self.root,
            BASE,
            self.vars["rscript"].get() or find_rscript(),
            find_cutadapt(BASE, self.saved.get("cutadapt", "")),
            self.setup_completed,
        )
        return self.setup_wizard

    def setup_completed(self, rscript, cutadapt):
        self.vars["rscript"].set(str(Path(rscript)))
        self.vars["cutadapt"].set(cutadapt)
        self.saved.update({
            "rscript": rscript,
            "cutadapt": cutadapt,
            "setup_completed": True,
        })
        try:
            save_settings(self.saved)
        except OSError as exc:
            messagebox.showerror(
                self.text("cannot_run"), self.text("settings_error", error=exc), parent=self.root
            )
        self.root.deiconify()

    def path_row(self, parent, row, label, key, types):
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w", padx=(0, 10))
        entry = ttk.Entry(parent, textvariable=self.vars[key])
        entry.grid(row=row, column=1, sticky="ew", pady=4)
        self.controls.append(entry)
        self.button(parent, self.text("choose"), partial(self.pick_path, key, types)).grid(
            row=row, column=2, sticky="ew", padx=(8, 0)
        )

    def option(self, parent, row, col, label, key, choices=None):
        ttk.Label(parent, text=label).grid(row=row, column=col, sticky="w", padx=5, pady=3)
        widget = (ttk.Entry(parent, textvariable=self.vars[key]) if choices is None else
                  ttk.Combobox(parent, textvariable=self.vars[key], values=choices))
        widget.grid(row=row, column=col + 1, sticky="ew", padx=5, pady=3)
        self.controls.append(widget)
        return widget

    def pick_path(self, key, types):
        try:
            if key == "fastq_dir" and types is None:
                path = self.choose_fastq_input()
            elif types is None:
                path = filedialog.askdirectory(
                    parent=self.root,
                    title=self.text("browse_output"),
                )
            else:
                path = filedialog.askopenfilename(
                    parent=self.root,
                    title=self.text("browse_key", label=self.text(key)),
                    filetypes=types,
                )
        except (OSError, tk.TclError) as exc:
            messagebox.showerror(
                self.text("metadata_title"),
                self.text("dialog_error", error=exc),
                parent=self.root,
            )
            return
        if path:
            self.vars[key].set(str(Path(path)))
            if key == "metadata_file":
                self.refresh_metadata()
            elif key == "fastq_dir":
                self.refresh_fastq_label()

    def choose_fastq_input(self):
        current_path = self.vars["fastq_dir"].get().strip()
        current = Path(current_path).expanduser() if current_path else Path.home()
        if current.is_file():
            current = current.parent
        if not current.is_dir():
            current = Path.home()
        initial_dir = str(current.resolve())

        if sys.platform == "win32":
            return _windows_fastq_input_dialog(
                self.root, "Select FASTQ folder", initial_dir
            )

        return filedialog.askdirectory(
            parent=self.root,
            title="Select FASTQ folder",
            initialdir=initial_dir,
            mustexist=True,
        )

    def select_fastq(self):
        files = filedialog.askopenfilenames(parent=self.root, title=self.text("browse_fastq"),
                                            filetypes=[("FASTQ / ZIP", "*.fastq *.fq *.fastq.gz *.fq.gz *.zip"), (self.text("all_files"), "*")])
        if files:
            self.fastq_files = list(files)
            self.refresh_fastq_label()

    def select_fastq_folder(self):
        folder = filedialog.askdirectory(parent=self.root, title=self.text("browse_fastq"))
        if not folder:
            return
        files = []
        excluded = {"result", "results", "runs", "filtered", "primer_trimmed", ".git"}
        output_text = self.vars["output_dir"].get().strip()
        output = Path(output_text).resolve() if output_text else None
        for directory, dirs, names in os.walk(folder):
            dirs[:] = [name for name in dirs if name not in excluded and
                       (output is None or (Path(directory) / name).resolve() != output)]
            files.extend(str(Path(directory) / name) for name in names
                         if name.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz")))
        self.fastq_files = sorted(files)
        self.refresh_fastq_label()
        metadata_files = sorted(Path(folder).glob("*.csv"))
        if len(metadata_files) == 1 and not self.vars["metadata_file"].get():
            self.vars["metadata_file"].set(str(metadata_files[0]))
            self.refresh_metadata()

    def refresh_fastq_label(self):
        folder = self.vars.get("fastq_dir")
        if folder and folder.get():
            input_path = folder.get().strip()
            self.fastq_label.configure(text="Counting FASTQ files…")
            threading.Thread(
                target=self._count_fastq_label,
                args=(input_path,),
                daemon=True,
            ).start()

    def _count_fastq_label(self, input_path):
        try:
            count = count_fastq_input_files(input_path)
            error = ""
        except (OSError, ValueError) as exc:
            count = 0
            error = str(exc)
        self.events.put(("fastq_count", (input_path, count, error)))

    def refresh_metadata(self, show_errors=True):
        try:
            columns, _ = read_metadata(self.vars["metadata_file"].get())
            self.metadata_columns = columns
            group_columns = [column for column in columns if column != "SampleID"]
            self.group_combo.configure(values=group_columns)
            if self.vars["group_column"].get() not in group_columns:
                self.vars["group_column"].set("")
            for control in self.control_widgets:
                self.control_widgets[control][1].configure(values=group_columns)
                if self.vars[f"{control}_column"].get() not in group_columns:
                    self.vars[f"{control}_column"].set("")
                self.refresh_control_levels(control)
            self.refresh_groups()
        except (OSError, ValueError, UnicodeError) as exc:
            if show_errors:
                messagebox.showerror(self.text("metadata_title"), str(exc), parent=self.root)

    def refresh_groups(self):
        try:
            _, rows = read_metadata(self.vars["metadata_file"].get())
            groups = sorted({row.get(self.vars["group_column"].get(), "") for row in rows} - {""})
            self.ref_combo.configure(values=groups)
            if self.vars["reference_group"].get() not in groups:
                self.vars["reference_group"].set("")
            for control in self.control_widgets:
                self.refresh_control_levels(control)
        except (OSError, ValueError, UnicodeError):
            pass

    def configuration(self):
        self.update_pipeline_steps()
        cfg = {key: var.get().strip() for key, var in self.vars.items()}
        for current, legacy in (
            ("reference_db_path", "silva_ref_path"),
            ("min_rank1_abundance_percent", "min_phylum_abundance_percent"),
            ("min_rank2_abundance_percent", "min_genus_abundance_percent"),
        ):
            if current not in cfg and legacy in cfg:
                cfg[current] = cfg.pop(legacy)
        required_settings = ("marker", "reference_db_path",
                             "min_rank1_abundance_percent", "min_rank2_abundance_percent")
        missing = [key for key in required_settings if key not in cfg]
        if missing:
            raise ValueError("Missing configuration settings: " + ", ".join(missing))
        cfg["marker"] = cfg["marker"].upper()
        profile = get_marker_profile(cfg["marker"])
        cfg.setdefault("abundance_plot_type", "barplot")
        cfg.setdefault("heatmap_rank", "Genus")
        cfg.setdefault("heatmap_top_n", "30")
        if cfg["abundance_plot_type"] not in {"barplot", "heatmap", "both"}:
            raise ValueError("Choose stacked bar plots, an abundance heatmap, or both.")
        if cfg["abundance_plot_type"] in {"heatmap", "both"}:
            rank1 = profile["relative_abundance_ranks"][0]
            if cfg["heatmap_rank"] in {"Phylum", "Division"}:
                cfg["heatmap_rank"] = rank1
            if cfg["heatmap_rank"] not in {"Genus", rank1}:
                raise ValueError("Choose Genus or the marker's Phylum/Division rank for the heatmap.")
            try:
                cfg["heatmap_top_n"] = int(cfg["heatmap_top_n"])
            except ValueError:
                raise ValueError("Heatmap Top N taxa must be an integer from 10 to 100.") from None
            if not 10 <= cfg["heatmap_top_n"] <= 100:
                raise ValueError("Heatmap Top N taxa must be an integer from 10 to 100.")
            if cfg["abundance_plot_type"] == "heatmap":
                # Hidden barplot controls do not constrain heatmap generation.
                cfg["min_rank1_abundance_percent"] = "1.0"
                cfg["min_rank2_abundance_percent"] = "2.0"
        else:
            cfg["heatmap_rank"] = "Genus"
            cfg["heatmap_top_n"] = 30
        cfg["asv_length_filter_enabled"] = self.asv_length_filter_enabled.get()
        if not any(variable.get() for variable in self.step_vars.values()):
            raise ValueError("Select at least one pipeline step.")
        required_inputs = {"fastq_dir": "FASTQ folder", "metadata_file": "Metadata file",
                           "output_dir": "Output folder"}
        if self.step_vars["step_2"].get() or self.step_vars["step_3"].get():
            required_inputs.update(reference_db_path="Reference database", group_column="Comparison group")
        if self.step_vars["step_3"].get():
            required_inputs["reference_group"] = "Reference group"
        if cfg["primer_mode"] == "sequences":
            required_inputs.update(primer_fwd="Forward primer", primer_rev="Reverse primer")
        else:
            required_inputs.update(trim_left_fwd="Forward primer length", trim_left_rev="Reverse primer length")
        if cfg["truncation_mode"] == "manual":
            required_inputs.update(trunc_len_fwd="Forward read length after primer removal",
                                   trunc_len_rev="Reverse read length after primer removal")
        if cfg["asv_length_filter_enabled"]:
            required_inputs.update(asv_min_len="Minimum ASV length", asv_max_len="Maximum ASV length")
        if self.step_vars["step_3"].get() and cfg["rarefaction_mode"] == "Custom depth":
            required_inputs["custom_rarefaction_depth"] = "Custom rarefaction depth"
        missing_inputs = [label for key, label in required_inputs.items() if not cfg.get(key)]
        if missing_inputs:
            raise ValueError("Fill in the required inputs:\n" + "\n".join(missing_inputs))
        # Inactive numeric fields are not user requirements; use zero sentinels.
        if cfg["primer_mode"] == "sequences":
            cfg["trim_left_fwd"] = cfg["trim_left_rev"] = "0"
        if cfg["truncation_mode"] != "manual":
            cfg["trunc_len_fwd"] = cfg["trunc_len_rev"] = "0"
        fastq_input = Path(cfg["fastq_dir"]) if cfg["fastq_dir"] else None
        if fastq_input is None or not (fastq_input.is_dir() or fastq_input.is_file()):
            raise ValueError("Select an existing FASTQ folder or ZIP archive.")
        output_dir = cfg["output_dir"]
        if not output_dir:
            raise ValueError(self.text("select_output"))
        cfg["output_dir"] = Path(output_dir).resolve().as_posix()
        fastq_dir, fastq_paths = resolve_fastq_input(fastq_input, cfg["output_dir"])
        cfg["fastq_dir"] = fastq_dir.as_posix()
        cfg["fastq_files"] = fastq_paths
        cfg["steps"] = {key: variable.get() for key, variable in self.step_vars.items()}

        primer_mode = cfg["primer_mode"]
        if primer_mode not in {"sequences", "length"}:
            raise ValueError("Choose how you want to provide primer information.")
        cfg["primer_method"] = "cutadapt" if primer_mode == "sequences" else "trim_left"
        cfg.pop("primer_mode")
        required_files = ["rscript", "metadata_file"]
        if self.step_vars["step_2"].get() or self.step_vars["step_3"].get():
            required_files.append("reference_db_path")
        if cfg["primer_method"] == "cutadapt":
            required_files.append("cutadapt")
        for key in required_files:
            if not cfg.get(key) or not Path(cfg[key]).is_file():
                if key == "cutadapt":
                    raise ValueError("Primer removal is unavailable. Run dependency setup and try again.")
                raise ValueError(self.text("invalid_file", key=key))
            cfg[key] = Path(cfg[key]).resolve().as_posix()
        if cfg["primer_method"] == "trim_left":
            cfg["cutadapt"] = ""
        pairs = pair_files(cfg["fastq_files"])
        columns, rows = read_metadata(cfg["metadata_file"])
        sample_id_map = match_metadata_sample_ids(pairs, (row["SampleID"] for row in rows))
        group = cfg["group_column"]
        needs_groups = self.step_vars["step_2"].get() or self.step_vars["step_3"].get()
        if needs_groups and (
            group not in columns
            or group == "SampleID"
            or not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", group)
        ):
            raise ValueError(self.text("invalid_group"))
        meta = {row["SampleID"]: row for row in rows}
        missing = {sample for sample in pairs if sample_id_map[sample] not in meta}
        if missing:
            raise ValueError(self.text("missing_metadata", samples=", ".join(sorted(missing)[:15])))

        control_samples = {}
        for control in ("mock_control", "negative_control"):
            mode = cfg[f"{control}_mode"]
            column = cfg[f"{control}_column"] if mode == "Metadata column" else ""
            level = cfg[f"{control}_level"] if column else ""
            if mode not in {"None", "Metadata column"}:
                raise ValueError(f"Choose a valid {control.replace('_', ' ')} option.")
            if column and (column not in columns or column == "SampleID" or not level):
                raise ValueError(f"Choose a metadata column and level for {control.replace('_', ' ')}.")
            if column and level not in {row[column] for row in rows}:
                raise ValueError(f"The selected {control.replace('_', ' ')} level is not in the metadata.")
            control_samples[control] = {
                sample: bool(column and meta[sample_id_map[sample]][column] == level)
                for sample in pairs
            }
            cfg[f"{control}_column"] = column
            cfg[f"{control}_level"] = level
            cfg.pop(f"{control}_mode")

        needs_diversity = self.step_vars["step_3"].get()
        is_negative = control_samples["negative_control"]
        is_mock = control_samples["mock_control"]
        run_decontam = bool(cfg["negative_control_column"]) and needs_diversity
        if needs_diversity and run_decontam and (
            not any(is_negative.values()) or all(is_negative.values())
        ):
            raise ValueError(self.text("invalid_decontam_controls"))
        if needs_diversity and cfg["mock_control_column"] and not any(is_mock.values()):
            raise ValueError("No selected FASTQ samples match the chosen mock-community level.")

        counts = {}
        if needs_groups:
            for sample in pairs:
                if needs_diversity and is_mock[sample]:
                    continue
                value = meta[sample_id_map[sample]][group]
                if not value or value != value.strip() or value == "NA":
                    raise ValueError(self.text("invalid_group_value", sample=sample))
                counts[value] = counts.get(value, 0) + 1
        if needs_diversity and (len(counts) < 2 or min(counts.values()) < 2):
            raise ValueError(self.text("minimum_groups"))
        if needs_diversity and cfg["reference_group"] not in counts:
            raise ValueError(self.text("invalid_reference"))

        for key in ("trunc_q", "min_len", "permutations"):
            try:
                cfg[key] = int(cfg[key])
            except ValueError:
                raise ValueError(self.text("integer", key=key)) from None
            if cfg[key] < 0 or (key != "trunc_q" and cfg[key] == 0):
                raise ValueError(self.text("positive", key=key))
        if cfg["asv_length_filter_enabled"]:
            for key in ("asv_min_len", "asv_max_len"):
                try:
                    cfg[key] = int(cfg[key])
                except ValueError:
                    raise ValueError(self.text("integer", key=key)) from None
                if cfg[key] <= 0:
                    raise ValueError(self.text("positive", key=key))
            if cfg["asv_min_len"] >= cfg["asv_max_len"]:
                raise ValueError("Minimum ASV length must be less than maximum ASV length.")
        else:
            cfg["asv_min_len"] = None
            cfg["asv_max_len"] = None
        for key in (
            "trim_left_fwd", "trim_left_rev", "max_n", "trunc_qmin",
            "trunc_len_fwd", "trunc_len_rev", "mergepairs_min_overlap",
            "mergepairs_max_mismatch", "read_step",
        ):
            try:
                cfg[key] = int(cfg[key])
            except ValueError:
                raise ValueError(self.text("integer", key=key)) from None
            if cfg[key] < 0 or (key in {"trunc_qmin", "mergepairs_min_overlap", "read_step"} and cfg[key] == 0):
                raise ValueError(self.text("nonnegative", key=key))
        for key in (
            "max_ee_fwd", "max_ee_rev", "max_sample_loss_percent",
            "plateau_gain_percent", "plateau_samples_percent",
            "min_rank1_abundance_percent", "min_rank2_abundance_percent",
        ):
            try:
                cfg[key] = float(cfg[key])
            except ValueError:
                raise ValueError(f"{key.replace('_', ' ')} must be a number.") from None
            if not math.isfinite(cfg[key]):
                raise ValueError(f"{key.replace('_', ' ')} must be a finite number.")
            if cfg[key] < 0:
                raise ValueError(f"{key.replace('_', ' ')} must be zero or greater.")
        if cfg["max_sample_loss_percent"] >= 100:
            raise ValueError("Maximum sample loss must be less than 100%.")
        if cfg["plateau_gain_percent"] > 100 or cfg["plateau_samples_percent"] > 100:
            raise ValueError("Rarefaction percentages must not exceed 100%.")
        if (
            cfg["min_rank1_abundance_percent"] > 100
            or cfg["min_rank2_abundance_percent"] > 100
        ):
            raise ValueError("Relative abundance thresholds must not exceed 100%.")
        cfg["figure_format"] = cfg["figure_format"].strip().lower()
        if cfg["figure_format"] not in {"png", "pdf"}:
            raise ValueError("Choose PNG or PDF for the figure format.")
        # Choose automatically, leaving one logical processor for other work.
        # Overwrite legacy thread selections rather than restoring a UI setting.
        cfg["threads"] = max(1, (os.cpu_count() or 1) - 1)

        cfg["max_ee"] = [cfg.pop("max_ee_fwd"), cfg.pop("max_ee_rev")]
        cfg["trunc_len"] = [cfg.pop("trunc_len_fwd"), cfg.pop("trunc_len_rev")]
        cfg["automatic_truncation"] = cfg.pop("truncation_mode") == "automatic"
        cfg["manual_read_length"] = None
        if not cfg["automatic_truncation"]:
            if any(length <= 0 for length in cfg["trunc_len"]):
                raise ValueError("Read lengths after primer removal must be positive integers.")
            cfg["manual_read_length"] = cfg["trunc_len"].copy()
            if cfg["primer_method"] == "trim_left":
                cfg["trunc_len"] = [
                    length + cfg[key] for length, key in zip(
                        cfg["trunc_len"], ("trim_left_fwd", "trim_left_rev")
                    )
                ]
        cfg["disable_truncation"] = False
        cfg["cutadapt_error_rate"] = 0.1
        cfg["discard_untrimmed"] = False
        cfg["max_sample_loss"] = cfg.pop("max_sample_loss_percent") / 100
        cfg["plateau_gain"] = cfg.pop("plateau_gain_percent") / 100
        cfg["plateau_samples"] = cfg.pop("plateau_samples_percent") / 100
        cfg["remove_phix"] = self.remove_phix.get()
        cfg["chimera_method"] = cfg["chimera_method"].casefold()
        cfg["mergepairs_strategy"] = (
            "merge" if cfg["mergepairs_strategy"] == "Merge paired reads" else "concatenate"
        )
        cfg["run_permanova"] = self.run_permanova.get()
        cfg["run_betadisper"] = self.run_betadisper.get()
        cfg["run_alpha_statistics"] = self.run_alpha_statistics.get()
        cfg["alpha_measures"] = [key for key, variable in self.alpha_vars.items() if variable.get()]
        if needs_diversity and not cfg["alpha_measures"]:
            raise ValueError("Select at least one alpha-diversity metric.")
        for key in ("enable_cache", "force_recompute", "cache_dir"):
            cfg.pop(key, None)
        cfg["beta_distances"] = [key for key, variable in self.beta_vars.items() if variable.get()]
        cfg["ordination"] = self.ordination_var.get()
        cfg["differential_method"] = self.differential_var.get()
        rarefaction_modes = {
            "Automatic": "auto_rarefy",
            "Minimum depth": "min",
            "Custom depth": "custom",
        }
        try:
            cfg["rarefaction_mode"] = rarefaction_modes.pop(cfg["rarefaction_mode"])
        except KeyError:
            raise ValueError(self.text("rarefaction")) from None
        if cfg["rarefaction_mode"] == "custom":
            try:
                cfg["custom_rarefaction_depth"] = int(cfg["custom_rarefaction_depth"])
            except ValueError:
                raise ValueError(self.text("integer", key="Depth")) from None
            if cfg["custom_rarefaction_depth"] <= 0:
                raise ValueError(self.text("positive", key="Depth"))
        try:
            cfg["permanova_permutations"] = int(cfg.pop("permutations"))
        except ValueError:
            raise ValueError("Permutations must be a positive integer.") from None
        if cfg["permanova_permutations"] < 1:
            raise ValueError("Permutations must be a positive integer.")

        cfg["pooling"] = cfg["pooling"] == "true" if cfg["pooling"] == "true" else (
            False if cfg["pooling"] == "false" else "pseudo"
        )
        if cfg["primer_method"] == "cutadapt":
            for key in ("primer_fwd", "primer_rev"):
                cfg[key] = cfg[key].upper()
                if not re.fullmatch(r"[ACGTRYSWKMBDHVN]+", cfg[key]):
                    raise ValueError(self.text("primer", key=key.replace("primer_", "").title()))
        cfg["fastq_files"] = [path for sample in sorted(pairs) for path in pairs[sample].values()]
        cfg["sample_id_map"] = sample_id_map
        cfg["negative_control_prefixes"] = []
        cfg["positive_control_prefixes"] = []
        cfg["run_decontam"] = run_decontam
        cfg["exclude_mock_controls"] = bool(cfg["mock_control_column"]) and needs_diversity
        cfg["exclude_controls"] = False
        cfg["project_dir"] = cfg["fastq_dir"]
        cfg["pipeline_dir"] = (BASE / "pipeline").as_posix()
        cfg["marker_profiles_path"] = MARKER_PROFILES_PATH.as_posix()
        cfg["run_ancombc"] = cfg["differential_method"] == "ANCOM-BC"
        if needs_diversity and cfg["run_permanova"] and not cfg["beta_distances"]:
            raise ValueError("Select at least one beta-diversity metric or turn off PERMANOVA.")
        scripts = sorted((BASE / "pipeline").glob("[0-9][0-9]_*.R"))
        if len(scripts) != 17 or [path.name[:2] for path in scripts] != [f"{i:02}" for i in range(17)]:
            raise ValueError(self.text("scripts"))
        return cfg

    def start(self, check_only=False):
        if self.busy:
            return
        try:
            cfg = self.configuration()
            prefix = "check" if check_only else "analysis"
            run_dir = Path(cfg["output_dir"]) / f"{prefix}_{datetime.now():%Y%m%d_%H%M%S_%f}"
            run_dir.mkdir(parents=True, exist_ok=False)
            for folder in ("rds", "tables", "figures", "logs"):
                (run_dir / folder).mkdir()
            self.temp_dir = create_app_temp_dir()
            cfg["run_dir"] = run_dir.as_posix()
            config_path = self.temp_dir / "biome_runtime.json"
            cfg["config_path"] = config_path.as_posix()
            config_path.write_text(json.dumps(cfg, ensure_ascii=False, indent=2), encoding="utf-8")
            write_configuration_report(cfg, run_dir / "logs" / "BIOME_config.txt")
            with (run_dir / "tables" / "input_path.csv").open("w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle)
                writer.writerow(["SampleID", "R1", "R2"])
                sample_id_map = cfg["sample_id_map"]
                for sample, mates in sorted(pair_files(cfg["fastq_files"]).items()):
                    writer.writerow([sample_id_map[sample], mates["1"], mates["2"]])
            try:
                saved = {**self.saved, **cfg, "setup_completed": True}
                saved.update({key: variable.get() for key, variable in self.vars.items()})
                saved.update({key: variable.get() for key, variable in self.step_vars.items()})
                saved.update({
                    "remove_phix": self.remove_phix.get(),
                    "run_permanova": self.run_permanova.get(),
                    "run_betadisper": self.run_betadisper.get(),
                    "asv_length_filter_enabled": self.asv_length_filter_enabled.get(),
                    "beta_bray": self.beta_vars["bray"].get(),
                    "beta_wunifrac": self.beta_vars["wunifrac"].get(),
                    "beta_unifrac": self.beta_vars["unifrac"].get(),
                    "beta_jaccard": self.beta_vars["jaccard"].get(),
                    "rarefaction_mode": self.rarefaction_var.get(),
                    "ordination": self.ordination_var.get(),
                    "differential_method": self.differential_var.get(),
                })
                save_settings(saved)
                self.saved = saved
            except OSError as exc:
                messagebox.showerror(
                    self.text("cannot_run"), self.text("settings_error", error=exc), parent=self.root
                )
        except (OSError, ValueError, UnicodeError) as exc:
            error = str(exc)
            if self.temp_dir is not None:
                try:
                    _remove_owned_temp_dir(self.temp_dir, root=self.temp_dir.parent)
                except (OSError, ValueError) as cleanup_exc:
                    error += f"\nTemporary cache cleanup failed: {cleanup_exc}"
                self.temp_dir = None
            messagebox.showerror(self.text("cannot_run"), error, parent=self.root)
            return
        self.last_run = run_dir
        self.busy = True
        self.cancelled.clear()
        self.set_running_controls(True)
        self.status.set("Running")
        self.append_log("\n" + self.text("run_folder", path=run_dir) + "\n")
        threading.Thread(target=self.worker, args=(cfg, check_only), daemon=True).start()

    def worker(self, cfg, check_only):
        run_dir = Path(cfg["run_dir"])
        log_dir = run_dir / "logs"
        log_dir.mkdir(parents=True, exist_ok=True)
        command = [cfg["rscript"], "--no-save", "--no-restore", "--encoding=UTF-8",
                   str(BASE / "run_pipeline.R"), cfg.get("config_path", str(run_dir / "rds" / "biome_runtime.json"))]
        if check_only:
            command.append("--check")
        try:
            env = os.environ.copy()
            temp_dir = getattr(self, "temp_dir", None)
            if temp_dir is not None:
                env.update({
                    "TMPDIR": str(temp_dir),
                    "TMP": str(temp_dir),
                    "TEMP": str(temp_dir),
                })
            with (log_dir / "pipeline.log").open("w", encoding="utf-8") as log:
                start_message = f"[{datetime.now().astimezone().isoformat(timespec='seconds')}] Starting R pipeline.\n"
                log.write(start_message)
                log.flush()
                self.events.put(("log", start_message))
                process = subprocess.Popen(command, cwd=str(BASE), stdout=subprocess.PIPE,
                                           stderr=subprocess.STDOUT, text=True, encoding="utf-8",
                                           errors="replace", env=env,
                                           creationflags=(subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0))
                self.process = process
                if self.cancelled.is_set():
                    self.kill_process(process)
                plot_warning_sent = False
                for line in process.stdout:
                    log.write(line)
                    log.flush()
                    self.events.put(("log", line))
                    if line.startswith("BIOME_PLOT_WARNING:") and not plot_warning_sent:
                        self.events.put(("plot_warning", line.partition(":")[2].strip()))
                        plot_warning_sent = True
                code = process.wait()
            result = {"exit_code": code, "cancelled": self.cancelled.is_set(), "check_only": check_only}
            (log_dir / "status.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
        except Exception as exc:
            error = "\n" + self.text("launch_error", error=exc) + "\n"
            result = {"exit_code": -1, "cancelled": self.cancelled.is_set(), "check_only": check_only}
            try:
                with (log_dir / "pipeline.log").open("a", encoding="utf-8") as log:
                    log.write(error)
                (log_dir / "status.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
            except OSError:
                pass
            self.events.put(("log", error))
        finally:
            self.process = None
            temp_dir = getattr(self, "temp_dir", None)
            if temp_dir is not None:
                try:
                    _remove_owned_temp_dir(temp_dir, root=temp_dir.parent)
                except (OSError, ValueError) as exc:
                    self.events.put(("cleanup_error", str(exc)))
                finally:
                    self.temp_dir = None
        self.events.put(("done", result))

    def poll(self):
        try:
            while True:
                kind, data = self.events.get_nowait()
                if kind == "log":
                    self.append_log(data)
                elif kind == "plot_warning":
                    messagebox.showwarning("Relative abundance plots", data, parent=self.root)
                elif kind == "fastq_count":
                    input_path, count, error = data
                    if self.vars["fastq_dir"].get().strip() == input_path:
                        label = "FASTQ input unavailable" if error else f"{count} FASTQ files"
                        self.fastq_label.configure(text=label)
                elif kind == "done":
                    self.busy = False
                    self.set_running_controls(False)
                    if data["cancelled"]:
                        self.status.set("Stopped")
                    elif data["exit_code"] == 0:
                        self.status.set("Completed")
                        messagebox.showinfo(self.text("complete"), self.status.get() + f"\n{self.last_run}", parent=self.root)
                    else:
                        self.status.set(f"Failed (code {data['exit_code']})")
                        messagebox.showerror(self.text("pipeline_stopped"), self.text("log_at", path=self.last_run), parent=self.root)
                elif kind == "cleanup_error":
                    messagebox.showerror(
                        self.text("cannot_run"),
                        self.text("temp_cache_error", error=data),
                        parent=self.root,
                    )
        except queue.Empty:
            pass
        self.root.after(100, self.poll)

    def append_log(self, text):
        self.log.configure(state="normal")
        self.log.insert("end", text)
        # Keep the interface responsive during verbose, long analyses.
        if int(self.log.index("end-1c").split(".")[0]) > 4000:
            self.log.delete("1.0", "1000.0")
        self.log.see("end")
        self.log.configure(state="disabled")

    @staticmethod
    def kill_process(process):
        if process.poll() is not None:
            return
        if os.name == "nt":
            subprocess.run(["taskkill", "/PID", str(process.pid), "/T", "/F"],
                           capture_output=True, creationflags=subprocess.CREATE_NO_WINDOW)
        else:
            process.terminate()

    def stop(self):
        self.cancelled.set()
        self.stop_button.configure(state="disabled")
        self.set_status("stopping")
        if self.process is not None:
            threading.Thread(target=self.kill_process, args=(self.process,), daemon=True).start()

    def open_results(self):
        path = self.last_run or Path(self.vars["output_dir"].get())
        if Path(path).is_dir() and os.name == "nt":
            os.startfile(str(path))
        else:
            messagebox.showinfo(self.text("results"), str(path), parent=self.root)

    def close(self):
        if self.busy:
            messagebox.showinfo(self.text("running"), self.text("stop_before_close"), parent=self.root)
            return
        try:
            clear_orphaned_temp_dirs()
        except (OSError, ValueError) as exc:
            messagebox.showerror(
                self.text("cannot_run"), self.text("temp_cache_error", error=exc), parent=self.root
            )
        self.root.destroy()

    def clear_temporary_cache(self):
        if self.busy:
            messagebox.showinfo(
                self.text("running"), self.text("temp_cache_running"), parent=self.root
            )
            return
        try:
            freed = clear_orphaned_temp_dirs()
        except (OSError, ValueError) as exc:
            messagebox.showerror(
                self.text("cannot_run"), self.text("temp_cache_error", error=exc), parent=self.root
            )
            return
        messagebox.showinfo(
            self.text("complete"),
            self.text("temp_cache_cleared", size=format_size(freed)),
            parent=self.root,
        )


def main():
    root = None
    try:
        root = tk.Tk()
        root.report_callback_exception = (
            lambda exc_type, exc_value, exc_traceback:
                report_startup_failure(exc_type, exc_value, exc_traceback, parent=root)
        )
        root.withdraw()
        app = PipelineApp(root)
        app.open_setup()
        root.mainloop()
    except Exception:
        report_startup_failure(*sys.exc_info(), parent=root)
    finally:
        if root is not None:
            try:
                root.destroy()
            except tk.TclError:
                pass


if __name__ == "__main__":
    main()
