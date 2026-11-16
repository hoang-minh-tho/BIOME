"""Runtime dependency discovery and the first-run setup wizard."""

import os
from pathlib import Path
import queue
import re
import shutil
import subprocess
import tempfile
import threading
import tkinter as tk
from tkinter import messagebox, ttk
import urllib.error
import urllib.parse
import urllib.request


WIZARD_TEXT = {
    "title": "Preparing BIOME",
    "loading": "Loading...",
    "progress": "Preparing setup: {percent}%",
    "prompt_title": "Additional components required",
    "prompt_message": "BIOME needs to install required components to run.\n\nAllow installation now?",
    "install": "Install",
    "cancel": "Cancel",
    "failure_title": "BIOME could not start",
    "failure_message": "Required components could not be installed. Please try again.\n\nDetails:\n",
}


def _default_r_roots():
    roots = [
        Path(os.environ.get("ProgramFiles", "C:/Program Files")) / "R",
        Path(os.environ.get("ProgramFiles(x86)", "C:/Program Files (x86)")) / "R",
        Path(os.environ.get("LOCALAPPDATA", Path.home() / "AppData/Local")) / "Programs/R",
        Path("C:/R"),
    ]
    return roots


def _registry_r_roots():
    if os.name != "nt":
        return []
    try:
        import winreg
    except ImportError:
        return []
    roots = []
    for hive in (winreg.HKEY_CURRENT_USER, winreg.HKEY_LOCAL_MACHINE):
        for registry_path in (r"SOFTWARE\R-core\R64", r"SOFTWARE\R-core\R"):
            try:
                with winreg.OpenKey(hive, registry_path) as key:
                    for index in range(winreg.QueryInfoKey(key)[0]):
                        version = winreg.EnumKey(key, index)
                        with winreg.OpenKey(key, version) as version_key:
                            for value_name in ("InstallPath", ""):
                                try:
                                    value, _ = winreg.QueryValueEx(version_key, value_name)
                                except OSError:
                                    continue
                                if value:
                                    roots.append(Path(value))
                                    break
            except OSError:
                continue
    return roots


def find_rscript(search_path=None, install_roots=None, registry_roots=None):
    """Find Rscript through PATH, common installation locations, or Windows registry."""
    found = shutil.which("Rscript", path=search_path)
    if found:
        return str(Path(found).resolve())

    roots = list(install_roots) if install_roots is not None else _default_r_roots()
    if registry_roots is None:
        roots.extend(_registry_r_roots())
    else:
        roots.extend(Path(root) for root in registry_roots)

    candidates = []
    for root in roots:
        root = Path(root)
        candidates.extend((
            root / "bin" / "Rscript.exe",
            root / "bin" / "Rscript",
            root / "bin" / "x64" / "Rscript.exe",
            root / "bin" / "i386" / "Rscript.exe",
        ))
        try:
            candidates.extend(root.glob("R-*/bin/Rscript.exe"))
            candidates.extend(root.glob("R-*/bin/Rscript"))
            candidates.extend(root.glob("R-*/bin/x64/Rscript.exe"))
            candidates.extend(root.glob("R-*/bin/i386/Rscript.exe"))
        except OSError:
            continue
    existing = []
    for candidate in candidates:
        if candidate.is_file():
            try:
                existing.append((candidate.stat().st_mtime, candidate))
            except OSError:
                continue
    return str(max(existing, key=lambda item: item[0])[1]) if existing else ""


def find_cutadapt(resource_root, preferred_path=""):
    """Find a user-selected or packaged Cutadapt executable."""
    candidates = []
    if preferred_path:
        candidates.append(Path(preferred_path))
    candidates.extend((
        Path(resource_root) / "tools" / "cutadapt.exe",
        Path(resource_root) / "cutadapt.exe",
    ))
    found = shutil.which("cutadapt")
    if found:
        candidates.append(Path(found))
    candidates.append(Path.home() / "Downloads" / "cutadapt.exe")
    for candidate in candidates:
        if candidate.is_file():
            return str(candidate.resolve())
    return ""


def verify_executable(executable, timeout=20):
    """Return version output, or raise a useful error if the command cannot run."""
    if not executable or not Path(executable).is_file():
        raise FileNotFoundError(executable)
    completed = subprocess.run(
        [str(executable), "--version"],
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=timeout,
        creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0),
    )
    output = (completed.stdout + "\n" + completed.stderr).strip()
    if completed.returncode != 0 or not output:
        raise RuntimeError(output or f"Process exited with code {completed.returncode}.")
    return output


def install_r_runtime():
    """Download and silently install the current per-user Windows runtime."""
    page_url = "https://cloud.r-project.org/bin/windows/base/"
    with urllib.request.urlopen(page_url, timeout=30) as response:
        page_url = response.geturl()
        page = response.read().decode("utf-8", errors="replace")
    installer_names = re.findall(
        r"""href=["']([^"']*R-\d+\.\d+\.\d+-win\.exe)["']""", page, re.IGNORECASE
    )
    if not installer_names:
        raise RuntimeError("Could not find the Windows runtime installer on the official download page.")

    installer_url = urllib.parse.urljoin(page_url, installer_names[0])
    parsed_url = urllib.parse.urlparse(installer_url)
    if parsed_url.scheme != "https" or parsed_url.hostname not in {
        "cloud.r-project.org",
        "cran.r-project.org",
    }:
        raise RuntimeError("The runtime installer download URL was not trusted.")

    version = re.search(r"R-(\d+\.\d+\.\d+)-win\.exe$", parsed_url.path, re.IGNORECASE)
    if version is None:
        raise RuntimeError("Could not identify the runtime installer version.")
    install_root = Path(os.environ.get("LOCALAPPDATA", Path.home() / "AppData/Local"))
    install_dir = install_root / "Programs" / "R" / f"R-{version.group(1)}"
    with tempfile.TemporaryDirectory(prefix="biome-setup-") as temp_dir:
        installer = Path(temp_dir) / "runtime-setup.exe"
        with urllib.request.urlopen(installer_url, timeout=60) as response, installer.open("wb") as output:
            shutil.copyfileobj(response, output)
        subprocess.run(
            [
                str(installer),
                "/VERYSILENT",
                "/SUPPRESSMSGBOXES",
                "/NORESTART",
                "/SP-",
                "/CURRENTUSER",
                f"/DIR={install_dir}",
            ],
            check=True,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=900,
            creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0),
        )
    return install_dir


class SetupWizard:
    def __init__(
        self,
        parent,
        resource_root,
        rscript_path,
        cutadapt_path,
        on_success,
    ):
        self.parent = parent
        self.resource_root = Path(resource_root)
        self.text = WIZARD_TEXT
        self.rscript_path = rscript_path
        self.on_success = on_success
        self.events = queue.Queue()
        self.cutadapt_path = cutadapt_path or find_cutadapt(self.resource_root)
        self.window = None
        self.prompt = None
        self.installing = False
        self.progress = 0
        self._show_loading()
        self.window.grab_set()
        self._start_worker(install=False)

    def _show_loading(self):
        self.window = tk.Toplevel(self.parent)
        if os.name == "nt":
            self.window.iconbitmap(str(self.resource_root / "assets" / "biome.ico"))
        self.window.title(self.text["title"])
        self.window.geometry("340x125")
        self.window.resizable(False, False)
        if self.parent.state() != "withdrawn":
            self.window.transient(self.parent)
        self.window.protocol("WM_DELETE_WINDOW", lambda: None)
        self.progress_label = ttk.Label(
            self.window,
            text=self.text["progress"].format(percent=0),
            padding=(18, 14, 18, 6),
        )
        self.progress_label.pack(fill="x")
        self.progress_bar = ttk.Progressbar(self.window, mode="determinate", maximum=100, value=0)
        self.progress_bar.pack(fill="x", padx=18, pady=(0, 16))
        self.window.after(100, self.poll)

    def _start_worker(self, install):
        self.installing = install
        self._set_progress(0, "Installing required components" if install else "Checking required components")
        threading.Thread(target=self._run_setup, args=(install,), daemon=True).start()

    def _set_progress(self, percent, message=None):
        """Update the wizard from the Tk event loop using queued worker events."""
        self.events.put(("progress", (percent, message)))

    def _run_setup(self, install):
        details = []
        try:
            self._set_progress(10, "Checking required components")
            result = self._check_setup(self.rscript_path)
            details.extend(result["details"])
            self._set_progress(40, "Checking installed packages")
            if install:
                if "runtime" in result["missing"]:
                    self._set_progress(50, "Installing R")
                    install_r_runtime()
                    result["rscript"] = find_rscript()
                if result["missing"]:
                    if not result["rscript"]:
                        raise RuntimeError("The runtime could not be found after installation.")
                    self._set_progress(75, "Installing required packages")
                    self._run_installer(result["rscript"], details)
                self._set_progress(95, "Verifying installation")
                result = self._check_setup(result["rscript"])
                details.extend(result["details"])
                if result["missing"]:
                    raise RuntimeError("Required components are still missing after installation.")
            self.events.put(("progress", (100, "Setup complete")))
            self.events.put(("complete", result))
        except (
            OSError,
            RuntimeError,
            subprocess.SubprocessError,
            urllib.error.URLError,
            ValueError,
        ) as exc:
            if isinstance(exc, subprocess.CalledProcessError):
                details.extend(
                    output.strip()
                    for output in (exc.stdout, exc.stderr)
                    if output and output.strip()
                )
            details.append(str(exc))
            self.events.put(("failure", "\n".join(details)))

    def _check_setup(self, preferred_rscript=""):
        required_scripts = [
            self.resource_root / "pipeline" / "install_packages.R",
            self.resource_root / "run_pipeline.R",
        ]
        required_scripts.extend(
            self.resource_root / "pipeline" / f"{number:02d}_{name}.R"
            for number, name in enumerate((
                "parameters", "helpers", "quality_profiles", "primer_removal",
                "filter_trim", "learn_errors", "denoise_merge",
                "chimera_length_filter", "read_tracking", "taxonomy",
                "phyloseq_build", "phylogenetic_tree", "rarefaction",
                "alpha_diversity", "beta_diversity", "relative_abundance", "ancombc",
            ))
        )
        missing_files = [path for path in required_scripts if not path.is_file()]
        if missing_files:
            raise FileNotFoundError("Bundled application files are missing: " + ", ".join(map(str, missing_files)))

        candidates = []
        if preferred_rscript and Path(preferred_rscript).is_file():
            candidates.append(preferred_rscript)
        discovered = find_rscript()
        if discovered and discovered not in candidates:
            candidates.append(discovered)
        cutadapt = find_cutadapt(self.resource_root, self.cutadapt_path)
        missing = []
        details = []
        rscript = ""
        for candidate in candidates:
            try:
                verify_executable(candidate)
                rscript = candidate
                break
            except (OSError, RuntimeError, subprocess.SubprocessError) as exc:
                details.append(str(exc))
        if not rscript:
            missing.append("runtime")
        else:
            completed = subprocess.run(
                [rscript, "--vanilla", str(self.resource_root / "pipeline" / "install_packages.R"), "--check"],
                cwd=str(self.resource_root),
                capture_output=True,
                text=True,
                encoding="utf-8",
                errors="replace",
                timeout=300,
                creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0),
            )
            if completed.returncode != 0:
                missing.append("packages")
                details.append((completed.stdout + "\n" + completed.stderr).strip())

        if not cutadapt:
            missing.append("bundled_tool")
        else:
            try:
                verify_executable(cutadapt)
            except (OSError, RuntimeError, subprocess.SubprocessError) as exc:
                missing.append("bundled_tool")
                details.append(str(exc))
        return {
            "rscript": rscript,
            "cutadapt": cutadapt,
            "missing": missing,
            "details": details,
        }

    def _run_installer(self, rscript, details):
        completed = subprocess.run(
            [rscript, "--vanilla", str(self.resource_root / "pipeline" / "install_packages.R")],
            cwd=str(self.resource_root),
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=3600,
            creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0),
        )
        output = (completed.stdout + "\n" + completed.stderr).strip()
        if output:
            details.append(output)
        if completed.returncode != 0:
            raise RuntimeError(f"Component installation failed with exit code {completed.returncode}.")

    def poll(self):
        try:
            while True:
                event, value = self.events.get_nowait()
                if event == "progress":
                    percent, message = value
                    self.progress = percent
                    if self.window is not None and self.window.winfo_exists():
                        self.progress_bar.configure(value=percent)
                        label = message or ("Installing required components" if self.installing else "Checking required components")
                        self.progress_label.configure(text=f"{label}: {percent}%")
                elif event == "failure":
                    self.window.grab_release()
                    self.window.destroy()
                    if self.parent.state() == "withdrawn":
                        self.parent.deiconify()
                    messagebox.showerror(
                        self.text["failure_title"],
                        self.text["failure_message"] + value,
                        parent=self.parent,
                    )
                    self.parent.destroy()
                elif event == "complete":
                    if value["missing"]:
                        if self.installing:
                            raise RuntimeError("Setup completed with missing components.")
                        self.window.grab_release()
                        self.window.destroy()
                        self._show_install_prompt()
                    else:
                        self.window.grab_release()
                        self.window.destroy()
                        self.on_success(value["rscript"], value["cutadapt"])
        except queue.Empty:
            pass
        if self.window is not None and self.window.winfo_exists():
            self.window.after(100, self.poll)

    def _show_install_prompt(self):
        self.prompt = tk.Toplevel(self.parent)
        if os.name == "nt":
            self.prompt.iconbitmap(str(self.resource_root / "assets" / "biome.ico"))
        self.prompt.title(self.text["prompt_title"])
        self.prompt.geometry("390x190")
        self.prompt.resizable(False, False)
        if self.parent.state() != "withdrawn":
            self.prompt.transient(self.parent)
        self.prompt.protocol("WM_DELETE_WINDOW", self.cancel)
        ttk.Label(
            self.prompt,
            text=self.text["prompt_message"],
            wraplength=340,
            justify="left",
            padding=(22, 20),
        ).pack(fill="both", expand=True)
        buttons = ttk.Frame(self.prompt, padding=(16, 0, 16, 14))
        buttons.pack(fill="x")
        ttk.Button(buttons, text=self.text["cancel"], command=self.cancel).pack(side="right")
        ttk.Button(buttons, text=self.text["install"], command=self.install).pack(
            side="right", padx=(0, 8)
        )
        self.prompt.deiconify()
        self.prompt.update_idletasks()
        self.prompt.grab_set()
        self.prompt.lift()
        self.prompt.focus_force()

    def install(self):
        self.prompt.grab_release()
        self.prompt.destroy()
        self._show_loading()
        self.window.grab_set()
        self._start_worker(install=True)

    def cancel(self):
        if self.prompt is not None and self.prompt.winfo_exists():
            self.prompt.grab_release()
            self.prompt.destroy()
        self.parent.destroy()
