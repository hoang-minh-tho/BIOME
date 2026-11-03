"""Tests for dependency discovery and the automatic setup flow."""

from pathlib import Path
import tempfile
import unittest
from unittest.mock import Mock, patch

import app_setup


class DependencyDiscoveryTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def make_executable(self, path):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"test")
        return path

    def create_bundled_files(self):
        (self.root / "run_pipeline.R").write_text("", encoding="utf-8")
        pipeline = self.root / "pipeline"
        pipeline.mkdir(exist_ok=True)
        (pipeline / "install_packages.R").write_text("", encoding="utf-8")
        for number, name in enumerate((
            "parameters", "helpers", "quality_profiles", "primer_removal",
            "filter_trim", "learn_errors", "denoise_merge",
            "chimera_length_filter", "read_tracking", "taxonomy",
            "phyloseq_build", "phylogenetic_tree", "rarefaction",
            "alpha_diversity", "beta_diversity", "relative_abundance", "ancombc",
        )):
            (pipeline / f"{number:02d}_{name}.R").write_text("", encoding="utf-8")

    def test_find_rscript_in_standard_install_root(self):
        executable = self.make_executable(self.root / "R-4.6.1" / "bin" / "Rscript.exe")
        found = app_setup.find_rscript(
            search_path="",
            install_roots=[self.root],
            registry_roots=[],
        )
        self.assertTrue(Path(found).samefile(executable))

    def test_find_rscript_in_registry_root(self):
        executable = self.make_executable(self.root / "bin" / "Rscript.exe")
        found = app_setup.find_rscript(
            search_path="",
            install_roots=[],
            registry_roots=[self.root],
        )
        self.assertEqual(Path(found), executable)

    def test_find_cutadapt_prefers_packaged_executable(self):
        bundled = self.make_executable(self.root / "tools" / "cutadapt.exe")
        with patch.object(app_setup.shutil, "which", return_value=""):
            found = app_setup.find_cutadapt(self.root)
        self.assertTrue(Path(found).samefile(bundled))

    def test_find_cutadapt_uses_saved_executable(self):
        preferred = self.make_executable(self.root / "custom" / "cutadapt.exe")
        with patch.object(app_setup.shutil, "which", return_value=""):
            found = app_setup.find_cutadapt(self.root, preferred)
        self.assertTrue(Path(found).samefile(preferred))

    def test_verify_executable_reports_version(self):
        executable = self.make_executable(self.root / "tool.exe")
        with patch.object(
            app_setup.subprocess,
            "run",
            return_value=app_setup.subprocess.CompletedProcess(
                ["tool", "--version"], 0, stdout="Tool 1.2\n", stderr=""
            ),
        ):
            self.assertEqual(app_setup.verify_executable(str(executable)), "Tool 1.2")

    def test_verify_executable_rejects_nonzero_exit(self):
        executable = self.make_executable(self.root / "tool.exe")
        with patch.object(
            app_setup.subprocess,
            "run",
            return_value=app_setup.subprocess.CompletedProcess(
                ["tool", "--version"], 1, stdout="", stderr="not runnable"
            ),
        ):
            with self.assertRaisesRegex(RuntimeError, "not runnable"):
                app_setup.verify_executable(str(executable))

    def test_startup_check_is_silent_and_checks_bundled_components(self):
        self.create_bundled_files()
        rscript = self.make_executable(self.root / "Rscript.exe")
        cutadapt = self.make_executable(self.root / "tools" / "cutadapt.exe")
        wizard = app_setup.SetupWizard.__new__(app_setup.SetupWizard)
        wizard.resource_root = self.root
        wizard.cutadapt_path = str(cutadapt)

        with patch.object(app_setup, "find_rscript", return_value=str(rscript)), \
                patch.object(app_setup, "verify_executable", return_value="1.0"), \
                patch.object(app_setup.subprocess, "run", return_value=app_setup.subprocess.CompletedProcess(
                    ["Rscript"], 0, stdout="ready", stderr=""
                )) as run:
            result = wizard._check_setup()

        self.assertEqual(result["missing"], [])
        self.assertEqual(result["rscript"], str(rscript))
        self.assertTrue(Path(result["cutadapt"]).samefile(cutadapt))
        self.assertIn("--check", run.call_args.args[0])

    def test_startup_check_reports_missing_runtime_without_showing_details(self):
        self.create_bundled_files()
        wizard = app_setup.SetupWizard.__new__(app_setup.SetupWizard)
        wizard.resource_root = self.root
        wizard.cutadapt_path = ""

        with patch.object(app_setup, "find_rscript", return_value=""):
            result = wizard._check_setup()

        self.assertIn("runtime", result["missing"])
        self.assertEqual(result["details"], [])

    def test_missing_runtime_shows_install_prompt(self):
        wizard = app_setup.SetupWizard.__new__(app_setup.SetupWizard)
        wizard.events = app_setup.queue.Queue()
        wizard.events.put(("complete", {"missing": ["runtime"], "rscript": "", "cutadapt": ""}))
        wizard.installing = False
        wizard.window = Mock()
        wizard.window.winfo_exists.return_value = False
        wizard.parent = Mock()
        wizard.text = app_setup.WIZARD_TEXT
        wizard._show_install_prompt = Mock()

        wizard.poll()

        wizard.window.grab_release.assert_called_once_with()
        wizard.window.destroy.assert_called_once_with()
        wizard._show_install_prompt.assert_called_once_with()

    def test_install_prompt_is_visible_when_parent_is_withdrawn(self):
        wizard = app_setup.SetupWizard.__new__(app_setup.SetupWizard)
        wizard.parent = Mock()
        wizard.parent.state.return_value = "withdrawn"
        wizard.text = app_setup.WIZARD_TEXT
        wizard.resource_root = self.root

        with patch.object(app_setup.tk, "Toplevel") as toplevel, \
                patch.object(app_setup.ttk, "Label"), \
                patch.object(app_setup.ttk, "Frame"), \
                patch.object(app_setup.ttk, "Button"):
            wizard._show_install_prompt()

        prompt = toplevel.return_value
        prompt.transient.assert_not_called()
        prompt.deiconify.assert_called_once_with()
        prompt.update_idletasks.assert_called_once_with()
        prompt.grab_set.assert_called_once_with()
        prompt.lift.assert_called_once_with()
        prompt.focus_force.assert_called_once_with()

    def test_loading_window_shows_percentage_progress(self):
        wizard = app_setup.SetupWizard.__new__(app_setup.SetupWizard)
        wizard.parent = Mock()
        wizard.parent.state.return_value = "withdrawn"
        wizard.text = app_setup.WIZARD_TEXT
        wizard.resource_root = self.root
        wizard.events = app_setup.queue.Queue()
        wizard.installing = True

        with patch.object(app_setup.tk, "Toplevel"), \
                patch.object(app_setup.ttk, "Label") as label, \
                patch.object(app_setup.ttk, "Progressbar") as progressbar:
            wizard._show_loading()
            wizard.window.winfo_exists.return_value = True
            wizard.events.put(("progress", (62, "Installing R")))
            wizard.poll()

        progressbar.assert_called_once_with(
            wizard.window, mode="determinate", maximum=100, value=0
        )
        wizard.progress_bar.configure.assert_called_once_with(value=62)
        label.return_value.configure.assert_called_once_with(text="Installing R: 62%")

    def test_first_run_copy_is_minimal(self):
        self.assertEqual(app_setup.WIZARD_TEXT["loading"], "Loading...")
        self.assertEqual(app_setup.WIZARD_TEXT["title"], "Preparing BIOME")
        self.assertEqual(app_setup.WIZARD_TEXT["prompt_title"], "Additional components required")
        self.assertEqual(
            app_setup.WIZARD_TEXT["prompt_message"],
            "BIOME needs to install required components to run.\n\nAllow installation now?",
        )
        self.assertEqual(app_setup.WIZARD_TEXT["install"], "Install")
        self.assertEqual(app_setup.WIZARD_TEXT["cancel"], "Cancel")
        self.assertFalse(any(
            word in " ".join(app_setup.WIZARD_TEXT.values()).casefold()
            for word in ("rscript", "cran", "bioconductor", "browse", "download", "console")
        ))


if __name__ == "__main__":
    unittest.main(verbosity=2)
