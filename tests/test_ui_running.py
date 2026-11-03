"""Check real Tk control locking and restoration without running an analysis."""
import unittest
from unittest.mock import patch

# Share the existing loader so startup exception hooks target the same module.
from test_launcher import launcher


class RunningControlsTests(unittest.TestCase):
    def setUp(self):
        try:
            self.root = launcher.tk.Tk()
        except launcher.tk.TclError as error:
            self.skipTest(f"Tk display unavailable: {error}")
        self.root.withdraw()
        self.app = launcher.PipelineApp(self.root)

    def tearDown(self):
        if hasattr(self, "root"):
            self.root.destroy()

    def interactive_widgets(self):
        result = []
        pending = [self.app.main_frame]
        while pending:
            widget = pending.pop()
            pending.extend(widget.winfo_children())
            if isinstance(widget, (launcher.ttk.Button, launcher.ttk.Checkbutton,
                                   launcher.ttk.Radiobutton, launcher.ttk.Entry,
                                   launcher.ttk.Combobox)):
                result.append(widget)
        return result

    def test_all_analysis_controls_lock_and_restore_current_states(self):
        self.app.vars["truncation_mode"].set("manual")
        self.app.update_truncation_fields()
        self.app.asv_length_filter_enabled.set(True)
        self.app.update_asv_length_fields()
        widgets = self.interactive_widgets()
        before = {widget: str(widget.cget("state")) for widget in widgets}
        radio_labels = {str(widget.cget("text")) for widget in widgets
                        if isinstance(widget, launcher.ttk.Radiobutton)}
        self.assertIn("I know the primer sequences", radio_labels)
        self.assertIn("Filter ASVs by length", radio_labels)
        self.assertEqual(before[self.app.trunc_qmin_entry], "disabled")
        self.assertEqual(before[self.app.asv_min_len_entry], "normal")
        self.app.set_running_controls(True)
        for widget in widgets:
            self.assertEqual(widget.instate(["disabled"]),
                             widget is not self.app.stop_button, str(widget))
        self.app.set_running_controls(False)
        for widget, state in before.items():
            self.assertEqual(str(widget.cget("state")), state, str(widget))

    def test_zero_abundance_threshold_shows_warning_without_changing_values(self):
        self.assertFalse(self.app.abundance_warning_label.winfo_manager())
        self.app.vars["min_rank2_abundance_percent"].set("0")
        self.assertEqual(self.app.abundance_warning_label.winfo_manager(), "grid")
        self.assertIn("0% keeps all assigned taxa", self.app.abundance_warning_label.cget("text"))
        self.assertEqual(self.app.vars["min_rank2_abundance_percent"].get(), "0")
        self.app.vars["min_rank2_abundance_percent"].set("2")
        self.assertFalse(self.app.abundance_warning_label.winfo_manager())
        self.app.vars["min_rank1_abundance_percent"].set("0.0")
        self.assertEqual(self.app.abundance_warning_label.winfo_manager(), "grid")

    def test_runtime_plot_warning_is_displayed(self):
        self.app.events.put(("plot_warning", "Increase the collapse percentage."))
        with patch.object(launcher.messagebox, "showwarning") as show_warning:
            self.app.poll()
        show_warning.assert_called_once_with(
            "Relative abundance plots", "Increase the collapse percentage.", parent=self.root
        )

    def test_completion_cancellation_and_failure_restore_controls(self):
        self.app.last_run = "example-output"
        widgets = self.interactive_widgets()
        for cancelled, exit_code in ((False, 0), (True, 1), (False, 1)):
            with self.subTest(cancelled=cancelled, exit_code=exit_code):
                before = {widget: str(widget.cget("state")) for widget in widgets}
                self.app.busy = True
                self.app.set_running_controls(True)
                self.app.events.put(("done", {"cancelled": cancelled,
                                              "exit_code": exit_code}))
                with patch.object(launcher.messagebox, "showinfo"), patch.object(
                    launcher.messagebox, "showerror"
                ):
                    self.app.poll()
                self.assertFalse(self.app.busy)
                for widget, state in before.items():
                    self.assertEqual(str(widget.cget("state")), state, str(widget))


if __name__ == "__main__":
    unittest.main()
