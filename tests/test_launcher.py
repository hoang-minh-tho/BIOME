"""Input validation and subprocess orchestration; no sequencing data needed."""
import importlib.util
import inspect
import json
from pathlib import Path
import queue
import tempfile
import threading
from types import SimpleNamespace
import unittest
import zipfile
from unittest.mock import Mock, patch

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("launcher", ROOT / "run_pipeline.py")
launcher = importlib.util.module_from_spec(spec)
spec.loader.exec_module(launcher)


class Value:
    def __init__(self, value):
        self.value = value

    def get(self):
        return self.value

    def set(self, value):
        self.value = value


class Widget:
    def __init__(self):
        self.visible = True
        self.options = {}

    def grid(self):
        self.visible = True

    def grid_remove(self):
        self.visible = False

    def configure(self, **options):
        self.options.update(options)


class LauncherTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.folder = Path(self.temp.name) / "Dữ liệu có dấu cách"
        self.folder.mkdir()

    def tearDown(self):
        self.temp.cleanup()

    def test_biome_branding_and_logo_are_available(self):
        self.assertEqual(
            launcher.MESSAGES["window_title"],
            "BIOME — Bioinformatic Interface for Optimized Multi-marker Ecology",
        )
        self.assertTrue((ROOT / "assets" / "biome.ico").is_file())
        ui_source = inspect.getsource(launcher.PipelineApp.build_ui)
        self.assertNotIn("PhotoImage", ui_source)
        self.assertNotIn("brand_image", ui_source)
        self.assertIn("self.notebook.grid(row=0, column=0", ui_source)
        self.assertIn('self.root.title(self.text("window_title"))', ui_source)

    def make_files(self, names):
        paths = []
        for name in names:
            path = self.folder / name
            path.write_bytes(b"placeholder")
            paths.append(str(path))
        return paths

    def make_configuration_app(self):
        samples = ("A1", "A2", "A3", "B1", "B2")
        fastq_files = self.make_files([
            filename
            for sample in samples
            for filename in (f"{sample}_R1.fastq.gz", f"{sample}_R2.fastq.gz")
        ])
        metadata = self.folder / "metadata.csv"
        metadata.write_text(
            "SampleID,Group,Control\n"
            "A1,A,Negative\nA2,A,Sample\nA3,A,Mock\n"
            "B1,B,Sample\nB2,B,Sample\n",
            encoding="utf-8",
        )
        rscript = self.folder / "Rscript.exe"
        cutadapt = self.folder / "primer-tool.exe"
        reference = self.folder / "training.fa.gz"
        for path in (rscript, cutadapt, reference):
            path.write_bytes(b"placeholder")
        app = launcher.PipelineApp.__new__(launcher.PipelineApp)
        app.vars = {
            "rscript": str(rscript),
            "cutadapt": str(cutadapt),
            "fastq_dir": str(self.folder),
            "metadata_file": str(metadata),
            "reference_db_path": str(reference),
            "marker": "16S",
            "output_dir": str(self.folder / "output"),
            "group_column": "Group",
            "reference_group": "A",
            "primer_fwd": "ACGT",
            "primer_rev": "TGCA",
            "primer_mode": "sequences",
            "trim_left_fwd": "19",
            "trim_left_rev": "20",
            "truncation_mode": "automatic",
            "trunc_qmin": "30",
            "trunc_len_fwd": "240",
            "trunc_len_rev": "200",
            "max_ee_fwd": "2",
            "max_ee_rev": "2",
            "pooling": "pseudo",
            "max_n": "0",
            "trunc_q": "2",
            "min_len": "50",
            "asv_min_len": "230",
            "asv_max_len": "280",
            "chimera_method": "Consensus",
            "mergepairs_strategy": "Merge paired reads",
            "mergepairs_min_overlap": "12",
            "mergepairs_max_mismatch": "0",
            "rarefaction_mode": "Automatic",
            "custom_rarefaction_depth": "10000",
            "max_sample_loss_percent": "10",
            "plateau_gain_percent": "2",
            "plateau_samples_percent": "90",
            "read_step": "1000",
            "permutations": "999",
            "ordination": "PCoA",
            "differential_method": "ANCOM-BC",
            "min_rank1_abundance_percent": "1.0",
            "min_rank2_abundance_percent": "2.0",
            "figure_format": "PNG",
            "mock_control_mode": "Metadata column",
            "mock_control_column": "Control",
            "mock_control_level": "Mock",
            "negative_control_mode": "Metadata column",
            "negative_control_column": "Control",
            "negative_control_level": "Negative",
        }
        app.vars = {key: Value(value) for key, value in app.vars.items()}
        app.step_vars = {f"step_{number}": Value(True) for number in (1, 2, 3)}
        app.step_checkboxes = {key: Widget() for key in app.step_vars}
        app.remove_phix = Value(True)
        app.asv_length_filter_enabled = Value(False)
        app.run_permanova = Value(True)
        app.run_betadisper = Value(True)
        app.beta_vars = {
            "bray": Value(True), "wunifrac": Value(False),
            "unifrac": Value(False), "jaccard": Value(False),
        }
        app.ordination_var = Value("PCoA")
        app.differential_var = Value("ANCOM-BC")
        app.rarefaction_var = app.vars["rarefaction_mode"]
        app.remove_phix = Value(True)
        app.run_permanova = Value(True)
        app.run_betadisper = Value(True)
        return app, fastq_files

    def test_thread_count_is_automatic_without_a_user_setting(self):
        app, _fastq_files = self.make_configuration_app()
        for cpu_count, expected in ((None, 1), (1, 1), (2, 1), (8, 7)):
            with self.subTest(cpu_count=cpu_count):
                with patch.object(launcher.os, "cpu_count", return_value=cpu_count):
                    self.assertEqual(app.configuration()["threads"], expected)

    def test_legacy_thread_setting_does_not_override_automatic_count(self):
        app, _fastq_files = self.make_configuration_app()
        app.vars["threads"] = Value("invalid legacy setting")
        with patch.object(launcher.os, "cpu_count", return_value=8):
            self.assertEqual(app.configuration()["threads"], 7)

    def test_pipeline_step_dependencies_are_enforced_and_restored(self):
        app = launcher.PipelineApp.__new__(launcher.PipelineApp)
        app.step_vars = {
            "step_1": Value(False),
            "step_2": Value(True),
            "step_3": Value(False),
        }
        app.step_checkboxes = {key: Widget() for key in app.step_vars}

        app.update_pipeline_steps()
        self.assertTrue(app.step_vars["step_1"].get())
        self.assertEqual(app.step_checkboxes["step_1"].options["state"], "disabled")
        self.assertEqual(app.step_checkboxes["step_2"].options["state"], "normal")

        app.step_vars["step_2"].set(False)
        app.update_pipeline_steps()
        self.assertEqual(app.step_checkboxes["step_1"].options["state"], "normal")

        app.step_vars["step_3"].set(True)
        app.update_pipeline_steps()
        self.assertTrue(app.step_vars["step_1"].get())
        self.assertTrue(app.step_vars["step_2"].get())
        self.assertEqual(app.step_checkboxes["step_1"].options["state"], "disabled")
        self.assertEqual(app.step_checkboxes["step_2"].options["state"], "disabled")

        app.step_vars["step_3"].set(False)
        app.update_pipeline_steps()
        self.assertEqual(app.step_checkboxes["step_1"].options["state"], "disabled")
        self.assertEqual(app.step_checkboxes["step_2"].options["state"], "normal")

    def test_configuration_repairs_invalid_pipeline_selection(self):
        app, _fastq_files = self.make_configuration_app()
        app.step_vars["step_1"].set(False)
        app.step_vars["step_2"].set(False)
        app.step_vars["step_3"].set(True)

        config = app.configuration()

        self.assertEqual(config["steps"], {
            "step_1": True,
            "step_2": True,
            "step_3": True,
        })

    def test_configuration_has_no_cache_options(self):
        app, _fastq_files = self.make_configuration_app()

        config = app.configuration()

        for key in ("enable_cache", "force_recompute", "cache_dir"):
            self.assertNotIn(key, config)
        ui_source = inspect.getsource(launcher.PipelineApp.build_ui)
        self.assertNotIn('text="Project cache"', ui_source)

    def test_marker_profile_selection_is_written_to_r_configuration(self):
        app, _fastq_files = self.make_configuration_app()
        app.vars["marker"].set("18S")

        config = app.configuration()

        self.assertEqual(config["marker"], "18S")
        self.assertEqual(config["marker_profile"]["target_domain"], "Eukaryota")
        self.assertEqual(config["marker_profile"]["target_domain_rank"], "Domain")
        self.assertEqual(config["marker_profile"]["relative_abundance_ranks"], ["Division", "Genus"])
        self.assertTrue(Path(config["marker_profiles_path"]).is_file())

    def test_marker_profiles_keep_taxonomy_settings_out_of_pipeline_code(self):
        sixteen_s = launcher.get_marker_profile("16S")
        eighteen_s = launcher.get_marker_profile("18S")
        self.assertEqual(sixteen_s["target_domain"], "Bacteria")
        self.assertEqual(eighteen_s["target_domain"], "Eukaryota")
        self.assertEqual(eighteen_s["taxonomy_ranks"][0], "Domain")
        self.assertEqual(eighteen_s["relative_abundance_ranks"][0], "Division")

        phyloseq_script = (ROOT / "pipeline" / "10_phyloseq_build.R").read_text(encoding="utf-8")
        abundance_script = (ROOT / "pipeline" / "15_relative_abundance.R").read_text(encoding="utf-8")
        differential_script = (ROOT / "pipeline" / "16_ancombc.R").read_text(encoding="utf-8")
        self.assertIn("tax_df[[target_domain_rank]]", phyloseq_script)
        self.assertIn("taxrank = taxrank", abundance_script)
        self.assertIn("tax_level = relative_abundance_ranks[2]", differential_script)
        self.assertNotIn("Bacteria", phyloseq_script)

    def test_asv_length_filter_validates_and_serializes_enabled_bounds(self):
        app, _fastq_files = self.make_configuration_app()
        app.asv_length_filter_enabled.set(True)

        config = app.configuration()

        self.assertTrue(config["asv_length_filter_enabled"])
        self.assertEqual(config["asv_min_len"], 230)
        self.assertEqual(config["asv_max_len"], 280)

        app.vars["asv_max_len"].set("230")
        with self.assertRaisesRegex(ValueError, "Minimum ASV length must be less"):
            app.configuration()

    def test_asv_length_pipeline_bypasses_filter_when_disabled(self):
        length_script = (ROOT / "pipeline" / "07_chimera_length_filter.R").read_text(encoding="utf-8")
        self.assertIn("if (!isTRUE(asv_length_filter_enabled))", length_script)
        self.assertIn("return(seqtab_nochim)", length_script)
        self.assertIn("asv_lengths >= asv_min_len", length_script)
        self.assertIn("asv_lengths <= asv_max_len", length_script)

    def test_asv_length_filter_migrates_legacy_settings(self):
        self.assertFalse(launcher.saved_asv_length_filter_enabled({}))
        self.assertTrue(launcher.saved_asv_length_filter_enabled({
            "asv_min_len": "230",
            "asv_max_len": "280",
        }))
        self.assertFalse(launcher.saved_asv_length_filter_enabled({
            "asv_length_filter_enabled": False,
            "asv_min_len": "230",
            "asv_max_len": "280",
        }))

    def test_legacy_reference_and_abundance_settings_are_accepted(self):
        app, _fastq_files = self.make_configuration_app()
        app.vars["silva_ref_path"] = app.vars.pop("reference_db_path")
        app.vars["min_phylum_abundance_percent"] = app.vars.pop(
            "min_rank1_abundance_percent"
        )
        app.vars["min_genus_abundance_percent"] = app.vars.pop(
            "min_rank2_abundance_percent"
        )
        app.vars.pop("marker")

        config = app.configuration()

        self.assertEqual(config["marker"], "16S")
        self.assertEqual(config["reference_db_path"], str((self.folder / "training.fa.gz").resolve()).replace("\\", "/"))
        self.assertEqual(config["min_rank1_abundance_percent"], 1.0)
        self.assertEqual(config["min_rank2_abundance_percent"], 2.0)

    def test_rarefaction_fields_follow_selected_mode(self):
        app = launcher.PipelineApp.__new__(launcher.PipelineApp)
        app.rarefaction_var = Value("Automatic")
        app.custom_depth_row = Widget()
        app.minimum_depth_description = Widget()
        app.rarefaction_advanced_toggle = Widget()
        app.rarefaction_advanced = Widget()
        app.rarefaction_advanced_entries = [Widget() for _ in range(4)]
        app.rarefaction_advanced_open = False

        app.update_rarefaction_fields()
        self.assertTrue(app.rarefaction_advanced_toggle.visible)
        self.assertFalse(app.rarefaction_advanced.visible)
        self.assertTrue(all(entry.options["state"] == "normal" for entry in app.rarefaction_advanced_entries))

        app.rarefaction_var.set("Minimum depth")
        app.update_rarefaction_fields()
        self.assertFalse(app.rarefaction_advanced_toggle.visible)
        self.assertFalse(app.rarefaction_advanced.visible)
        self.assertTrue(app.minimum_depth_description.visible)
        self.assertFalse(app.custom_depth_row.visible)
        self.assertTrue(all(entry.options["state"] == "disabled" for entry in app.rarefaction_advanced_entries))

        app.rarefaction_var.set("Custom depth")
        app.update_rarefaction_fields()
        self.assertFalse(app.rarefaction_advanced_toggle.visible)
        self.assertFalse(app.minimum_depth_description.visible)
        self.assertTrue(app.custom_depth_row.visible)

        app.rarefaction_var.set("Automatic")
        app.rarefaction_advanced_open = True
        app.update_rarefaction_fields()
        self.assertTrue(app.rarefaction_advanced.visible)
        self.assertFalse(app.custom_depth_row.visible)
        self.assertFalse(app.minimum_depth_description.visible)

    def test_pairing_preserves_full_sample_id_and_zero_padding(self):
        files = self.make_files(["001_R1_001.fastq.gz", "001_R2_001.fastq.gz",
                                 "site_A_1.fq", "site_A_2.fq"])
        pairs = launcher.pair_files(files)
        self.assertEqual(set(pairs), {"001", "site_A"})
        self.assertEqual(set(pairs["001"]), {"1", "2"})

    def test_illumina_lane_suffix_matches_metadata_sample_id(self):
        mapping = launcher.match_metadata_sample_ids(
            ["14SB_S18_L001", "14WCB_S104_L001"],
            ["14SB", "14WCB"],
        )
        self.assertEqual(mapping, {
            "14SB_S18_L001": "14SB",
            "14WCB_S104_L001": "14WCB",
        })

    def test_barcode_lane_suffix_matches_metadata_sample_id(self):
        mapping = launcher.match_metadata_sample_ids(
            ["JCAEUKECH_AGTGCT_L001"], ["JCAEUKECH_AGTGCT"]
        )
        self.assertEqual(mapping, {"JCAEUKECH_AGTGCT_L001": "JCAEUKECH_AGTGCT"})

    def test_lane_only_exact_match_takes_precedence(self):
        mapping = launcher.match_metadata_sample_ids(
            ["sample_L001"], ["sample_L001", "sample"]
        )
        self.assertEqual(mapping, {"sample_L001": "sample_L001"})

    def test_multiple_lanes_for_one_metadata_sample_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "Multiple FASTQ sample names"):
            launcher.match_metadata_sample_ids(
                ["sample_L001", "sample_L002"], ["sample"]
            )

    def test_exact_metadata_sample_id_takes_precedence(self):
        mapping = launcher.match_metadata_sample_ids(
            ["14SB_S18_L001"],
            ["14SB_S18_L001", "14SB"],
        )
        self.assertEqual(mapping, {"14SB_S18_L001": "14SB_S18_L001"})

    def test_duplicate_lane_suffix_matches_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "Multiple FASTQ sample names"):
            launcher.match_metadata_sample_ids(
                ["sample_S1_L001", "sample_S2_L001"],
                ["sample"],
            )

    def test_control_prefixes_are_trimmed_deduplicated_and_case_insensitive(self):
        prefixes = launcher.parse_control_prefixes(" Blank, Neg, Blank ,, Water ")
        self.assertEqual(prefixes, ["Blank", "Neg", "Water"])
        self.assertTrue(launcher.is_control_sample("bLaNk_001", prefixes))
        self.assertTrue(launcher.is_control_sample("NEG-2", prefixes))
        self.assertFalse(launcher.is_control_sample("sample_blank", prefixes))
        self.assertFalse(launcher.is_control_sample("sample", []))

    def test_missing_reverse_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "Missing R1/R2"):
            launcher.pair_files(self.make_files(["sample_R1.fastq.gz"]))

    def test_duplicate_sample_across_folders_is_rejected(self):
        files = self.make_files(["sample_R1.fastq.gz", "sample_R2.fastq.gz"])
        other = self.folder / "other"
        other.mkdir()
        duplicate = other / "sample_R1.fastq.gz"
        duplicate.write_bytes(b"placeholder")
        with self.assertRaisesRegex(ValueError, "Duplicate sample"):
            launcher.pair_files(files + [str(duplicate)])

    def test_non_fastq_name_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "FASTQ names"):
            launcher.pair_files(self.make_files(["metadata.csv"]))

    def test_empty_file_is_rejected(self):
        path = self.folder / "empty_R1.fastq"
        path.touch()
        with self.assertRaisesRegex(ValueError, "empty"):
            launcher.pair_files([str(path)])

    def test_metadata_bom_leading_zeros_and_unicode(self):
        path = self.folder / "metadata.csv"
        path.write_text("SampleID,Group\n001,Đối chứng\n002,Điều trị\n", encoding="utf-8-sig")
        columns, rows = launcher.read_metadata(path)
        self.assertEqual(columns, ["SampleID", "Group"])
        self.assertEqual(rows[0]["SampleID"], "001")
        self.assertEqual(rows[0]["Group"], "Đối chứng")

    def test_metadata_must_be_csv_with_sample_id_column(self):
        tsv_path = self.folder / "metadata.tsv"
        tsv_path.write_text("SampleID\tGroup\nsample1\tA\nsample2\tB\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, r"\.csv extension"):
            launcher.read_metadata(tsv_path)

        missing_sample_id = self.folder / "no_sample_id.csv"
        missing_sample_id.write_text("ID,Group\nsample1,A\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "SampleID column"):
            launcher.read_metadata(missing_sample_id)

    def test_metadata_csv_rejects_invalid_csv_syntax(self):
        path = self.folder / "project_data.csv"
        path.write_text('SampleID,Group\nsample1,"unterminated\n', encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "Could not read metadata CSV"):
            launcher.read_metadata(path)

    def test_metadata_csv_can_have_arbitrary_filename_and_sample_id_column_position(self):
        path = self.folder / "coral_samples.csv"
        path.write_text("Group,SampleID,Control\nA,sample1,No\nB,sample2,Yes\n", encoding="utf-8")

        columns, rows = launcher.read_metadata(path)

        self.assertEqual(columns, ["Group", "SampleID", "Control"])
        self.assertEqual([row["SampleID"] for row in rows], ["sample1", "sample2"])

    def test_configuration_uses_selected_arbitrary_metadata_csv_path(self):
        app, _fastq_files = self.make_configuration_app()
        metadata_path = self.folder / "experiment1.csv"
        metadata_path.write_text(
            "Group,SampleID,Control\n"
            "A,A1,Negative\nA,A2,Sample\nA,A3,Mock\n"
            "B,B1,Sample\nB,B2,Sample\n",
            encoding="utf-8",
        )
        app.vars["metadata_file"].set(str(metadata_path))

        config = app.configuration()

        self.assertEqual(config["metadata_file"], str(metadata_path.resolve()).replace("\\", "/"))
        self.assertEqual(config["sample_id_map"]["A1"], "A1")

    def test_configuration_requires_fastq_folder(self):
        app, _files = self.make_configuration_app()
        app.vars["fastq_dir"].set("")
        with self.assertRaisesRegex(ValueError, "FASTQ folder"):
            app.configuration()

    def test_configuration_rejects_folder_without_fastq_files(self):
        app, _files = self.make_configuration_app()
        empty_folder = self.folder / "empty"
        empty_folder.mkdir()
        (empty_folder / "notes.txt").write_text("not sequencing data", encoding="utf-8")
        app.vars["fastq_dir"].set(str(empty_folder))
        with self.assertRaisesRegex(ValueError, "contains no FASTQ"):
            app.configuration()

    def test_configuration_extracts_fastq_zip_for_r_backend(self):
        app, fastq_files = self.make_configuration_app()
        archive = self.folder / "miseqsopdata.zip"
        with zipfile.ZipFile(archive, "w") as zipped:
            for path in fastq_files:
                zipped.write(path, f"miseqsopdata/reads/{Path(path).name}")
        app.vars["fastq_dir"].set(str(archive))

        config = app.configuration()

        extraction_root = (self.folder / "output" / "extracted_fastq").resolve()
        self.assertTrue(Path(config["fastq_dir"]).is_dir())
        self.assertTrue(Path(config["fastq_dir"]).is_relative_to(extraction_root))
        self.assertEqual(config["project_dir"], config["fastq_dir"])
        self.assertEqual(len(config["fastq_files"]), len(fastq_files))
        self.assertTrue(all(Path(path).is_file() for path in config["fastq_files"]))
        self.assertTrue(all(path.endswith(".fastq.gz") for path in config["fastq_files"]))

    def test_folder_zip_extraction_preserves_fastq_gz_compression(self):
        source = self.folder / "source"
        source.mkdir()
        archive = source / "miseqsopdata.zip"
        payloads = {
            "nested/reads/S1_R1.fastq.gz": b"compressed-r1",
            "nested/reads/S1_R2.fastq.gz": b"compressed-r2",
        }
        with zipfile.ZipFile(archive, "w") as zipped:
            for name, payload in payloads.items():
                zipped.writestr(name, payload)

        read_dir, fastq_paths = launcher.resolve_fastq_input(source, self.folder / "output")

        self.assertEqual(read_dir.name, "reads")
        self.assertEqual({Path(path).name for path in fastq_paths}, {
            "S1_R1.fastq.gz", "S1_R2.fastq.gz",
        })
        for path in fastq_paths:
            self.assertEqual(Path(path).read_bytes(), payloads[f"nested/reads/{Path(path).name}"])
        for name in ("S_R1.fastq", "S_R1.fastq.gz", "S_R1.fq", "S_R1.fq.gz"):
            self.assertTrue(launcher.is_fastq_path(name), name)

    def test_parent_folder_finds_supported_fastq_files_in_subfolders(self):
        source = self.folder / "project"
        reads = source / "nested" / "reads"
        reads.mkdir(parents=True)
        names = ("S_R1.fastq", "S_R2.fastq.gz", "T_1.fq", "T_2.fq.gz")
        for name in names:
            (reads / name).write_bytes(b"placeholder")

        read_dir, fastq_paths = launcher.resolve_fastq_input(source, self.folder / "output")

        self.assertEqual(read_dir, reads.resolve())
        self.assertEqual({Path(path).name for path in fastq_paths}, set(names))

    def test_folder_with_fastq_uses_reads_without_extracting_zip(self):
        source = self.folder / "source"
        source.mkdir()
        for name in ("S1_R1.fastq", "S1_R2.fastq"):
            (source / name).write_bytes(b"read")
        archive = source / "additional.zip"
        with zipfile.ZipFile(archive, "w") as zipped:
            zipped.writestr("S2_R1.fastq", b"archive read")
            zipped.writestr("S2_R2.fastq", b"archive read")

        read_dir, fastq_paths = launcher.resolve_fastq_input(source, self.folder / "output")

        self.assertEqual(read_dir, source.resolve())
        self.assertEqual({Path(path).name for path in fastq_paths}, {
            "S1_R1.fastq", "S1_R2.fastq",
        })
        self.assertFalse((self.folder / "output" / "extracted_fastq").exists())

    def test_configuration_maps_task_controls_to_existing_pipeline_parameters(self):
        app, _fastq_files = self.make_configuration_app()

        config = app.configuration()

        self.assertEqual(config["primer_method"], "cutadapt")
        self.assertTrue(config["automatic_truncation"])
        self.assertEqual(config["trunc_qmin"], 30)
        self.assertEqual(config["mergepairs_strategy"], "merge")
        self.assertEqual(config["mergepairs_min_overlap"], 12)
        self.assertEqual(config["max_sample_loss"], 0.10)
        self.assertEqual(config["plateau_gain"], 0.02)
        self.assertEqual(config["plateau_samples"], 0.90)
        self.assertTrue(config["run_decontam"])
        self.assertTrue(config["exclude_mock_controls"])
        self.assertEqual(config["negative_control_column"], "Control")
        self.assertEqual(config["negative_control_level"], "Negative")
        self.assertEqual(config["mock_control_level"], "Mock")
        self.assertNotIn("language", config)
        self.assertEqual(config["marker"], "16S")
        self.assertEqual(config["marker_profile"]["target_domain"], "Bacteria")
        self.assertEqual(config["reference_db_path"], str((self.folder / "training.fa.gz").resolve()).replace("\\", "/"))
        self.assertFalse(config["asv_length_filter_enabled"])
        self.assertIsNone(config["asv_min_len"])
        self.assertIsNone(config["asv_max_len"])
        self.assertEqual(config["min_rank1_abundance_percent"], 1.0)
        self.assertEqual(config["min_rank2_abundance_percent"], 2.0)
        self.assertEqual(config["figure_format"], "png")

    def test_configuration_accepts_independent_decimal_taxon_thresholds_and_format(self):
        app, _fastq_files = self.make_configuration_app()
        app.vars["min_rank1_abundance_percent"].set("0.5")
        app.vars["min_rank2_abundance_percent"].set("2.75")
        app.vars["figure_format"].set("PDF")

        config = app.configuration()

        self.assertEqual(config["min_rank1_abundance_percent"], 0.5)
        self.assertEqual(config["min_rank2_abundance_percent"], 2.75)
        self.assertEqual(config["figure_format"], "pdf")

    def test_configuration_accepts_zero_threshold_and_rejects_invalid_plot_settings(self):
        app, _fastq_files = self.make_configuration_app()
        app.vars["min_rank1_abundance_percent"].set("0")
        app.vars["min_rank2_abundance_percent"].set("0")
        config = app.configuration()
        self.assertEqual(config["min_rank1_abundance_percent"], 0)
        self.assertEqual(config["min_rank2_abundance_percent"], 0)

        app.vars["min_rank1_abundance_percent"].set("100.1")
        with self.assertRaisesRegex(ValueError, "must not exceed 100"):
            app.configuration()
        app.vars["min_rank1_abundance_percent"].set("1")
        app.vars["figure_format"].set("SVG")
        with self.assertRaisesRegex(ValueError, "Choose PNG or PDF"):
            app.configuration()


    def test_configuration_supports_length_primers_manual_truncation_and_concatenation(self):
        app, _fastq_files = self.make_configuration_app()
        app.vars["primer_mode"].value = "length"
        app.vars["truncation_mode"].value = "manual"
        app.vars["mergepairs_strategy"].value = "Concatenate reads"
        app.vars["mock_control_mode"].value = "None"
        app.vars["negative_control_mode"].value = "None"

        config = app.configuration()

        self.assertEqual(config["primer_method"], "trim_left")
        self.assertEqual(config["cutadapt"], "")
        self.assertFalse(config["automatic_truncation"])
        self.assertEqual(config["trunc_len"], [240, 200])
        self.assertEqual(config["mergepairs_strategy"], "concatenate")
        self.assertFalse(config["run_decontam"])
        self.assertFalse(config["exclude_mock_controls"])

    def test_duplicate_metadata_sample_is_rejected(self):
        path = self.folder / "metadata.csv"
        path.write_text("SampleID,Group\nsample,A\nsample,B\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "unique"):
            launcher.read_metadata(path)

    def test_malformed_csv_is_rejected(self):
        path = self.folder / "metadata.csv"
        path.write_text("SampleID,Group\nsample,A,unexpected\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "different number of fields"):
            launcher.read_metadata(path)

    def test_interface_messages_are_english_only(self):
        self.assertEqual(launcher.message("selected_files", count=3), "3 FASTQ files selected")
        self.assertEqual(launcher.message("clear_temp_cache"), "🧹 Clear Temporary Cache")
        with self.assertRaisesRegex(ValueError, "Select paired-end FASTQ"):
            launcher.pair_files([])

    def test_quality_profiles_are_limited_to_one_representative_sample(self):
        quality_script = (ROOT / "pipeline" / "02_quality_profiles.R").read_text(encoding="utf-8")
        self.assertEqual(quality_script.count("plotQualityProfile("), 2)
        self.assertIn("plotQualityProfile(fastq_info$fnF[1])", quality_script)
        self.assertIn("plotQualityProfile(fastq_info$fnR[1])", quality_script)

    def test_relative_abundance_uses_taxon_mean_and_separate_thresholds(self):
        abundance_script = (ROOT / "pipeline" / "15_relative_abundance.R").read_text(encoding="utf-8")
        self.assertIn("make_relative_abundance <- function(taxrank, threshold)", abundance_script)
        self.assertIn("relative_abundance_ranks[1]", abundance_script)
        self.assertIn("relative_abundance_ranks[2]", abundance_script)
        self.assertIn("mean_relative_abundance = mean(Abundance, na.rm = TRUE)", abundance_script)
        self.assertIn("mean_relative_abundance < threshold", abundance_script)
        self.assertNotIn("Phylum", abundance_script)
        self.assertNotIn("Genus", abundance_script)

    def test_all_figure_saving_uses_central_format_helper(self):
        helpers = (ROOT / "pipeline" / "01_helpers.R").read_text(encoding="utf-8")
        self.assertIn('paste0(filename, ".", figure_format)', helpers)
        self.assertNotIn('paste0(filename, ".pdf")', helpers)
        self.assertNotIn('paste0(filename, ".png")', helpers)

    def test_pipeline_uses_run_local_output_folders(self):
        helpers = (ROOT / "pipeline" / "01_helpers.R").read_text(encoding="utf-8")
        launcher_r = (ROOT / "run_pipeline.R").read_text(encoding="utf-8")
        denoise = (ROOT / "pipeline" / "06_denoise_merge.R").read_text(encoding="utf-8")
        rarefaction = (ROOT / "pipeline" / "12_rarefaction.R").read_text(encoding="utf-8")

        self.assertNotIn("project_cache_dir", launcher_r + helpers)
        self.assertNotIn("cfg$force_recompute", launcher_r)
        self.assertIn('dir_results <- "."', helpers)
        for folder in ("rds", "tables", "figures", "logs"):
            self.assertIn(f'file.path(dir_results, "{folder}")', helpers)
        self.assertIn("get_derep_fwd <- local({", denoise)
        self.assertIn('file.path(dir_rds, "rarefaction_depth.rds")', rarefaction)

    def test_run_local_rds_intermediates_are_not_deleted(self):
        launcher_r = (ROOT / "run_pipeline.R").read_text(encoding="utf-8")

        self.assertNotIn("unlink(", launcher_r)
        self.assertNotIn("save_checkpoints", launcher_r)

    def test_pipeline_has_no_cross_run_result_reuse(self):
        for path in (ROOT / "pipeline").glob("[0-9][0-9]_*.R"):
            script = path.read_text(encoding="utf-8")
            self.assertNotIn("load_or_compute", script, path.name)
            self.assertNotIn("fingerprint", script, path.name)

    def test_filtering_is_sequential_and_reports_sample_timestamps(self):
        filter_script = (ROOT / "pipeline" / "04_filter_trim.R").read_text(encoding="utf-8")
        self.assertIn("dada2::filterAndTrim(", filter_script)
        self.assertIn("multithread = FALSE", filter_script)
        self.assertIn("verbose     = TRUE", filter_script)
        self.assertIn("Filtering sample %d/%d: %s", filter_script)
        self.assertIn("Finished sample %d/%d: %s", filter_script)
        self.assertNotIn("parallel::makeCluster", filter_script)
        self.assertNotIn("parallel::parLapply", filter_script)

    def test_main_interface_uses_task_oriented_english_labels_and_fixed_log(self):
        ui_source = inspect.getsource(launcher.PipelineApp.build_ui)
        self.assertNotIn('text="Threads"', ui_source)
        self.assertNotIn('self.vars["threads"]', ui_source)

        for implementation_detail in ("Cutadapt", "trimLeft", "Toggle Console", "Language"):
            self.assertNotIn(implementation_detail, ui_source)
        for page in ("Project", "Read processing", "Downstream analysis"):
            self.assertIn(f'scrollable_tab("{page}")', ui_source)
        self.assertIn('("FASTQ folder:", "fastq_dir", None)', ui_source)
        self.assertNotIn("Choose ZIP", ui_source)
        self.assertNotIn("pick_fastq_zip", ui_source)
        picker_source = inspect.getsource(launcher.PipelineApp.choose_fastq_input)
        self.assertIn("Select FASTQ folder", picker_source)
        self.assertIn("_windows_fastq_input_dialog(", picker_source)
        self.assertNotIn("ZIP archive (cancel", picker_source)
        self.assertNotIn("tk.Toplevel(", picker_source)
        self.assertNotIn("Treeview", picker_source)
        self.assertNotIn("Listbox", picker_source)
        self.assertIn('text="Progress"', ui_source)
        self.assertNotIn("console_frame.grid_remove()", ui_source)

    def test_orphan_cleanup_removes_only_dead_biome_temp_sessions(self):
        temp_root = self.folder / "BIOME" / "temp"
        temp_root.mkdir(parents=True)
        stale = temp_root / "session_stale"
        stale.mkdir()
        (stale / ".owner").write_text("2147483647", encoding="ascii")
        (stale / "RtmpABC" / "worker.tmp").parent.mkdir(parents=True)
        (stale / "RtmpABC" / "worker.tmp").write_bytes(b"temporary")
        active = temp_root / "session_active"
        active.mkdir()
        (active / ".owner").write_text(str(launcher.os.getpid()), encoding="ascii")
        (active / "RtmpXYZ").mkdir()
        (active / "RtmpXYZ" / "keep.tmp").write_bytes(b"active")
        results = self.folder / "results"
        (results / "filtered").mkdir(parents=True)
        output_fastq = results / "filtered" / "sample.fastq.gz"
        output_fastq.write_bytes(b"user result")

        def process_is_running(pid):
            return pid == launcher.os.getpid()

        with patch.object(launcher, "_process_is_running", side_effect=process_is_running):
            freed = launcher.clear_orphaned_temp_dirs(temp_root)

        self.assertGreater(freed, 0)
        self.assertFalse(stale.exists())
        self.assertTrue((active / "RtmpXYZ" / "keep.tmp").exists())
        self.assertEqual(output_fastq.read_bytes(), b"user result")

    def test_temp_cleanup_rejects_directories_outside_app_temp_root(self):
        temp_root = self.folder / "BIOME" / "temp"
        results = self.folder / "results"
        unowned = results / "session_keep"
        unowned.mkdir(parents=True)
        (unowned / ".owner").write_text(str(launcher.os.getpid()), encoding="ascii")
        preserved = unowned / "important.csv"
        preserved.write_text("user result", encoding="utf-8")

        with self.assertRaisesRegex(ValueError, "outside the BIOME temp folder"):
            launcher._remove_owned_temp_dir(unowned, root=temp_root)

        self.assertEqual(preserved.read_text(encoding="utf-8"), "user result")

    def test_temp_cleanup_preserves_rds_files(self):
        temp_root = self.folder / "BIOME" / "temp"
        session = temp_root / "session_rds"
        nested = session / "RtmpABC" / "analysis"
        nested.mkdir(parents=True)
        (session / ".owner").write_text(str(launcher.os.getpid()), encoding="ascii")
        temporary_file = nested / "scratch.tmp"
        temporary_file.write_bytes(b"temporary")
        rds_file = nested / "intermediate.RDS"
        rds_file.write_bytes(b"keep this analysis object")

        freed = launcher._remove_owned_temp_dir(session, root=temp_root)

        self.assertGreater(freed, 0)
        self.assertFalse(temporary_file.exists())
        self.assertEqual(rds_file.read_bytes(), b"keep this analysis object")
        self.assertTrue(session.exists())
        self.assertTrue((session / ".owner").exists())

    def test_path_picker_routes_files_and_output_folder_to_dialogs(self):
        app = launcher.PipelineApp.__new__(launcher.PipelineApp)
        app.root = object()
        metadata_var = Mock()
        output_var = Mock()
        app.vars = {"metadata_file": metadata_var, "output_dir": output_var}
        app.refresh_metadata = Mock()
        csv_types = [("CSV", "*.csv"), ("All files", "*")]
        with patch.object(launcher.filedialog, "askopenfilename", return_value="metadata.csv") as open_file:
            app.pick_path("metadata_file", csv_types)
        open_file.assert_called_once_with(
            parent=app.root,
            title="Select Metadata CSV",
            filetypes=csv_types,
        )
        metadata_var.set.assert_called_once_with("metadata.csv")
        app.refresh_metadata.assert_called_once_with()

        with patch.object(launcher.filedialog, "askdirectory", return_value="results") as open_folder:
            app.pick_path("output_dir", None)
        open_folder.assert_called_once_with(parent=app.root, title="Select output folder")
        output_var.set.assert_called_once_with("results")

        fastq_var = Mock()
        app.vars["fastq_dir"] = fastq_var
        app.refresh_fastq_label = Mock()
        with patch.object(app, "choose_fastq_input", return_value="C:\\reads.zip") as choose_input:
            app.pick_path("fastq_dir", None)
        choose_input.assert_called_once_with()
        fastq_var.set.assert_called_once_with("C:\\reads.zip")
        app.refresh_fastq_label.assert_called_once_with()

    def test_fastq_picker_uses_one_windows_shell_dialog(self):
        app = launcher.PipelineApp.__new__(launcher.PipelineApp)
        app.root = object()
        app.vars = {"fastq_dir": Value(str(self.folder))}
        with patch.object(launcher.sys, "platform", "win32"), patch.object(
            launcher, "_windows_fastq_input_dialog", return_value="reads.zip"
        ) as shell_dialog, patch.object(launcher.filedialog, "askopenfilename") as open_file, patch.object(
            launcher.filedialog, "askdirectory"
        ) as open_folder:
            selected = app.choose_fastq_input()
        self.assertEqual(selected, "reads.zip")
        shell_dialog.assert_called_once_with(
            app.root,
            "Select FASTQ folder",
            str(self.folder.resolve()),
        )
        open_file.assert_not_called()
        open_folder.assert_not_called()

    def test_windows_fastq_picker_initializes_display_name_pointer(self):
        shell32 = Mock()
        shell32.SHBrowseForFolderW.return_value = None
        windll = SimpleNamespace(user32=Mock(), shell32=shell32, ole32=Mock())
        windll.ole32.CoInitializeEx.return_value = 0
        with patch.object(launcher.ctypes, "windll", windll, create=True), patch.object(
            launcher.ctypes, "WINFUNCTYPE", launcher.ctypes.CFUNCTYPE, create=True
        ):
            selected = launcher._windows_fastq_input_dialog(
                Mock(winfo_id=Mock(return_value=123)),
                "Select FASTQ folder",
                str(self.folder),
            )

        self.assertEqual(selected, "")
        shell32.SHBrowseForFolderW.assert_called_once()
        browse_info = shell32.SHBrowseForFolderW.call_args.args[0]._obj
        self.assertEqual(browse_info.pszDisplayName, "")
        windll.ole32.CoUninitialize.assert_called_once_with()

    def test_fastq_picker_uses_non_windows_folder_dialog(self):
        app = launcher.PipelineApp.__new__(launcher.PipelineApp)
        app.root = object()
        app.vars = {"fastq_dir": Value(str(self.folder))}
        with patch.object(launcher.sys, "platform", "linux"), patch.object(
            launcher.filedialog, "askdirectory", return_value=str(self.folder)
        ) as open_folder:
            selected = app.choose_fastq_input()
        self.assertEqual(selected, str(self.folder))
        open_folder.assert_called_once_with(
            parent=app.root, title="Select FASTQ folder",
            initialdir=str(self.folder.resolve()), mustexist=True,
        )

    def test_fastq_input_count_includes_nested_files_and_zip_members(self):
        nested = self.folder / "nested"
        nested.mkdir()
        (self.folder / "sample_R1.fastq.gz").write_bytes(b"read")
        (nested / "sample_R2.fq").write_bytes(b"read")
        archive = nested / "reads.zip"
        with zipfile.ZipFile(archive, "w") as zipped:
            zipped.writestr("archive/sample_R1.fastq", b"read")
        self.assertEqual(launcher.count_fastq_input_files(self.folder), 3)
        self.assertEqual(launcher.count_fastq_input_files(archive), 1)

    def test_fastq_label_count_runs_in_background_and_updates_via_poll(self):
        app = launcher.PipelineApp.__new__(launcher.PipelineApp)
        app.vars = {"fastq_dir": Value(str(self.folder))}
        app.fastq_label = Widget()
        app.events = queue.Queue()
        app.root = Mock()
        (self.folder / "sample_R1.fastq").write_bytes(b"read")

        with patch.object(launcher.threading, "Thread") as thread:
            app.refresh_fastq_label()
            self.assertEqual(app.fastq_label.options["text"], "Counting FASTQ files…")
            thread.return_value.start.assert_called_once_with()
            worker = thread.call_args.kwargs["target"]
            args = thread.call_args.kwargs["args"]
            worker(*args)

        app.poll()
        self.assertEqual(app.fastq_label.options["text"], "1 FASTQ files")
        app.root.after.assert_called_once_with(100, app.poll)

    def test_path_picker_surfaces_tk_dialog_errors(self):
        app = launcher.PipelineApp.__new__(launcher.PipelineApp)
        app.root = object()
        app.vars = {"metadata_file": Mock()}
        with patch.object(
            launcher.filedialog,
            "askopenfilename",
            side_effect=launcher.tk.TclError("dialog unavailable"),
        ), patch.object(launcher.messagebox, "showerror") as showerror:
            app.pick_path("metadata_file", [("CSV", "*.csv")])
        showerror.assert_called_once()
        self.assertIn("dialog unavailable", showerror.call_args.args[1])

    def test_startup_failure_is_logged_and_shown(self):
        log_dir = self.folder / "settings"
        try:
            raise RuntimeError("startup failed")
        except RuntimeError as exc:
            with patch.object(launcher, "SETTINGS_DIR", log_dir), patch.object(
                launcher.messagebox, "showerror"
            ) as showerror:
                launcher.report_startup_failure(type(exc), exc, exc.__traceback__)

        log_path = log_dir / "startup.log"
        self.assertIn("RuntimeError: startup failed", log_path.read_text(encoding="utf-8"))
        showerror.assert_called_once()
        self.assertIn(str(log_path), showerror.call_args.args[1])

    def test_uncaught_import_failure_is_logged_and_shown_without_tk(self):
        log_dir = self.folder / "import-settings"
        try:
            raise ImportError("startup import failed")
        except ImportError as exc:
            with patch.object(launcher, "SETTINGS_DIR", log_dir), patch.object(
                launcher, "_show_native_startup_error"
            ) as show_error:
                launcher._uncaught_startup_exception(type(exc), exc, exc.__traceback__)

        log_path = log_dir / "startup.log"
        self.assertIn("ImportError: startup import failed", log_path.read_text(encoding="utf-8"))
        show_error.assert_called_once()
        self.assertIn(str(log_path), show_error.call_args.args[0])

    def test_uncaught_background_thread_failure_is_logged_and_shown(self):
        log_dir = self.folder / "thread-settings"
        try:
            raise RuntimeError("background setup failed")
        except RuntimeError as exc:
            args = SimpleNamespace(
                exc_type=type(exc),
                exc_value=exc,
                exc_traceback=exc.__traceback__,
            )
            with patch.object(launcher, "SETTINGS_DIR", log_dir), patch.object(
                launcher, "_show_native_startup_error"
            ) as show_error:
                launcher.threading.excepthook(args)

        log_path = log_dir / "startup.log"
        self.assertIn("RuntimeError: background setup failed", log_path.read_text(encoding="utf-8"))
        show_error.assert_called_once()
        self.assertIn(str(log_path), show_error.call_args.args[0])

    def test_metadata_browse_button_invokes_path_picker(self):
        app = launcher.PipelineApp.__new__(launcher.PipelineApp)
        app.root = object()
        app.controls = []
        app.vars = {"metadata_file": Mock()}
        app.refresh_metadata = Mock()
        button = Mock()
        app.button = Mock(return_value=button)
        csv_types = [("CSV", "*.csv")]
        with patch.object(launcher.ttk, "Label", return_value=Mock()), patch.object(
            launcher.ttk, "Entry", return_value=Mock()
        ), patch.object(launcher.filedialog, "askopenfilename", return_value="metadata.csv") as open_file:
            app.path_row(object(), 0, "Metadata CSV", "metadata_file", csv_types)
            app.button.call_args.args[2]()
        open_file.assert_called_once()
        app.vars["metadata_file"].set.assert_called_once_with("metadata.csv")
        app.refresh_metadata.assert_called_once_with()

    def test_worker_uses_one_r_process_and_preserves_failure_status(self):
        app = launcher.PipelineApp.__new__(launcher.PipelineApp)
        app.events = queue.Queue()
        app.cancelled = threading.Event()
        cfg = {"rscript": "C:/R path/Rscript.exe", "run_dir": self.folder.as_posix()}
        app.temp_dir = launcher.create_app_temp_dir(self.folder / "app_temp")
        scratch_file = app.temp_dir / "RtmpABC" / "buffers.tmp"
        scratch_file.parent.mkdir()
        scratch_file.write_bytes(b"temporary R data")
        fake_process = unittest.mock.Mock()
        sample_progress = "[2026-10-22 20:00:00] Filtering sample 1/20: sample1\n"
        fake_process.stdout = iter(["[STEP 03]\n", sample_progress, "Cutadapt failed\n"])
        fake_process.wait.return_value = 1
        with patch.object(launcher.subprocess, "Popen", return_value=fake_process) as popen:
            app.worker(cfg, check_only=False)
        popen.assert_called_once()
        command = popen.call_args.args[0]
        self.assertEqual(command[0], cfg["rscript"])
        self.assertEqual(command[-2], str(ROOT / "run_pipeline.R"))
        self.assertEqual(popen.call_args.kwargs["env"]["TMPDIR"], str(scratch_file.parents[1]))
        self.assertFalse(scratch_file.parents[1].exists())
        self.assertIsNone(app.temp_dir)
        self.assertEqual(json.loads((self.folder / "logs" / "status.json").read_text())["exit_code"], 1)
        pipeline_log = (self.folder / "logs" / "pipeline.log").read_text()
        self.assertIn("Starting R pipeline.", pipeline_log)
        self.assertIn(sample_progress, pipeline_log)
        self.assertIn("Cutadapt failed", pipeline_log)
        emitted_logs = []
        while not app.events.empty():
            event, line = app.events.get_nowait()
            if event == "log":
                emitted_logs.append(line)
        self.assertIn(sample_progress, emitted_logs)

    def test_start_writes_input_path_csv(self):
        start_source = inspect.getsource(launcher.PipelineApp.start)
        self.assertIn('"input_path.csv"', start_source)
        self.assertNotIn("manifest", start_source)

    def test_worker_notifies_once_and_logs_all_large_plot_warnings(self):
        app = launcher.PipelineApp.__new__(launcher.PipelineApp)
        app.events = queue.Queue()
        app.cancelled = threading.Event()
        app.temp_dir = None
        lines = ["BIOME_PLOT_WARNING: Phylum is crowded.\n",
                 "BIOME_PLOT_WARNING: Genus is crowded.\n"]
        process = Mock(stdout=iter(lines))
        process.wait.return_value = 0
        cfg = {"rscript": "Rscript", "run_dir": self.folder.as_posix()}
        with patch.object(launcher.subprocess, "Popen", return_value=process):
            app.worker(cfg, check_only=False)
        events = []
        while not app.events.empty():
            events.append(app.events.get_nowait())
        self.assertEqual([data for kind, data in events if kind == "plot_warning"],
                         ["Phylum is crowded."])
        log = (self.folder / "logs" / "pipeline.log").read_text()
        for line in lines:
            self.assertIn(line, log)

    def test_start_creates_only_four_run_folders(self):
        app, _files = self.make_configuration_app()
        config = app.configuration()
        app.configuration = Mock(return_value=config)
        app.root = Mock()
        app.busy = False
        app.temp_dir = None
        app.cancelled = threading.Event()
        app.saved = {}
        app.status = Value("")
        app.append_log = Mock()
        app.set_running_controls = Mock()
        create_temp = launcher.create_app_temp_dir
        with patch.object(launcher, "save_settings"), patch.object(
            launcher.threading, "Thread"
        ), patch.object(launcher, "create_app_temp_dir",
                        side_effect=lambda: create_temp(self.folder / "temp")), patch.object(
            launcher.messagebox, "showerror", side_effect=AssertionError("Run startup failed")
        ):
            app.start()
        self.assertEqual({path.name for path in app.last_run.iterdir()},
                         {"rds", "tables", "figures", "logs"})
        self.assertTrue((app.last_run / "tables" / "input_path.csv").is_file())
        self.assertTrue((app.last_run / "logs" / "BIOME_config.txt").is_file())
        self.assertFalse(list((app.last_run / "logs").glob("*config*.json")))
        self.assertEqual(Path(config["config_path"]).parent, app.temp_dir)
        self.assertEqual(json.loads(Path(config["config_path"]).read_text(encoding="utf-8"))["marker"], "16S")
        app.set_running_controls.assert_called_once_with(True)
        launcher._remove_owned_temp_dir(app.temp_dir, root=app.temp_dir.parent)

    def test_configuration_report_is_readable_and_keeps_user_values(self):
        path = self.folder / "BIOME_config.txt"
        config = {"marker": "18S", "threads": 7, "max_sample_loss": 0.1,
                  "min_rank1_abundance_percent": 1.25, "min_rank2_abundance_percent": 2.75,
                  "asv_length_filter_enabled": False, "reference_group": "Control",
                  "steps": {"step_1": True, "step_2": True, "step_3": False},
                  "fastq_files": [str(self.folder / "001_R1.fastq.gz")],
                  "config_path": "private/runtime.json"}
        launcher.write_configuration_report(config, path)
        text = path.read_text(encoding="utf-8")
        self.assertIn("BIOME - Run configuration", text)
        self.assertIn("Marker: 18S", text)
        self.assertIn("Max sample loss: 10%", text)
        self.assertIn("First rank collapse threshold: 1.25%", text)
        self.assertIn("Second rank collapse threshold: 2.75%", text)
        self.assertIn("Asv length filter enabled: No", text)
        self.assertIn("Downstream analysis: No", text)
        self.assertIn(config["fastq_files"][0], text)
        self.assertNotIn("private/runtime.json", text)

    def test_legacy_settings_are_read_and_new_settings_take_precedence(self):
        current = self.folder / "BIOME" / "settings.json"
        legacy = self.folder / "old-settings.json"
        fallback = self.folder / "last_selection.json"
        legacy.write_text(json.dumps({"rscript": "legacy Rscript"}), encoding="utf-8")
        with patch.object(launcher, "SETTINGS", current), patch.object(
            launcher, "SETTINGS_DIR", current.parent
        ), patch.object(launcher, "LEGACY_APP_SETTINGS", legacy), patch.object(
            launcher, "LEGACY_SETTINGS", fallback
        ):
            self.assertEqual(launcher.load_settings()["rscript"], "legacy Rscript")
            launcher.save_settings({"rscript": "BIOME Rscript"})
            self.assertEqual(launcher.load_settings()["rscript"], "BIOME Rscript")
        self.assertTrue(legacy.is_file())


if __name__ == "__main__":
    unittest.main(verbosity=2)
