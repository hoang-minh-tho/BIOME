# Required settings must fail before the launcher loads pipeline parameters.
launcher <- as.list(parse("run_pipeline.R"))
definition <- Filter(function(expr) {
  is.call(expr) && identical(expr[[1]], as.name("<-")) &&
    identical(expr[[2]], as.name("validate_required_settings"))
}, launcher)[[1]]
eval(definition)

for (marker in c("16S", "18S")) {
  for (method in c("cutadapt", "trim_left")) {
    cfg <- list(
      marker = marker, primer_method = method, marker_profiles_path = "profiles.json",
      pipeline_dir = "pipeline", run_dir = "run", fastq_dir = "reads",
      metadata_file = "metadata.csv", trim_left_fwd = 0L, trim_left_rev = 0L,
      min_len = 50L, mergepairs_min_overlap = 12L, mergepairs_max_mismatch = 0L,
      mergepairs_strategy = "merge", chimera_method = "consensus", trunc_qmin = 30L,
      threads = 2L, max_n = 0L, trunc_q = 2L, max_sample_loss = 0.1,
      plateau_gain = 0.02, plateau_samples = 0.9, read_step = 1000L,
      cutadapt_error_rate = 0.1, ordination = "PCoA", permanova_permutations = 999L,
      differential_method = "ANCOM-BC", rarefaction_mode = "custom",
      custom_rarefaction_depth = 100L, figure_format = "png",
      min_rank1_abundance_percent = 1, min_rank2_abundance_percent = 2,
      pooling = "pseudo", asv_length_filter_enabled = TRUE,
      asv_min_len = 300L, asv_max_len = 500L,
      disable_truncation = FALSE, automatic_truncation = TRUE, remove_phix = TRUE,
      discard_untrimmed = FALSE, run_betadisper = FALSE, run_permanova = FALSE,
      run_decontam = FALSE, exclude_controls = FALSE, exclude_mock_controls = FALSE,
      max_ee = c(2, 2), trunc_len = c(0, 0),
      negative_control_column = "", negative_control_level = "",
      mock_control_column = "", mock_control_level = "", beta_distances = character(),
      negative_control_prefixes = character(), positive_control_prefixes = character(),
      fastq_files = c("sample_R1.fastq.gz", "sample_R2.fastq.gz"),
      sample_id_map = c(sample = "sample"),
      steps = list(step_1 = TRUE, step_2 = TRUE, step_3 = TRUE),
      reference_db_path = "reference.fa.gz", group_column = "Group", reference_group = "A"
    )
    if (method == "cutadapt") {
      cfg$primer_fwd <- "ACGT"
      cfg$primer_rev <- "TGCA"
      cfg$cutadapt <- "cutadapt.exe"
    }
    validate_required_settings(cfg)
    inactive <- cfg
    inactive$asv_length_filter_enabled <- FALSE
    inactive$asv_min_len <- NULL
    inactive$asv_max_len <- NULL
    inactive$rarefaction_mode <- "min"
    inactive$custom_rarefaction_depth <- NULL
    validate_required_settings(inactive)
    for (name in names(cfg)) {
      for (invalid in list(NULL)) {
        broken <- cfg
        broken[name] <- list(invalid)
        error <- tryCatch(validate_required_settings(broken), error = identity)
        stopifnot(inherits(error, "error"),
                  grepl(name, conditionMessage(error), fixed = TRUE))
      }
    }
  }
}
cat("Missing required settings rejected; explicit FALSE, zero and inactive options accepted for 16S/18S.\n")
