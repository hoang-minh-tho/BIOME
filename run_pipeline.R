# Called once from Python; all steps share the same R environment.
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("Usage: Rscript run_pipeline.R config.json [--check]")
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Missing jsonlite. Run install_packages.R in the selected R installation.")
}
cfg <- jsonlite::fromJSON(args[1], simplifyVector = TRUE)
marker <- if (is.null(cfg$marker) || !nzchar(cfg$marker)) "16S" else toupper(cfg$marker)
marker_profiles_path <- if (!is.null(cfg$marker_profiles_path) &&
                            nzchar(cfg$marker_profiles_path)) {
  cfg$marker_profiles_path
} else {
  file.path(cfg$pipeline_dir, "marker_profiles.json")
}
if (!file.exists(marker_profiles_path)) {
  stop("Marker profile file not found: ", marker_profiles_path)
}
marker_profiles <- jsonlite::fromJSON(marker_profiles_path, simplifyVector = FALSE)
if (!(marker %in% names(marker_profiles))) {
  stop(
    "Unsupported marker: ", marker, ". Configured markers: ",
    paste(names(marker_profiles), collapse = ", ")
  )
}
marker_profile <- marker_profiles[[marker]]
taxonomy_ranks <- unlist(marker_profile$taxonomy_ranks, use.names = FALSE)
target_domain_rank <- marker_profile$target_domain_rank
target_domain <- marker_profile$target_domain
relative_abundance_ranks <- unlist(
  marker_profile$relative_abundance_ranks, use.names = FALSE
)
contaminant_rank <- marker_profile$contaminant_rank
contaminant_taxa <- unlist(marker_profile$contaminants, use.names = FALSE)
organelle_taxa <- as.character(unlist(marker_profile$organelle_taxa, use.names = FALSE))
run_stage_max <- max(which(vapply(cfg$steps, isTRUE, logical(1))))
step_files <- sort(list.files(cfg$pipeline_dir, pattern = "^[0-9]{2}_.*\\.R$",
                              full.names = TRUE))
if (length(step_files) != 17L ||
    !identical(substr(basename(step_files), 1, 2), sprintf("%02d", 0:16))) {
  stop("The pipeline directory must contain exactly one script per step, 00 to 16.")
}
source(step_files[1], local = .GlobalEnv, encoding = "UTF-8")
if (!is.null(cfg$min_rank1_abundance_percent)) {
  min_rank1_abundance <- as.numeric(cfg$min_rank1_abundance_percent) / 100
} else if (!is.null(cfg$min_phylum_abundance_percent)) {
  min_rank1_abundance <- as.numeric(cfg$min_phylum_abundance_percent) / 100
}
if (!is.null(cfg$min_rank2_abundance_percent)) {
  min_rank2_abundance <- as.numeric(cfg$min_rank2_abundance_percent) / 100
} else if (!is.null(cfg$min_genus_abundance_percent)) {
  min_rank2_abundance <- as.numeric(cfg$min_genus_abundance_percent) / 100
}
figure_format <- if (is.null(cfg$figure_format)) "png" else tolower(cfg$figure_format)
if (!figure_format %in% c("png", "pdf")) {
  stop("Figure format must be PNG or PDF.")
}
fastq_dir <- if (is.null(cfg$fastq_dir)) cfg$project_dir else cfg$fastq_dir
project_dir <- fastq_dir
metadata_file <- cfg$metadata_file
reference_db_path <- if (!is.null(cfg$reference_db_path) &&
                         nzchar(cfg$reference_db_path)) {
  cfg$reference_db_path
} else if (!is.null(cfg$silva_ref_path) && nzchar(cfg$silva_ref_path)) {
  cfg$silva_ref_path
} else {
  reference_db_path
}
group_column <- cfg$group_column
reference_group <- cfg$reference_group
primer_method <- if (is.null(cfg$primer_method)) "cutadapt" else cfg$primer_method
if (!primer_method %in% c("cutadapt", "trim_left")) {
  stop("Unknown primer_method: ", primer_method)
}
primer_fwd <- cfg$primer_fwd
primer_rev <- cfg$primer_rev
trim_values <- c(
  if (is.null(cfg$trim_left_fwd)) 19 else cfg$trim_left_fwd,
  if (is.null(cfg$trim_left_rev)) 20 else cfg$trim_left_rev
)
trim_values <- suppressWarnings(as.numeric(trim_values))
if (anyNA(trim_values) || any(!is.finite(trim_values)) ||
    any(trim_values < 0) || any(trim_values != floor(trim_values))) {
  stop("trim_left_fwd and trim_left_rev must be non-negative integers.")
}
trim_left_fwd <- as.integer(trim_values[1])
trim_left_rev <- as.integer(trim_values[2])
asv_length_filter_enabled <- if (is.null(cfg$asv_length_filter_enabled)) {
  if (is.null(cfg$marker)) TRUE else
    isTRUE(marker_profile$legacy_asv_length_filter$enabled)
} else {
  isTRUE(cfg$asv_length_filter_enabled)
}
asv_min_len <- if (asv_length_filter_enabled) {
  if (is.null(cfg$asv_min_len)) marker_profile$legacy_asv_length_filter$minimum
  else as.integer(cfg$asv_min_len)
} else {
  NULL
}
asv_max_len <- if (asv_length_filter_enabled) {
  if (is.null(cfg$asv_max_len)) marker_profile$legacy_asv_length_filter$maximum
  else as.integer(cfg$asv_max_len)
} else {
  NULL
}
if (asv_length_filter_enabled &&
    (is.null(asv_min_len) || is.null(asv_max_len) ||
     anyNA(c(asv_min_len, asv_max_len)) ||
     asv_min_len <= 0L || asv_min_len >= asv_max_len)) {
  stop("ASV length filtering requires a positive minimum length smaller than the maximum length.")
}
trunc_qmin <- if (is.null(cfg$trunc_qmin)) 30L else as.integer(cfg$trunc_qmin)
N_THREADS <- as.integer(cfg$threads)
message("[CPU] Automatically selected ", N_THREADS,
        " thread(s) for supported operations.")
trunc_len <- as.integer(unlist(cfg$trunc_len, use.names = FALSE))
disable_truncation <- isTRUE(cfg$disable_truncation)
automatic_truncation <- isTRUE(cfg$automatic_truncation)
max_ee <- as.numeric(unlist(cfg$max_ee, use.names = FALSE))
max_n <- as.integer(cfg$max_n)
trunc_q <- as.integer(cfg$trunc_q)
min_len <- as.integer(cfg$min_len)
chimera_method <- if (is.null(cfg$chimera_method)) chimera_method else cfg$chimera_method
remove_phix <- if (is.null(cfg$remove_phix)) remove_phix else isTRUE(cfg$remove_phix)
mergepairs_strategy <- if (is.null(cfg$mergepairs_strategy)) {
  mergepairs_strategy
} else {
  cfg$mergepairs_strategy
}
mergepairs_min_overlap <- if (is.null(cfg$mergepairs_min_overlap)) {
  mergepairs_min_overlap
} else {
  as.integer(cfg$mergepairs_min_overlap)
}
mergepairs_max_mismatch <- if (is.null(cfg$mergepairs_max_mismatch)) {
  mergepairs_max_mismatch
} else {
  as.integer(cfg$mergepairs_max_mismatch)
}
max_sample_loss <- if (is.null(cfg$max_sample_loss)) max_sample_loss else as.numeric(cfg$max_sample_loss)
plateau_gain <- if (is.null(cfg$plateau_gain)) plateau_gain else as.numeric(cfg$plateau_gain)
plateau_samples <- if (is.null(cfg$plateau_samples)) plateau_samples else as.numeric(cfg$plateau_samples)
read_step <- if (is.null(cfg$read_step)) read_step else as.integer(cfg$read_step)
pooling <- if (identical(cfg$pooling, "pseudo")) "pseudo" else isTRUE(cfg$pooling)
cutadapt_error_rate <- if (is.null(cfg$cutadapt_error_rate)) 0.1 else as.numeric(cfg$cutadapt_error_rate)
discard_untrimmed <- isTRUE(cfg$discard_untrimmed)
run_betadisper <- isTRUE(cfg$run_betadisper)
run_permanova <- if (is.null(cfg$run_permanova)) TRUE else isTRUE(cfg$run_permanova)
ordination_method <- cfg$ordination
permanova_permutations <- as.integer(cfg$permanova_permutations)
beta_distances <- unlist(cfg$beta_distances, use.names = FALSE)
differential_method <- cfg$differential_method
rarefaction_mode <- cfg$rarefaction_mode
custom_rarefaction_depth <- cfg$custom_rarefaction_depth
CUTADAPT_COMMAND <- if (is.null(cfg$cutadapt)) "" else cfg$cutadapt
pipeline_input_files <- unlist(cfg$fastq_files, use.names = FALSE)
pipeline_sample_id_map <- unlist(cfg$sample_id_map, use.names = TRUE)
if (!length(pipeline_sample_id_map)) {
  stop("The run configuration is missing the FASTQ-to-metadata SampleID mapping.")
}
negative_control_column <- cfg$negative_control_column
negative_control_level <- cfg$negative_control_level
mock_control_column <- cfg$mock_control_column
mock_control_level <- cfg$mock_control_level
negative_control_prefixes <- unlist(cfg$negative_control_prefixes, use.names = FALSE)
positive_control_prefixes <- unlist(cfg$positive_control_prefixes, use.names = FALSE)
run_decontam <- isTRUE(cfg$run_decontam)
exclude_controls <- isTRUE(cfg$exclude_controls)
exclude_mock_controls <- isTRUE(cfg$exclude_mock_controls)
is_control <- function(sample_ids, prefixes) {
  if (!length(prefixes)) return(rep(FALSE, length(sample_ids)))
  vapply(sample_ids, function(sample_id) {
    any(startsWith(tolower(sample_id), tolower(prefixes)))
  }, logical(1))
}
set.seed(SEED)

required <- c("tidyverse", "dada2", "ShortRead", "Biostrings", "phyloseq",
              "DECIPHER", "phangorn", "vegan", "Polychrome", "rlang")
if (run_stage_max >= 3L && differential_method == "ANCOM-BC") required <- c(required, "ANCOMBC")
if (run_stage_max >= 3L && differential_method == "DESeq2") required <- c(required, "DESeq2")
if (run_stage_max >= 3L && differential_method == "ALDEx2") required <- c(required, "ALDEx2")
if (run_decontam) required <- c(required, "decontam")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
message("R: ", R.version.string)
message("R libraries: ", paste(.libPaths(), collapse = "; "))
if (length(missing)) {
  stop("Missing R packages: ", paste(missing, collapse = ", "),
       ". Run install_packages.R in this R installation.")
}
if (!all(c("hash", "hash_file") %in% getNamespaceExports("rlang"))) {
  stop("Update rlang: this pipeline needs rlang::hash and rlang::hash_file.")
}
if (primer_method == "cutadapt") {
  if (!file.exists(CUTADAPT_COMMAND) && Sys.which(CUTADAPT_COMMAND) == "") {
    stop("The primer-removal component is missing. Run dependency setup and try again.")
  }
  version <- system2(CUTADAPT_COMMAND, "--version", stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(version, "status")) && attr(version, "status") != 0L) {
    stop("The primer-removal component could not be started: ", paste(version, collapse = "\n"))
  }
  message("Primer-removal component is ready.")
} else {
  message("Primer removal: removing ", trim_left_fwd, " forward and ",
          trim_left_rev, " reverse bases.")
}
if (!file.exists(metadata_file) ||
    (run_stage_max >= 2L && !file.exists(reference_db_path)) ||
    !length(pipeline_input_files) || !all(file.exists(pipeline_input_files))) {
  stop("A selected input file is missing.")
}
if (tolower(tools::file_ext(metadata_file)) != "csv") {
  stop("Metadata file must have a .csv extension: ", metadata_file)
}
meta <- read.table(metadata_file, header = TRUE, sep = ",",
                   colClasses = "character", check.names = FALSE,
                   stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM",
                   quote = "\"", comment.char = "")
if (!"SampleID" %in% names(meta)) {
  stop("Column 'SampleID' was not found in metadata file: ", metadata_file)
}
sample_ids <- meta[["SampleID"]]
if (anyNA(sample_ids) || any(!nzchar(sample_ids)) || anyDuplicated(sample_ids)) {
  stop("Column 'SampleID' must contain unique, nonempty values.")
}
if (run_stage_max >= 2L && !group_column %in% names(meta)) stop("Group column missing in metadata.")
fwd_pattern <- "(_R1|_1)(_001)?\\.(fastq|fq)(\\.gz)?$"
rev_pattern <- "(_R2|_2)(_001)?\\.(fastq|fq)(\\.gz)?$"
fwd <- sort(pipeline_input_files[grepl(fwd_pattern, basename(pipeline_input_files))])
rev <- sort(pipeline_input_files[grepl(rev_pattern, basename(pipeline_input_files))])
raw_sample_ids <- sub(fwd_pattern, "", basename(fwd))
raw_reverse_ids <- sub(rev_pattern, "", basename(rev))
if (!length(fwd) || anyDuplicated(raw_sample_ids) ||
    !identical(raw_sample_ids, raw_reverse_ids)) {
  stop("FASTQ pairs must match, with unique sample names; use the same input folder.")
}
sample_ids <- unname(pipeline_sample_id_map[raw_sample_ids])
if (anyNA(sample_ids) || anyDuplicated(sample_ids)) {
  stop("FASTQ sample names do not have unique metadata SampleID mappings.")
}
idx <- match(sample_ids, meta[["SampleID"]])
if (anyNA(idx)) stop("FASTQ SampleIDs missing from the selected metadata file.")
groups <- if (run_stage_max >= 2L) meta[[group_column]][idx] else rep("", length(idx))
metadata_control <- function(column, level) {
  if (is.null(column)) return(NULL)
  if (!length(column) || !nzchar(column)) return(rep(FALSE, length(idx)))
  if (!column %in% names(meta)) stop("Control metadata column missing: ", column)
  if (is.null(level) || !length(level) || !nzchar(level)) {
    stop("A control metadata level was not selected for column: ", column)
  }
  values <- meta[[column]][idx]
  !is.na(values) & values == level
}
is_neg_control <- if (is.null(negative_control_column)) {
  is_control(sample_ids, negative_control_prefixes)
} else {
  metadata_control(negative_control_column, negative_control_level)
}
is_pos_control <- if (is.null(mock_control_column)) {
  is_control(sample_ids, positive_control_prefixes)
} else {
  metadata_control(mock_control_column, mock_control_level)
}
names(is_neg_control) <- sample_ids
names(is_pos_control) <- sample_ids
if (run_stage_max >= 3L && run_decontam && (!any(is_neg_control) || all(is_neg_control))) {
  stop("decontam requires at least one Negative Control and one non-negative sample.")
}
use_for_diversity <- if (exclude_controls) {
  !(is_neg_control | is_pos_control)
} else if (exclude_mock_controls) {
  !is_pos_control
} else {
  rep(TRUE, length(sample_ids))
}
analysis_groups <- groups[use_for_diversity]
if (run_stage_max >= 3L && (anyNA(analysis_groups) || any(!nzchar(trimws(analysis_groups))) ||
    length(unique(analysis_groups)) < 2L ||
    any(table(analysis_groups) < 2L))) {
  stop("Full diversity comparisons need at least two groups and two samples per group.")
}
if (run_stage_max >= 3L && !reference_group %in% analysis_groups) {
  stop("Reference group is absent from samples retained for analysis.")
}
message("Inputs validated: ", length(sample_ids), " paired-end samples.")
if ("--check" %in% args) {
  message("CHECK_OK")
  quit(save = "no", status = 0)
}
setwd(cfg$run_dir)
started <- Sys.time()
for (i in 2:length(step_files)) {
  script_number <- i - 1L
  if (run_stage_max == 1L && script_number > 8L) next
  if (run_stage_max == 2L && script_number %in% 11:14) next
  if (run_stage_max < 3L && script_number == 16L) next
  if (script_number == 11L && !any(beta_distances %in% c("wunifrac", "unifrac"))) next
  if (i == 4L && primer_method == "trim_left") {
    message("\n[STEP ", basename(step_files[i]), "] skipped: primer removal is not needed.")
    next
  }
  message("\n[STEP ", basename(step_files[i]), "] ", format(Sys.time()))
  source(step_files[i], local = .GlobalEnv, encoding = "UTF-8")
  if (i == 11L && run_stage_max >= 3L) {
    if (any(phyloseq::sample_sums(physeq) == 0)) {
      stop("A sample has no reads after taxonomic filtering. ",
           "Review primers, filtering and selected reference database.")
    }
    if (any(beta_distances %in% c("wunifrac", "unifrac")) &&
        phyloseq::ntaxa(physeq) < 3L) {
      stop("UniFrac tree construction requires at least three ASVs after taxonomic filtering. ",
           "Review filtering or select a non-phylogenetic distance metric.")
    }
  }
}
capture.output(sessionInfo(), file = file.path(dir_logs, "sessionInfo.txt"))
message("\nPIPELINE_OK: ", normalizePath(dir_results, winslash = "/"))
message("Elapsed: ", round(as.numeric(difftime(Sys.time(), started, units = "mins")), 2),
        " minutes")
