# Load the shared notation without running helper setup or metadata loading.
for (expr in as.list(parse("pipeline/01_helpers.R"))) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      identical(expr[[2]], as.name("significance_label"))) eval(expr)
}
# Regression tests for ASV bounds, sparse groups, missing ranks, and tree guards.
# Run from the repository root: Rscript tests/test_pipeline_edge_cases.R
suppressPackageStartupMessages({
  library(phyloseq)
  library(tidyverse)
  library(vegan)
})
set.seed(1234)
launcher <- as.list(parse("run_pipeline.R"))
assignment_to <- function(expr, name) {
  is.call(expr) && identical(expr[[1]], as.name("<-")) &&
    identical(expr[[2]], as.name(name))
}
start <- which(vapply(launcher, assignment_to, logical(1),
                      name = "asv_length_filter_enabled"))
end <- which(vapply(launcher, assignment_to, logical(1), name = "trunc_qmin")) - 1L
profiles <- jsonlite::fromJSON("pipeline/marker_profiles.json", simplifyVector = FALSE)
for (marker in c("16S", "18S")) {
  for (enabled in c(TRUE, FALSE)) {
    context <- new.env(parent = globalenv())
    context$marker_profile <- profiles[[marker]]
    sys.source("pipeline/00_parameters.R", context)
    context$cfg <- list(marker = marker, asv_length_filter_enabled = enabled,
                        asv_min_len = 300L, asv_max_len = 500L)
    for (expr in launcher[start:end]) eval(expr, context)
    if (enabled) {
      stopifnot(context$asv_min_len == 300L, context$asv_max_len == 500L)
      context$cfg$asv_max_len <- 280L
      error <- tryCatch({
        for (expr in launcher[start:end]) eval(expr, context)
        NULL
      }, error = identity)
      stopifnot(inherits(error, "error"))
    } else {
      stopifnot(is.null(context$asv_min_len), is.null(context$asv_max_len))
    }
  }
}
cat("16S/18S ASV length bounds passed.\n")

make_physeq <- function(counts, genus = c("G1", "G2", "G3", "G4")) {
  dimnames(counts) <- list(paste0("Sample", seq_len(nrow(counts))),
                           paste0("ASV", seq_len(ncol(counts))))
  taxonomy <- cbind(Kingdom = rep("Bacteria", ncol(counts)),
                    Phylum = rep("Firmicutes", ncol(counts)), Genus = genus)
  rownames(taxonomy) <- colnames(counts)
  metadata <- data.frame(Group = rep(c("Control", "Treatment"), each = 3),
                          row.names = rownames(counts))
  phyloseq(otu_table(counts, taxa_are_rows = FALSE), tax_table(taxonomy),
            sample_data(metadata))
}
counts <- rbind(c(40, 30, 20, 10), c(30, 40, 10, 20), c(20, 10, 40, 30),
                c(10, 30, 20, 40), c(20, 30, 40, 10), c(30, 10, 20, 40))
root <- tempfile("biome-edge-tests-")
dir.create(root)
helpers <- as.list(parse("pipeline/01_helpers.R"))
new_context <- function(name, physeq) {
  context <- new.env(parent = globalenv())
  for (folder in c("rds", "tables", "figures")) {
    path <- file.path(root, name, folder)
    dir.create(path, recursive = TRUE)
    assign(paste0("dir_", folder), path, context)
  }
  for (expr in helpers) {
    if (any(vapply(c("compute_and_save", "save_plot", "save_base_plot"),
                    function(name) assignment_to(expr, name), logical(1)))) {
      eval(expr, context)
    }
  }
  context$SEED <- 1234L
  context$figure_format <- "png"
  context$group_column <- "Group"
  context$alpha_measures <- c("Observed", "Chao1", "Shannon", "Simpson")
  context$run_alpha_statistics <- TRUE
  context$rarefaction_mode <- "custom"
  context$custom_rarefaction_depth <- 50L
  context$ordination_method <- "PCoA"
  context$beta_distances <- "bray"
  context$run_permanova <- TRUE
  context$run_betadisper <- TRUE
  context$permanova_permutations <- 9L
  context$relative_abundance_ranks <- c("Phylum", "Genus")
  context$min_rank1_abundance <- 0.01
  context$min_rank2_abundance <- 0.02
  saveRDS(physeq, file.path(context$dir_rds, "physeq.rds"))
  context
}
run_step <- function(path, context) {
  warnings <- character()
  withCallingHandlers(sys.source(path, context), warning = function(warning) {
    warnings <<- c(warnings, conditionMessage(warning))
    invokeRestart("muffleWarning")
  })
  warnings
}
tryCatch({
  for (groups in c(1L, 2L, 3L)) {
    for (enabled in c(FALSE, TRUE)) {
      physeq <- make_physeq(counts)
      sample_data(physeq)$Group <- rep(LETTERS[seq_len(groups)], length.out = 6)
      context <- new_context(paste0("alpha_", groups, "_", enabled), physeq)
      context$alpha_measures <- c("Shannon", "Simpson")
      context$run_alpha_statistics <- enabled
      saveRDS(physeq, file.path(context$dir_rds, "physeq_rarefied.rds"))
      run_step("pipeline/13_alpha_diversity.R", context)
      metrics <- read.csv(file.path(context$dir_tables, "alpha_diversity.csv"))
      stopifnot(all(c("Shannon", "Simpson") %in% names(metrics)),
                !any(c("Observed", "Chao1") %in% names(metrics)),
                file.exists(file.path(context$dir_figures, "alpha_diversity_boxplots.png")))
      results <- list.files(context$dir_tables, pattern = "alpha_diversity_(overall|assumptions|pairwise)[.]csv")
      stopifnot(length(results) == if (enabled) 3L else 0L)
      if (length(results)) {
        for (result in results) {
          table <- read.csv(file.path(context$dir_tables, result))
          stopifnot("p_value" %in% names(table))
        }
        stopifnot(all(grepl("p [=<]|not tested", unlist(context$p_labels))),
                  if (groups %in% 2:3) grepl("BH-adjusted", context$alpha_plot$labels$caption) else TRUE)
      } else stopifnot(length(context$p_labels) == 0L)
    }
  }
  cat("Alpha metric selection, statistics toggle, and automatic group tests passed.\n")

  for (case in c("lost_group", "two_samples", "singleton_group", "balanced", "low_asv")) {
    input <- counts
    if (case == "lost_group") input[1:3, ] <- 5L
    if (case == "two_samples") input[1:4, ] <- 5L
    if (case == "singleton_group") input[2:3, ] <- 5L
    physeq <- if (case == "low_asv") {
      make_physeq(input[, 1:2] * 3L, c("G1", "G2"))
    } else {
      make_physeq(input)
    }
    context <- new_context(case, physeq)
    run_step("pipeline/12_rarefaction.R", context)
    run_step("pipeline/13_alpha_diversity.R", context)
    for (dispersion in c(TRUE, FALSE)) {
      context$run_betadisper <- dispersion
      warnings <- run_step("pipeline/14_beta_diversity.R", context)
      status <- read.csv(file.path(context$dir_tables, "beta_statistics_status.csv"))
      stopifnot(file.exists(file.path(context$dir_rds, "dist_bray.rds")))
      expected <- if (case %in% c("balanced", "low_asv")) "run" else "skipped"
      stopifnot(status$status == expected)
      if (expected == "skipped") {
        stopifnot(nzchar(status$reason))
        stopifnot(any(grepl("Skipping PERMANOVA", warnings)))
        stopifnot(!file.exists(file.path(context$dir_tables, "permanova_results.csv")))
      } else {
        stopifnot(file.exists(file.path(context$dir_tables, "permanova_results.csv")))
      }
    }
  }
  cat("Lost groups, sparse replicates, and valid beta tests passed.\n")

  for (case in c("all_missing", "blank", "partial", "assigned")) {
    genus <- switch(case, all_missing = rep(NA_character_, 4),
                     blank = c("", " ", "\t", ""),
                     partial = c("G1", NA, NA, NA), assigned = paste0("G", 1:4))
    context <- new_context(case, make_physeq(counts, genus))
    run_step("pipeline/15_relative_abundance.R", context)
    data <- read.csv(file.path(context$dir_tables, "relative_abundance_genus.csv"))
    status <- read.csv(file.path(context$dir_tables, "relative_abundance_status.csv"))
    expected <- if (case %in% c("all_missing", "blank")) "skipped" else "run"
    stopifnot(status$status[status$rank == "Genus"] == expected)
    stopifnot(file.exists(file.path(context$dir_figures, "relative_abundance_phylum.png")))
    if (expected == "skipped") {
      stopifnot(nrow(data) == 0L)
      stopifnot(!file.exists(file.path(context$dir_figures, "relative_abundance_genus.png")))
      context$reference_group <- "Control"
      context$differential_method <- "ANCOM-BC"
      run_step("pipeline/16_ancombc.R", context)
      da_status <- read.csv(file.path(context$dir_tables, "differential_abundance_status.csv"))
      stopifnot(da_status$status == "skipped")
    } else {
      stopifnot(nrow(data) > 0L, all(is.finite(data$Abundance)))
      stopifnot(all(abs(tapply(data$Abundance, data$Sample, sum) - 1) < 1e-8))
    }
  }
  cat("Missing, blank, partial, and complete rank assignments passed.\n")

  loop <- Filter(function(expr) is.call(expr) && identical(expr[[1]], as.name("for")),
                   launcher)[[1]]
  guard <- tail(Filter(function(expr) is.call(expr) && identical(expr[[1]], as.name("if")),
                         as.list(loop[[4]])[-1]), 1)[[1]]
  context <- new.env(parent = globalenv())
  context$i <- 11L
  context$run_stage_max <- 3L
  sparse <- make_physeq(counts[, 1:2], c("G1", "G2"))
  for (metric in c("bray", "jaccard", "unifrac", "wunifrac")) {
    context$physeq <- sparse
    context$beta_distances <- metric
    error <- tryCatch({ eval(guard, context); NULL }, error = identity)
    stopifnot(inherits(error, "error") == (metric %in% c("unifrac", "wunifrac")))
  }
  context$physeq <- make_physeq(counts)
  eval(guard, context)
  context$beta_distances <- "bray"
  zero <- counts
  zero[1, ] <- 0L
  context$physeq <- make_physeq(zero)
  error <- tryCatch({ eval(guard, context); NULL }, error = identity)
  stopifnot(inherits(error, "error"))
  cat("Bray/Jaccard low-ASV acceptance and UniFrac/empty-sample guards passed.\n")
}, finally = unlink(root, recursive = TRUE))
