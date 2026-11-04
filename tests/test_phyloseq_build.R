# Regression for step 10 using real phyloseq objects and synthetic ASV counts.
# Run from the repository root: Rscript tests/test_phyloseq_build.R
suppressPackageStartupMessages({
  library(phyloseq)
  library(tidyverse)
  library(Biostrings)
})

profiles <- jsonlite::fromJSON("pipeline/marker_profiles.json",
                               simplifyVector = FALSE)
launcher <- as.list(parse("run_pipeline.R"))
profile_variables <- c("taxonomy_ranks", "target_domain_rank", "target_domain",
                      "relative_abundance_ranks", "contaminant_rank",
                      "contaminant_taxa", "organelle_taxa")
assignments <- Filter(function(expr) {
  is.call(expr) && identical(expr[[1]], as.name("<-")) && is.symbol(expr[[2]]) &&
    as.character(expr[[2]]) %in% profile_variables
}, launcher)
helpers <- as.list(parse("pipeline/01_helpers.R"))
compute_assignment <- Filter(function(expr) {
  is.call(expr) && identical(expr[[1]], as.name("<-")) &&
    identical(expr[[2]], as.name("compute_and_save"))
}, helpers)[[1]]

run_root <- tempfile("biome-phyloseq-test-")
dir.create(run_root)
tryCatch({
  for (test_case in c("16S", "18S", "18S_without_organelle_field")) {
    marker <- if (test_case == "16S") "16S" else "18S"
    context <- new.env(parent = globalenv())
    context$marker_profile <- profiles[[marker]]
    if (test_case == "18S_without_organelle_field") {
      context$marker_profile$organelle_taxa <- NULL
    }
    for (expr in assignments) eval(expr, context)
    eval(compute_assignment, context)
    context$dir_rds <- file.path(run_root, test_case, "rds")
    context$dir_tables <- file.path(run_root, test_case, "tables")
    dir.create(context$dir_rds, recursive = TRUE)
    dir.create(context$dir_tables, recursive = TRUE)
    samples <- paste0("Sample", 1:4)
    sequences <- c("ACGTACGT", "ACGTACGA", "ACGTACGG", "ACGTACGC", "ACGTACCC")
    counts <- matrix(10L, nrow = 4, ncol = 5,
                     dimnames = list(samples, sequences))
    taxonomy <- matrix("Unclassified", nrow = 5,
                       ncol = length(context$taxonomy_ranks),
                       dimnames = list(sequences, context$taxonomy_ranks))
    taxonomy[, context$target_domain_rank] <- context$target_domain
    taxonomy[3, "Family"] <- " Mitochondria "
    taxonomy[4, "Class"] <- "cHlOrOpLaSt"
    taxonomy[5, context$target_domain_rank] <-
      if (marker == "16S") "Archaea" else "Bacteria"
    context$metadata <- data.frame(SampleID = samples,
                                    Group = c("Control", "Control", "Treatment", "Treatment"))
    context$is_neg_control <- setNames(rep(FALSE, 4), samples)
    context$run_decontam <- FALSE
    context$remove_unassigned <- TRUE
    context$exclude_controls <- FALSE
    context$exclude_mock_controls <- FALSE
    saveRDS(counts, file.path(context$dir_rds, "seqtab_final.rds"))
    saveRDS(taxonomy, file.path(context$dir_rds, "taxonomy.rds"))
    saveRDS(data.frame(SampleID = samples),
            file.path(context$dir_rds, "fastq_info.rds"))

    source("pipeline/10_phyloseq_build.R", local = context)
    expected_taxa <- if (marker == "16S") 2L else 4L
    stopifnot(ntaxa(context$physeq) == expected_taxa)
    stopifnot(nsamples(context$physeq) == 4L)
    stopifnot(file.exists(file.path(context$dir_rds, "physeq.rds")))
    stopifnot(file.exists(file.path(context$dir_tables, "asv_sequences.fasta")))
    cat(test_case, ": phyloseq build and taxonomic filtering passed.\n")
  }
}, finally = {
  unlink(run_root, recursive = TRUE)
})
