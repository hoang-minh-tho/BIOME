# Read counts must be captured before taxonomic filtering and control exclusion.
suppressPackageStartupMessages({ library(phyloseq); library(tidyverse); library(Biostrings) })
root <- tempfile("biome-decontam-tracking-")
dir.create(root)
helpers <- as.list(parse("pipeline/01_helpers.R"))
compute <- Filter(function(expr) {
  is.call(expr) && identical(expr[[1]], as.name("<-")) &&
    identical(expr[[2]], as.name("compute_and_save"))
}, helpers)[[1]]
tryCatch({
  for (marker in c("16S", "18S")) {
    context <- new.env(parent = globalenv())
    eval(compute, context)
    context$dir_rds <- file.path(root, marker, "rds")
    context$dir_tables <- file.path(root, marker, "tables")
    dir.create(context$dir_rds, recursive = TRUE)
    dir.create(context$dir_tables, recursive = TRUE)
    samples <- paste0("Sample", 1:8)
    sequences <- c("ACGTACGT", "ACGTACGA", "ACGTACGG")
    counts <- cbind(rep(100L, 8), c(rep(50L, 4), rep(0L, 4)), rep(25L, 8))
    dimnames(counts) <- list(samples, sequences)
    rank1 <- if (marker == "16S") "Kingdom" else "Domain"
    target <- if (marker == "16S") "Bacteria" else "Eukaryota"
    taxonomy <- cbind(rep(target, 3), c("Genuine", "Contaminant", "Nontarget"))
    colnames(taxonomy) <- c(rank1, "Genus")
    taxonomy[3, rank1] <- "ExcludedDomain"
    rownames(taxonomy) <- sequences
    context$metadata <- data.frame(SampleID = samples,
                                   Group = rep(c("Negative", "Treatment"), each = 4))
    context$is_neg_control <- setNames(rep(c(TRUE, FALSE), each = 4), samples)
    context$is_pos_control <- setNames(rep(FALSE, 8), samples)
    context$run_decontam <- TRUE
    context$remove_unassigned <- TRUE
    context$exclude_controls <- TRUE
    context$exclude_mock_controls <- FALSE
    context$target_domain_rank <- rank1
    context$target_domain <- target
    context$contaminant_taxa <- character()
    context$contaminant_rank <- "Genus"
    context$organelle_taxa <- character()
    saveRDS(counts, file.path(context$dir_rds, "seqtab_final.rds"))
    saveRDS(taxonomy, file.path(context$dir_rds, "taxonomy.rds"))
    saveRDS(data.frame(SampleID = samples), file.path(context$dir_rds, "fastq_info.rds"))
    tracking <- data.frame(input = rowSums(counts) + 20L, length_filtered = rowSums(counts),
                           row.names = samples)
    tracking <- tracking[rev(samples), , drop = FALSE]
    saveRDS(tracking, file.path(context$dir_rds, "read_tracking.rds"))
    sys.source("pipeline/10_phyloseq_build.R", context)
    updated <- readRDS(file.path(context$dir_rds, "read_tracking.rds"))
    decontaminated <- readRDS(file.path(context$dir_rds, "seqtab_decontam.rds"))
    stopifnot(identical(rownames(updated), rev(samples)),
              identical(updated$input, tracking$input),
              identical(updated$length_filtered, tracking$length_filtered),
              all(updated$decontaminated == rowSums(decontaminated)[rownames(updated)]),
              all(updated$decontaminated == 125),
              all(updated$decontam_removed == ifelse(rownames(updated) %in% samples[1:4], 50, 0)),
              all(updated$mitochondria_removed == 0),
              all(updated$chloroplast_removed == 0),
              all(updated$other_taxa_removed == 25),
              all(updated$taxonomy_filtered == 100),
              all(updated$decontaminated == updated$mitochondria_removed +
                  updated$chloroplast_removed + updated$other_taxa_removed +
                  updated$taxonomy_filtered),
              nsamples(context$physeq) == 4L,
              all(sample_sums(context$physeq) == 100))
    csv <- read.csv(file.path(context$dir_tables, "read_tracking.csv"), row.names = 1)
    stopifnot(identical(names(csv), names(updated)),
              all(csv$decontaminated == updated$decontaminated))
    # No contaminants removed: still report the post-decontam counts explicitly.
    result <- context$update_decontam_read_tracking(counts, context$dir_rds, context$dir_tables)
    stopifnot(all(result$decontam_removed == 0))
    error <- tryCatch(context$update_decontam_read_tracking(
      counts[-1, , drop = FALSE], context$dir_rds, context$dir_tables), error = identity)
    stopifnot(inherits(error, "error"), grepl("Sample IDs", conditionMessage(error)))
  }
}, finally = unlink(root, recursive = TRUE))
cat("16S/18S post-decontam read tracking, sample alignment, CSV/RDS, and filtering boundaries passed.\n")
