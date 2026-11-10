# Run from the repository root: Rscript tests/test_abundance_heatmap.R
suppressPackageStartupMessages({ library(phyloseq); library(tidyverse) })
stopifnot(requireNamespace("pheatmap", quietly = TRUE))
set.seed(1234)
stage <- as.list(parse("pipeline/15_relative_abundance.R"))
helpers <- as.list(parse("pipeline/01_helpers.R"))
definition <- function(expr, name) {
  is.call(expr) && identical(expr[[1]], as.name("<-")) && identical(expr[[2]], as.name(name))
}
for (expr in stage) {
  if (definition(expr, "prepare_abundance_heatmap") || definition(expr, "save_abundance_heatmap")) {
    eval(expr)
  }
}
root <- tempfile("biome-heatmap-")
dir.create(root)
make_input <- function(rank1) {
  counts <- matrix(sample.int(100, 44 * 6, replace = TRUE), nrow = 44,
                   dimnames = list(paste0("ASV", 1:44), paste0("Sample", 1:6)))
  counts[44, ] <- 0
  counts[1, 1:3] <- 2000
  counts[2, 4:6] <- 1800
  taxonomy <- cbind(Domain = rep("Target", 44), Parent = paste0("P", 1:44),
                    Genus = paste0("G", 1:44))
  colnames(taxonomy)[2] <- rank1
  taxonomy[1:2, "Genus"] <- "SharedName" # Distinct parent lineages.
  taxonomy[3:5, "Genus"] <- c(NA, "", "  ")
  # Two ASVs with identical taxonomy must aggregate to a single row.
  taxonomy[7, ] <- taxonomy[6, ]
  rownames(taxonomy) <- rownames(counts)
  metadata <- data.frame(arm = rep(c("Control", "Treatment"), each = 3),
                          row.names = colnames(counts))
  colnames(metadata) <- "Study arm"
  phyloseq(otu_table(counts, taxa_are_rows = TRUE), tax_table(taxonomy), sample_data(metadata))
}
new_context <- function(name, physeq) {
  context <- new.env(parent = globalenv())
  for (folder in c("rds", "tables", "figures")) {
    path <- file.path(root, name, folder)
    dir.create(path, recursive = TRUE)
    assign(paste0("dir_", folder), path, context)
  }
  context$group_column <- "Study arm"
  context$relative_abundance_ranks <- c(rank_names(physeq)[2], "Genus")
  context$SEED <- 1234L
  context$min_rank1_abundance <- 0.01
  context$min_rank2_abundance <- 0.02
  context$abundance_plot_type <- "heatmap"
  context$heatmap_rank <- "Genus"
  context$heatmap_top_n <- 30L
  context$figure_format <- "png"
  for (expr in helpers) {
    if (definition(expr, "save_base_plot") || definition(expr, "save_plot")) eval(expr, context)
  }
  saveRDS(physeq, file.path(context$dir_rds, "physeq.rds"))
  # There is deliberately no rarefied RDS in this run.
  context
}
tryCatch({
  for (rank1 in c("Phylum", "Division")) {
    physeq <- make_input(rank1)
    full <- prepare_abundance_heatmap(physeq, "Genus", 100L, "Study arm")
    stopifnot(nrow(full$abundance) == 42L, ncol(full$abundance) == 6L,
              !anyDuplicated(rownames(full$abundance)),
              sum(grepl("^SharedName", rownames(full$abundance))) == 2L,
              sum(grepl("^Unclassified", rownames(full$abundance))) == 3L,
              !any(c("Other", "Others", "", NA) %in% rownames(full$abundance)),
              !"ASV44" %in% full$taxon_ids,
              all(abs(colSums(full$abundance) - 1) < 1e-10),
              all(is.finite(full$z_scores)),
              all(abs(rowMeans(full$z_scores)) < 1e-10),
              identical(rownames(full$annotation), colnames(full$abundance)),
              identical(as.character(full$annotation[["Study arm"]]),
                        as.character(sample_data(physeq)[["Study arm"]])))
    for (top_n in c(10L, 20L, 25L, 30L, 40L, 50L, 100L)) {
      data <- prepare_abundance_heatmap(physeq, "Genus", top_n, "Study arm")
      stopifnot(nrow(data$abundance) == min(top_n, 42L),
                identical(data$taxon_ids, head(full$taxon_ids, top_n)),
                isTRUE(all.equal(data$abundance, full$abundance[seq_len(nrow(data$abundance)), , drop = FALSE])))
    }
    # No normalization to the selected subset; same global taxa in every column.
    data <- prepare_abundance_heatmap(physeq, "Genus", 10L, "Study arm")
    stopifnot(all(colSums(data$abundance) < 1))
    # Reordering metadata/samples must not change group colors or Top N selection.
    reordered <- prepare_abundance_heatmap(prune_samples(rev(sample_names(physeq)), physeq),
                                           "Genus", 10L, "Study arm")
    stopifnot(identical(data$taxon_ids, reordered$taxon_ids),
              identical(data$annotation_colors, reordered$annotation_colors))
    for (rank in c("Genus", rank1)) {
      for (format in c("png", "pdf")) {
        context <- new_context(paste(rank1, rank, format, sep = "_"), physeq)
        context$heatmap_rank <- rank
        context$figure_format <- format
        sys.source("pipeline/15_relative_abundance.R", context)
        expected <- paste0("heatmap_top30_", tolower(rank), ".", format)
        stopifnot(identical(list.files(context$dir_figures), expected),
                  file.info(file.path(context$dir_figures, expected))$size > 1000,
                  !file.exists(file.path(context$dir_tables, "relative_abundance_status.csv")),
                  read.csv(file.path(context$dir_tables, "heatmap_status.csv"))$status == "run")
        if (format == "png") {
          bytes <- readBin(file.path(context$dir_figures, expected), "raw", n = 1000)
          position <- grepRaw("pHYs", bytes, fixed = TRUE)[1]
          # PNG resolution: 300 DPI = 11811 pixels per metre.
          value <- sum(as.integer(bytes[position + 4:7]) * 256^(3:0))
          stopifnot(value == 11811)
        } else {
          header <- readChar(file.path(context$dir_figures, expected), 4L)
          stopifnot(header == "%PDF")
        }
      }
    }
    # Barplots remain the default for direct callers without new settings.
    context <- new_context(paste0(rank1, "_legacy"), physeq)
    rm("abundance_plot_type", envir = context)
    context$group_column <- "Group"
    sample_data(physeq)$Group <- sample_data(physeq)[["Study arm"]]
    saveRDS(physeq, file.path(context$dir_rds, "physeq.rds"))
    sys.source("pipeline/15_relative_abundance.R", context)
    stopifnot(length(list.files(context$dir_figures, pattern = "^relative_abundance_")) == 2L,
              !length(list.files(context$dir_figures, pattern = "^heatmap_")))
    context <- new_context(paste0(rank1, "_both"), physeq)
    context$group_column <- "Group"
    context$abundance_plot_type <- "both"
    sys.source("pipeline/15_relative_abundance.R", context)
    stopifnot(length(list.files(context$dir_figures, pattern = "^relative_abundance_")) == 2L,
              length(list.files(context$dir_figures, pattern = "^heatmap_")) == 1L)
  }
  # Constant rows and a single sample produce finite zeros, without clustering errors.
  single <- prune_taxa("ASV1", prune_samples("Sample1", make_input("Phylum")))
  names(sample_data(single)) <- "Study arm"
  data <- prepare_abundance_heatmap(single, "Genus", 30L, "Study arm")
  stopifnot(all(data$z_scores == 0))
  context <- new_context("single", single)
  sys.source("pipeline/15_relative_abundance.R", context)
  stopifnot(read.csv(file.path(context$dir_tables, "heatmap_status.csv"))$status == "run")
  constant <- prune_taxa("ASV1", make_input("Phylum"))
  names(sample_data(constant)) <- "Study arm"
  data <- prepare_abundance_heatmap(constant, "Genus", 30L, "Study arm")
  stopifnot(all(data$z_scores == 0))
  # A real device failure must be logged and return control to the caller.
  context <- new_context("failure", make_input("Phylum"))
  context$save_base_plot <- function(...) stop("Injected heatmap device error")
  messages <- character()
  withCallingHandlers({
    sys.source("pipeline/15_relative_abundance.R", context)
    context$next_step_reached <- TRUE
  }, message = function(message) {
    messages <<- c(messages, conditionMessage(message))
    invokeRestart("muffleMessage")
  })
  status <- read.csv(file.path(context$dir_tables, "heatmap_status.csv"))
  stopifnot(context$next_step_reached, status$status == "failed",
            grepl("Injected heatmap device error", status$reason),
            any(grepl("ERROR: Abundance heatmap failed", messages)),
            !length(list.files(context$dir_figures)))
  physeq <- make_input("Phylum")
  sample_data(physeq)$Group <- sample_data(physeq)[["Study arm"]]
  context <- new_context("both_failure", physeq)
  context$group_column <- "Group"
  context$abundance_plot_type <- "both"
  context$save_base_plot <- function(...) stop("Injected heatmap device error")
  sys.source("pipeline/15_relative_abundance.R", context)
  stopifnot(read.csv(file.path(context$dir_tables, "heatmap_status.csv"))$status == "failed",
            length(list.files(context$dir_figures, pattern = "^relative_abundance_")) == 2L)
}, finally = unlink(root, recursive = TRUE))
cat("Heatmaps: global Top N, lineage labels, scaling, metadata, 16S/18S, PNG/PDF, exclusivity, and error recovery passed.\n")
