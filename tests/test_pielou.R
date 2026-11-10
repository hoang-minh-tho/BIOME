# Load the shared notation without running helper setup or metadata loading.
for (expr in as.list(parse("pipeline/01_helpers.R"))) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      identical(expr[[2]], as.name("significance_label"))) eval(expr)
}
suppressPackageStartupMessages({library(phyloseq); library(tidyverse)})
root <- tempfile("biome-pielou-")
dir.create(root)
samples <- paste0("Sample", 1:6)
counts <- rbind(c(10, 10), c(19, 1), c(20, 0), c(8, 12), c(15, 5), c(20, 0))
dimnames(counts) <- list(samples, c("ASV1", "ASV2"))
helpers <- as.list(parse("pipeline/01_helpers.R"))
new_context <- function(name, data, groups, measures, statistics) {
  context <- new.env(parent = globalenv())
  for (folder in c("rds", "tables", "figures")) {
    path <- file.path(root, name, folder)
    dir.create(path, recursive = TRUE)
    assign(paste0("dir_", folder), path, context)
  }
  context$figure_format <- "png"
  context$alpha_measures <- measures
  context$run_alpha_statistics <- statistics
  context$group_column <- "Group"
  metadata <- data.frame(Group = groups, row.names = samples)
  physeq <- phyloseq(otu_table(data, taxa_are_rows = FALSE), sample_data(metadata))
  saveRDS(physeq, file.path(context$dir_rds, "physeq_rarefied.rds"))
  for (expr in helpers) {
    if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
        identical(expr[[2]], as.name("save_plot"))) eval(expr, context)
  }
  context
}
tryCatch({
  for (groups in c(1L, 2L, 3L)) {
    for (statistics in c(FALSE, TRUE)) {
      context <- new_context(paste(groups, statistics), counts,
                              rep(LETTERS[seq_len(groups)], each = 6 / groups),
                              "Pielou", statistics)
      sys.source("pipeline/13_alpha_diversity.R", context)
      result <- read.csv(file.path(context$dir_tables, "alpha_diversity.csv"))
      stopifnot("Pielou" %in% names(result), !any(c("Shannon", "Observed") %in% names(result)),
                abs(result$Pielou[1] - 1) < 1e-10,
                abs(result$Pielou[2] - (-0.95 * log(0.95) - 0.05 * log(0.05)) / log(2)) < 1e-10,
                all(is.na(result$Pielou[c(3, 6)])),
                all(result$Pielou[!is.na(result$Pielou)] >= 0),
                all(result$Pielou[!is.na(result$Pielou)] <= 1 + 1e-10),
                file.exists(file.path(context$dir_figures, "alpha_diversity_boxplots.png")))
      if (statistics && groups > 1L) {
        filename <- "alpha_diversity_overall.csv"
        tests <- read.csv(file.path(context$dir_tables, filename))
        stopifnot(tests$Metric == "Pielou", is.finite(tests$p_value))
        if (groups == 3L) {
          tests <- read.csv(file.path(context$dir_tables, "alpha_diversity_pairwise.csv"))
          stopifnot(nrow(tests) == 3L, all(tests$Metric == "Pielou"))
        }
      }
    }
  }
  for (groups in c(2L, 3L)) {
    singleton <- counts
    singleton[, 1] <- 20
    singleton[, 2] <- 0
    context <- new_context(paste0("undefined", groups), singleton,
                            rep(LETTERS[seq_len(groups)], each = 6 / groups),
                            c("Observed", "Shannon", "Pielou"), TRUE)
    sys.source("pipeline/13_alpha_diversity.R", context)
    result <- read.csv(file.path(context$dir_tables, "alpha_diversity.csv"))
    stopifnot(all(is.na(result$Pielou)), all(result$Observed == 1), all(result$Shannon == 0))
    filename <- "alpha_diversity_overall.csv"
    tests <- read.csv(file.path(context$dir_tables, filename))
    stopifnot(is.na(tests$p_value[tests$Metric == "Pielou"]),
              tests$Status[tests$Metric == "Pielou"] == "not testable")
  }
  context <- new_context("undefined_only", singleton, rep(c("A", "B"), each = 3), "Pielou", TRUE)
  sys.source("pipeline/13_alpha_diversity.R", context)
  result <- read.csv(file.path(context$dir_tables, "alpha_diversity.csv"))
  stopifnot(all(is.na(result$Pielou)),
            file.exists(file.path(context$dir_figures, "alpha_diversity_boxplots.png")))
}, finally = unlink(root, recursive = TRUE))
cat("Pielou formula, independent selection, undefined values, plots, and group comparisons passed.\n")
