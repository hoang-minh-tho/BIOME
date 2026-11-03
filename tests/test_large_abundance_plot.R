# Run from the repository root: Rscript tests/test_large_abundance_plot.R
suppressPackageStartupMessages({ library(phyloseq); library(tidyverse) })
set.seed(1234)
root <- tempfile("biome-large-plot-")
dir.create(root)
tryCatch({
  context <- new.env(parent = globalenv())
  for (folder in c("rds", "tables", "figures")) {
    path <- file.path(root, folder)
    dir.create(path)
    assign(paste0("dir_", folder), path, context)
  }
  context$SEED <- 1234L
  context$figure_format <- "pdf"
  context$group_column <- "Group"
  context$relative_abundance_ranks <- c("Phylum", "Genus")
  context$min_rank1_abundance <- 0
  context$min_rank2_abundance <- 0
  for (expr in as.list(parse("pipeline/01_helpers.R"))) {
    if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
        identical(expr[[2]], as.name("save_plot"))) eval(expr, context)
  }
  # The previous width formula exceeded 50 inches for this input.
  counts <- matrix(rep(c(10L, 20L, 30L, 40L), each = 110L), nrow = 110L,
                    dimnames = list(paste0("S", 1:110), paste0("ASV", 1:4)))
  taxonomy <- cbind(Kingdom = rep("Bacteria", 4),
                    Phylum = paste0("P", 1:4), Genus = paste0("G", 1:4))
  rownames(taxonomy) <- colnames(counts)
  metadata <- data.frame(Group = rep(c("Control", "Treatment"), each = 55),
                          row.names = rownames(counts))
  physeq <- phyloseq(otu_table(counts, taxa_are_rows = FALSE),
                      tax_table(taxonomy), sample_data(metadata))
  saveRDS(physeq, file.path(context$dir_rds, "physeq.rds"))
  for (format in c("png", "pdf")) {
    context$figure_format <- format
    messages <- character()
    withCallingHandlers(sys.source("pipeline/15_relative_abundance.R", context),
      message = function(message) {
        text <- conditionMessage(message)
        if (startsWith(text, "BIOME_PLOT_WARNING:")) {
          messages <<- c(messages, text)
          invokeRestart("muffleMessage")
        }
      })
    status <- read.csv(file.path(context$dir_tables, "relative_abundance_status.csv"))
    stopifnot(all(status$status == "run"), all(status$width_inches <= 24),
               all(status$height_inches <= 20), length(messages) == 2L)
    stopifnot(all(grepl("Collapse", status$plot_warning)))
    stopifnot(context$min_rank1_abundance == 0, context$min_rank2_abundance == 0)
    for (rank in c("phylum", "genus")) {
      data <- read.csv(file.path(context$dir_tables, paste0("relative_abundance_", rank, ".csv")))
      stopifnot(nrow(data) == 440L)
      stopifnot(!any(data[[tools::toTitleCase(rank)]] == "Others"))
      stopifnot(all(abs(tapply(data$Abundance, data$Sample, sum) - 1) < 1e-8))
      plot <- file.path(context$dir_figures, paste0("relative_abundance_", rank, ".", format))
      stopifnot(file.exists(plot), file.info(plot)$size > 0)
    }
  }
}, finally = unlink(root, recursive = TRUE))
cat("Large zero-threshold plots export within size limits; all taxa and CSV values are retained.\n")
