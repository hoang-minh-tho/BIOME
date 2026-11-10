# Load the shared notation without running helper setup or metadata loading.
for (expr in as.list(parse("pipeline/01_helpers.R"))) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      identical(expr[[2]], as.name("significance_label"))) eval(expr)
}
# Run from the repository root: Rscript tests/test_alpha_plot_labels.R
suppressPackageStartupMessages({ library(phyloseq); library(tidyverse) })
set.seed(1234)
root <- tempfile("biome-alpha-labels-")
dir.create(root)
tryCatch({
  context <- new.env(parent = globalenv())
  for (folder in c("rds", "tables", "figures")) {
    path <- file.path(root, folder)
    dir.create(path)
    assign(paste0("dir_", folder), path, context)
  }
  for (expr in as.list(parse("pipeline/01_helpers.R"))) {
    if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
        identical(expr[[2]], as.name("save_plot"))) eval(expr, context)
  }
  context$group_column <- "Group"
  context$alpha_measures <- c("Observed", "Simpson", "Pielou")
  context$run_alpha_statistics <- TRUE
  counts <- matrix(sample(1:100, 120, replace = TRUE), nrow = 12,
                   dimnames = list(paste0("S", 1:12), paste0("ASV", 1:10)))
  for (groups in c(2L, 3L, 4L, 7L)) {
    group_names <- c("G_formosa", "G_millepora", "H_formosa", "H_fungia",
                     "H_millepora", "Seawater", "Snail_gut")[seq_len(groups)]
    group_counts <- counts[rep(seq_len(nrow(counts)), length.out = groups * 3L), ]
    rownames(group_counts) <- paste0("S", seq_len(nrow(group_counts)))
    metadata <- data.frame(Group = rep(group_names, each = 3L),
                           row.names = rownames(group_counts))
    physeq <- phyloseq(otu_table(group_counts, taxa_are_rows = FALSE), sample_data(metadata))
    saveRDS(physeq, file.path(context$dir_rds, "physeq_rarefied.rds"))
    for (format in c("png", "pdf")) {
      context$figure_format <- format
      suppressWarnings(sys.source("pipeline/13_alpha_diversity.R", context))
      stopifnot(identical(context$format_alpha_p(c(1e-9, 0.0234, 0.0041, 1, NA)),
                          c("p = 1.000e-09", "p = 0.0234", "p = 0.0041", "p = 1.0000", "not tested")))
      labels <- unname(context$p_labels)
      tested_labels <- labels[!grepl("not tested", labels)]
      stopifnot(length(labels) == 3L,
                length(tested_labels) > 0L,
                all(lengths(strsplit(tested_labels, "\n")) == 2L),
                !any(grepl("Shapiro|Levene|normality|variance|\\*|\\bns\\b", labels)),
                all(grepl("ANOVA|Kruskal-Wallis", tested_labels)),
                if (groups >= 4L) nrow(context$alpha_annotations) == 0L else TRUE,
                if (groups >= 4L) is.null(context$alpha_plot$labels$caption) else TRUE,
                context$alpha_ncol == if (groups > 4L) 1L else 3L,
                context$alpha_nrow == if (groups > 4L) 3L else 1L,
                identical(context$alpha_facet_labels("Observed", NA_real_, "Kruskal-Wallis"),
                          "Observed\nnot tested"))
      output <- file.path(context$dir_figures, paste0("alpha_diversity_boxplots.", format))
      stopifnot(file.exists(output), file.info(output)$size > 0)
      # Check header text fits the available per-facet width on the export device.
      grDevices::pdf(file.path(root, "text-width.pdf"),
                    width = context$alpha_plot_width, height = context$alpha_plot_height)
      widths <- strwidth(unlist(strsplit(labels, "\n")), units = "inches",
                         cex = 11 / 12, font = 2)
      grDevices::dev.off()
      stopifnot(all(widths < 2.8))
      preview_dir <- Sys.getenv("BIOME_ALPHA_PREVIEW_DIR")
      if (groups == 7L && nzchar(preview_dir)) {
        dir.create(preview_dir, recursive = TRUE, showWarnings = FALSE)
        file.copy(output, file.path(preview_dir, basename(output)), overwrite = TRUE)
      }
    }
    context$run_alpha_statistics <- FALSE
    suppressWarnings(sys.source("pipeline/13_alpha_diversity.R", context))
    stopifnot(is.null(context$alpha_test_caption), length(context$p_labels) == 0L,
              is.null(context$alpha_plot$labels$caption))
    context$run_alpha_statistics <- TRUE
  }
}, finally = unlink(root, recursive = TRUE))
cat("Three-metric alpha plots export in PNG/PDF with automatic tests, clean headers, and group-aware annotations.\n")
