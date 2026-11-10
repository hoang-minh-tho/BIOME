# Device-specific regression checks for stacked legends and long group labels.
suppressPackageStartupMessages({ library(phyloseq); library(tidyverse) })
for (expr in as.list(parse("pipeline/15_relative_abundance.R"))) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      as.character(expr[[2]]) %in% c("prepare_abundance_heatmap", "stack_heatmap_legends",
                                     "save_abundance_heatmap")) eval(expr)
}
for (expr in as.list(parse("pipeline/01_helpers.R"))) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      identical(expr[[2]], as.name("save_base_plot"))) eval(expr)
}
root <- tempfile("biome-legend-tests-")
dir.create(root)
dir_figures <- root
set.seed(1234)
tryCatch({
  for (format in c("png", "pdf")) {
    figure_format <- format
    for (labels in list(c("Coastal_Bay", "Open_Ocean"),
                        c("Coastal_Bay_Outer_Estuary_Seasonal_Monitoring_Station_WWWW",
                          "Open_Ocean_Offshore_Reference_Station"),
                        paste0("Coastal_Bay_Station_", seq_len(20)))) {
      groups <- rep(labels, each = 2)
      abundance <- matrix(runif(12 * length(groups)), nrow = 12,
                          dimnames = list(paste0("Taxon", 1:12), paste0("Sample", seq_along(groups))))
      annotation <- data.frame(Group = factor(groups), row.names = colnames(abundance))
      colors <- setNames(grDevices::hcl.colors(length(labels), "Dark 3"), labels)
      data <- list(z_scores = t(scale(t(abundance))), annotation = annotation,
                   annotation_colors = list(Group = colors))
      # Intercept drawing while the actual PNG/PDF device is active. Measure
      # the resulting viewport bounds as well as rendering the exported figure.
      checked <- FALSE
      trace("stack_heatmap_legends", exit = quote({
        result <- returnValue()
        index <- match("stacked_legends", result$layout$name)
        stopifnot(!is.na(index),
                  !any(c("legend", "annotation_legend") %in% result$layout$name))
        stacked <- result$grobs[[index]]
        stopifnot(stacked$layout$t[1] < stacked$layout$t[2])
        available_width <- grid::convertWidth(sum(stacked$widths), "inches", valueOnly = TRUE)
        for (legend in stacked$grobs) {
          for (child in as.list(legend$children)) {
            if (inherits(child, "text")) {
              right <- grid::convertWidth(child$x + grid::grobWidth(child), "inches", valueOnly = TRUE)
              stopifnot(right <= available_width - 0.2)
            }
          }
        }
        available_height <- grid::convertHeight(
          sum(result$heights[result$layout$t[index]:result$layout$b[index]]), "inches", valueOnly = TRUE
        )
        legend_height <- grid::convertHeight(sum(stacked$heights), "inches", valueOnly = TRUE)
        stopifnot(legend_height < available_height)
        matrix_column <- result$layout$l[match("matrix", result$layout$name)]
        stopifnot(grid::convertWidth(result$widths[matrix_column], "inches", valueOnly = TRUE) > 3)
        checked <<- TRUE
      }), print = FALSE)
      save_abundance_heatmap(data, "Genus", 30L)
      untrace("stack_heatmap_legends")
      stopifnot(checked)
      filename <- file.path(root, paste0("heatmap_top30_genus.", format))
      stopifnot(file.exists(filename), file.info(filename)$size > 1000)
    }
  }
}, finally = {
  if (inherits(stack_heatmap_legends, "functionWithTrace")) untrace("stack_heatmap_legends")
  unlink(root, recursive = TRUE)
})
cat("Stacked legends fit PNG/PDF bounds for Coastal_Bay, long labels, and 20 groups.\n")
