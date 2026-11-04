# Run from the repository root: Rscript tests/test_abundance_legend.R
suppressPackageStartupMessages(library(tidyverse))
expressions <- as.list(parse("pipeline/15_relative_abundance.R"))
definition <- Filter(function(expr) {
  is.call(expr) && identical(expr[[1]], as.name("<-")) &&
    identical(expr[[2]], as.name("abundance_barplot"))
}, expressions)[[1]]
context <- new.env(parent = globalenv())
context$group_column <- "Group"
eval(definition, context)
for (rank in c("Phylum", "Genus")) {
  data <- data.frame(Sample = c("S1", "S1"), Group = c("Control", "Control"),
                     Abundance = c(0.9, 0.1))
  data[[rank]] <- c("Assigned", "Others")
  palette <- c(Assigned = "#3377EE", Others = "#D3D3D3")
  for (percentage in c(1, 2, 2.75)) {
    plot <- context$abundance_barplot(data, rank, palette, percentage / 100)
    labels <- plot$scales$get_scales("fill")$get_labels(names(palette))
    stopifnot(identical(unname(labels), c("Assigned", paste0("Other (<", percentage, "%)"))))
    stopifnot(identical(data[[rank]], c("Assigned", "Others")))
  }
}
cat("Per-rank legend percentages and stable CSV category names verified.\n")
