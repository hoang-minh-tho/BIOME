# Verify the shared notation on the production differential-abundance plot.
suppressPackageStartupMessages(library(tidyverse))
for (expr in as.list(parse("pipeline/01_helpers.R"))) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      identical(expr[[2]], as.name("significance_label"))) eval(expr)
}
find_contrast_loop <- function(expr) {
  if (missing(expr)) return(NULL)
  if (!is.call(expr)) return(NULL)
  if (identical(expr[[1]], as.name("for")) &&
      identical(expr[[2]], as.name("contrast_index"))) return(expr)
  for (child in as.list(expr)[-1L]) {
    found <- find_contrast_loop(child)
    if (!is.null(found)) return(found)
  }
  NULL
}
loop <- NULL
for (expr in as.list(parse("pipeline/16_ancombc.R"))) {
  loop <- find_contrast_loop(expr)
  if (!is.null(loop)) break
}
stopifnot(!is.null(loop))
results <- data.frame(taxon = paste0("ASV", 1:7), lfc_GroupB = seq(-3, 3),
                      q_GroupB = c(0.0001, 0.001, 0.01, 0.05, 0.05001, 1, NA))
lfc_cols <- "lfc_GroupB"
reference_group <- "A"
save_plot <- function(plot, filename) {
  stopifnot(identical(as.character(plot$data$Significance),
                      c("****", "***", "**", "*", "ns", "ns", "not tested")),
            plot$labels$colour == "Significance (adjusted p)")
  suppressWarnings(ggplot_build(plot))
}
eval(loop)
for (path in list.files("pipeline", pattern = "[.]R$", full.names = TRUE)) parse(path)
cat("Differential plot significance boundaries and pipeline R syntax passed.\n")
