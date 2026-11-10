# 1. Helper functions

suppressPackageStartupMessages(library(tidyverse))

# Shared by statistical plots and result tables; thresholds are inclusive.
significance_label <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "not tested",
    p <= 0.0001 ~ "****",
    p <= 0.001 ~ "***",
    p <= 0.01 ~ "**",
    p <= 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

# Each analysis owns these four output folders.
dir_results <- "."
dir_figures <- file.path(dir_results, "figures")
dir_tables <- file.path(dir_results, "tables")
dir_logs <- file.path(dir_results, "logs")
dir_rds <- file.path(dir_results, "rds")
for (directory in c(dir_figures, dir_tables, dir_logs, dir_rds)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

# Always calculate a fresh result and retain its RDS intermediate for inspection.
compute_and_save <- function(path, compute_fun) {
  message("[Processing] ", basename(path))
  result <- compute_fun()
  saveRDS(result, path)
  result
}

# Set plot format --------------------------------------------------------------
save_plot <- function(plot, filename, width = 8, height = 6) {
  output_path <- file.path(dir_figures, paste0(filename, ".", figure_format))
  if (identical(figure_format, "png")) {
    ggplot2::ggsave(
      output_path, plot = plot, width = width, height = height, dpi = 300
    )
  } else {
    ggplot2::ggsave(output_path, plot = plot, width = width, height = height)
  }
}

save_base_plot <- function(plot_fun, filename, width = 8, height = 6) {
  output_path <- file.path(dir_figures, paste0(filename, ".", figure_format))
  if (identical(figure_format, "png")) {
    grDevices::png(output_path, width = width, height = height, units = "in", res = 300)
  } else {
    grDevices::pdf(output_path, width = width, height = height)
  }
  on.exit(grDevices::dev.off(), add = TRUE)
  plot_fun()
}

# Extract samples IDs from FASTQ -----------------------------------------------
sample_id_from_filename <- function(filenames) {
  make.unique(sapply(strsplit(basename(filenames), "_"), `[`, 1))
}


# Metadata loading -------------------------------------------------------------
if (tolower(tools::file_ext(metadata_file)) != "csv") {
  stop("Metadata file must have a .csv extension: ", metadata_file)
}
if (!file.exists(metadata_file)) {
  stop("Metadata file not found: ", metadata_file)
}

metadata <- read.table(
  metadata_file, header = TRUE, sep = ",",
  colClasses = "character", fileEncoding = "UTF-8-BOM",
  stringsAsFactors = FALSE, check.names = FALSE,
  quote = "\"", comment.char = ""
)

if (!"SampleID" %in% colnames(metadata)) {
  stop("Column 'SampleID' was not found in metadata file: ", metadata_file)
}

if (run_stage_max >= 2L && !group_column %in% colnames(metadata)) {
  stop("Column '", group_column, "' was not found in metadata file: ", metadata_file)
}

if (run_stage_max >= 3L && !reference_group %in% metadata[[group_column]]) {
  stop("Reference group '", reference_group,
       "' was not found in column '", group_column, "'.")
}
