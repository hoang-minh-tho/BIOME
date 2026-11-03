# 15. RELATIVE ABUNDANCE
# Computed on the original (non-rarefied) phyloseq object

suppressPackageStartupMessages({
  library(phyloseq)
  library(tidyverse)
  library(Polychrome)
})

physeq <- readRDS(file.path(dir_rds, "physeq.rds"))

metadata_df <- as(sample_data(physeq), "data.frame")
if (!(group_column %in% colnames(metadata_df))) {
  stop(
    "Column '", group_column,
    "' was not found in the selected metadata file.\nAvailable columns: ",
    paste(colnames(metadata_df), collapse = ", ")
  )
}

make_relative_abundance <- function(taxrank, threshold) {
  if (!(taxrank %in% rank_names(physeq))) {
    stop("The selected taxonomy profile does not contain rank '", taxrank, "'.")
  }
  assignments <- as.character(tax_table(physeq)[, taxrank])
  if (!any(!is.na(assignments) & nzchar(trimws(assignments)))) {
    reason <- paste("No assigned taxa at rank", taxrank)
    warning("Skipping relative abundance plot: ", reason)
    empty <- data.frame(Sample = character(), Abundance = numeric())
    empty[[group_column]] <- character()
    empty[[taxrank]] <- character()
    empty <- empty[, c("Sample", group_column, taxrank, "Abundance")]
    attr(empty, "skip_reason") <- reason
    return(empty)
  }
  other_label <- "Others"
  physeq %>%
    tax_glom(taxrank = taxrank, NArm = TRUE) %>%
    transform_sample_counts(function(x) x / sum(x)) %>%
    psmelt() %>%
    group_by(.data[[taxrank]]) %>%
    mutate(mean_relative_abundance = mean(Abundance, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(
      !!taxrank := ifelse(
        mean_relative_abundance < threshold,
        other_label,
        as.character(.data[[taxrank]])
      )
    ) %>%
    group_by(Sample, .data[[group_column]], .data[[taxrank]]) %>%
    summarise(Abundance = sum(Abundance), .groups = "drop")
}

build_taxon_palette <- function(taxon_names, other_label, other_color = "#D3D3D3") {
  unique_taxa <- unique(as.character(taxon_names))
  if (other_label %in% unique_taxa) {
    regular_taxa <- sort(setdiff(unique_taxa, other_label))
    set.seed(SEED)
    colors <- Polychrome::createPalette(
      max(3, length(regular_taxa)),
      seedcolors = c("#EE3333", "#3377EE", "#33AA33"),
      range = c(20, 60)
    )[seq_along(regular_taxa)]
    names(colors) <- regular_taxa
    c(colors, setNames(other_color, other_label))
  } else {
    ordered_taxa <- sort(unique_taxa)
    set.seed(SEED)
    palette <- Polychrome::createPalette(
      max(3, length(ordered_taxa)),
      seedcolors = c("#EE3333", "#3377EE", "#33AA33"),
      range = c(15, 85)
    )[seq_along(ordered_taxa)]
    names(palette) <- ordered_taxa
    palette
  }
}

abundance_barplot <- function(data, taxrank, palette, threshold) {
  legend_ncol <- min(6L, ceiling(length(palette) / 25))
  percentage <- format(threshold * 100, scientific = FALSE, trim = TRUE, digits = 12)
  other_legend_label <- paste0("Other (<", percentage, "%)")
  ggplot(
    data,
    aes(x = Sample, y = Abundance, fill = .data[[taxrank]])
  ) +
    geom_bar(stat = "identity", width = 0.9) +
    facet_wrap(
      as.formula(paste0("~", group_column)),
      nrow = 1,
      scales = "free_x"
    ) +
    scale_fill_manual(
      values = palette,
      breaks = names(palette),
      labels = function(values) ifelse(values == "Others", other_legend_label, values),
      name = taxrank
    ) +
    scale_y_continuous(
      labels = scales::percent,
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(x = "Sample", y = "Relative abundance") +
    guides(fill = guide_legend(ncol = legend_ncol, byrow = FALSE)) +
    theme_bw(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.text = element_text(size = 9),
      legend.key.size = unit(0.4, "cm"),
      plot.margin = margin(10, 20, 10, 10)
    )
}

rank1 <- relative_abundance_ranks[1]
rank2 <- relative_abundance_ranks[2]
rank1_abundance <- make_relative_abundance(rank1, min_rank1_abundance)
rank2_abundance <- make_relative_abundance(rank2, min_rank2_abundance)
rank1_label <- tolower(gsub("[^A-Za-z0-9]+", "_", rank1))
rank2_label <- tolower(gsub("[^A-Za-z0-9]+", "_", rank2))

write.csv(
  rank1_abundance,
  file.path(dir_tables, paste0("relative_abundance_", rank1_label, ".csv")),
  row.names = FALSE
)
write.csv(
  rank2_abundance,
  file.path(dir_tables, paste0("relative_abundance_", rank2_label, ".csv")),
  row.names = FALSE
)

abundance_outputs <- list(rank1_abundance, rank2_abundance)
abundance_status <- lapply(seq_along(abundance_outputs), function(index) {
  data <- abundance_outputs[[index]]
  taxrank <- relative_abundance_ranks[index]
  reason <- attr(data, "skip_reason")
  if (!nrow(data)) {
    return(data.frame(rank = taxrank, status = "skipped", reason = reason))
  }
  palette <- build_taxon_palette(data[[taxrank]], "Others")
  data[[taxrank]] <- factor(data[[taxrank]], levels = names(palette))
  n_samples <- n_distinct(data$Sample)
  n_taxa <- length(palette)
  legend_columns <- min(6L, ceiling(n_taxa / 25))
  requested_width <- max(if (index == 1L) 10 else 11,
                         n_samples * 0.5 + legend_columns * 2.8)
  requested_height <- max(6.5, ceiling(n_taxa / legend_columns) * 0.18)
  # Keep raster memory use bounded and remain below ggsave's 50-inch limit.
  width <- min(24, requested_width)
  height <- min(20, requested_height)
  plot_warning <- ""
  if (requested_width > width || requested_height > height || n_taxa > 100L) {
    plot_warning <- sprintf(
      paste0("%s relative abundance contains %d taxa across %d samples. ",
             "The plot may be crowded; export size is limited to %.1f x %.1f inches. ",
             "Increase the 'Collapse ... taxa below' percentage to simplify the legend. ",
             "0%% keeps all assigned taxa. Your selected thresholds and CSV values are unchanged."),
      taxrank, n_taxa, n_samples, width, height
    )
    message("BIOME_PLOT_WARNING: ", plot_warning)
  }
  label <- tolower(gsub("[^A-Za-z0-9]+", "_", taxrank))
  threshold <- c(min_rank1_abundance, min_rank2_abundance)[index]
  save_plot(abundance_barplot(data, taxrank, palette, threshold),
            paste0("relative_abundance_", label), width = width, height = height)
  data.frame(rank = taxrank, status = "run", reason = "",
             width_inches = width, height_inches = height,
             taxa = n_taxa, samples = n_samples, plot_warning = plot_warning)
})
write.csv(bind_rows(abundance_status),
          file.path(dir_tables, "relative_abundance_status.csv"), row.names = FALSE)
