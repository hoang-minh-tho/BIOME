# 15. RELATIVE ABUNDANCE
# Computed on the original (non-rarefied) phyloseq object

suppressPackageStartupMessages({
  library(phyloseq)
  library(tidyverse)
  library(Polychrome)
})

physeq <- readRDS(file.path(dir_rds, "physeq.rds"))

selected_plot_type <- if (exists("abundance_plot_type")) abundance_plot_type else "barplot"
metadata_df <- as(sample_data(physeq), "data.frame")
if (selected_plot_type %in% c("barplot", "both") && !(group_column %in% colnames(metadata_df))) {
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

prepare_abundance_heatmap <- function(physeq, taxrank, top_n, group_column) {
  if (!taxrank %in% rank_names(physeq)) {
    stop("The selected taxonomy profile does not contain rank '", taxrank, "'.")
  }
  metadata <- data.frame(sample_data(physeq), check.names = FALSE)
  if (!group_column %in% colnames(metadata)) stop("Missing comparison group: ", group_column)
  # Preserve lineage during aggregation, including unclassified taxa. Equal
  # terminal names under different parents must remain distinct heatmap rows.
  taxonomy <- as(tax_table(physeq), "matrix")
  taxonomy[] <- trimws(taxonomy)
  taxonomy[is.na(taxonomy) | taxonomy == ""] <- "Unclassified"
  tax_table(physeq) <- tax_table(taxonomy)
  physeq <- prune_taxa(taxa_sums(physeq) > 0, physeq)
  if (!ntaxa(physeq)) stop("No nonzero taxa are available for the heatmap.")
  aggregated <- tax_glom(physeq, taxrank = taxrank, NArm = FALSE)
  aggregated <- transform_sample_counts(aggregated, function(x) {
    total <- sum(x)
    if (total > 0) x / total else rep(0, length(x))
  })
  abundance <- as(otu_table(aggregated), "matrix")
  if (!taxa_are_rows(aggregated)) abundance <- t(abundance)
  # Rank once across all samples; do not renormalize the selected subset.
  totals <- rowSums(abundance)
  selected <- head(order(-totals, rownames(abundance)), min(top_n, nrow(abundance)))
  all_taxon_ids <- rownames(abundance)
  labels <- as.character(tax_table(aggregated)[all_taxon_ids, taxrank])
  duplicated_labels <- duplicated(labels) | duplicated(labels, fromLast = TRUE)
  labels[duplicated_labels] <- paste0(labels[duplicated_labels], " [", all_taxon_ids[duplicated_labels], "]")
  labels <- make.unique(labels)
  abundance <- abundance[selected, , drop = FALSE]
  taxon_ids <- rownames(abundance)
  rownames(abundance) <- labels[selected]
  z_scores <- sweep(abundance, 1, rowMeans(abundance), "-")
  row_sd <- apply(abundance, 1, stats::sd)
  row_sd[!is.finite(row_sd) | row_sd == 0] <- 1
  z_scores <- sweep(z_scores, 1, row_sd, "/")
  groups <- as.character(metadata[colnames(abundance), group_column])
  groups[is.na(groups) | !nzchar(trimws(groups))] <- "Unclassified"
  annotation <- data.frame(factor(groups, levels = sort(unique(groups))),
                           row.names = colnames(abundance), check.names = FALSE)
  colnames(annotation) <- group_column
  levels <- levels(annotation[[group_column]])
  colors <- setNames(grDevices::hcl.colors(length(levels), "Dark 3"), levels)
  list(abundance = abundance, z_scores = z_scores, annotation = annotation,
       annotation_colors = setNames(list(colors), group_column), taxon_ids = taxon_ids)
}

stack_heatmap_legends <- function(heatmap) {
  layout <- heatmap$layout
  legend_indices <- match(c("legend", "annotation_legend"), layout$name)
  if (anyNA(legend_indices)) return(heatmap)
  legends <- heatmap$grobs[legend_indices]
  # Measure the actual text on the active output device, rather than estimating
  # widths from character counts (wide glyphs and bold headings need more room).
  text_widths <- unlist(lapply(legends, function(legend) {
    vapply(Filter(function(child) inherits(child, "text"), as.list(legend$children)),
           function(child) grid::convertWidth(child$x + grid::grobWidth(child),
                                              "inches", valueOnly = TRUE), numeric(1))
  }))
  legend_width <- max(text_widths) + 0.25
  legend_heights <- vapply(legends, function(legend) {
    max(vapply(as.list(legend$children), function(child) {
      # Legend children are positioned downwards from the top of their viewport.
      offset <- grid::unit(1, "npc") - min(child$y)
      if (inherits(child, "text")) {
        # Measure glyph height without the labels' different vertical positions.
        offset <- offset + grid::grobHeight(grid::textGrob(child$label, gp = child$gp))
      }
      if (inherits(child, "rect")) offset <- offset + max(child$height)
      grid::convertHeight(offset, "inches", valueOnly = TRUE)
    }, numeric(1)))
  }, numeric(1))
  stacked <- gtable::gtable(
    widths = grid::unit(legend_width, "inches"),
    heights = grid::unit(c(0.12, legend_heights[1], 0.25, legend_heights[2], 0.12), "inches")
  )
  stacked <- gtable::gtable_add_grob(stacked, legends[[1]], t = 2, l = 1, clip = "off")
  stacked <- gtable::gtable_add_grob(stacked, legends[[2]], t = 4, l = 1, clip = "off")
  legend_column <- layout$l[legend_indices[1]]
  annotation_column <- layout$l[legend_indices[2]]
  matrix_column <- layout$l[match("matrix", layout$name)]
  heatmap$grobs <- heatmap$grobs[-legend_indices]
  heatmap$layout <- layout[-legend_indices, , drop = FALSE]
  heatmap$widths[legend_column] <- grid::unit(legend_width, "inches")
  heatmap$widths[annotation_column] <- grid::unit(0, "inches")
  heatmap$widths[matrix_column] <- grid::unit(1, "npc") - sum(heatmap$widths[-matrix_column])
  # Put the stacked legends in the same viewport; leave room for top/bottom text.
  heatmap <- gtable::gtable_add_grob(
    heatmap, stacked, t = min(layout$t[legend_indices]),
    b = max(layout$b[legend_indices]), l = legend_column,
    clip = "off", name = "stacked_legends"
  )
  heatmap
}

save_abundance_heatmap <- function(data, taxrank, top_n) {
  limit <- max(abs(data$z_scores))
  if (!is.finite(limit) || limit == 0) limit <- 1
  n_samples <- ncol(data$z_scores)
  n_taxa <- nrow(data$z_scores)
  annotation_labels <- c(names(data$annotation_colors),
                         unlist(lapply(data$annotation_colors, names), use.names = FALSE))
  # Reserve space before opening the device; the layout then uses exact metrics.
  width <- min(24, max(10, n_samples * 0.18 + 5,
                      max(nchar(annotation_labels)) * 10 / 72 + 6))
  n_groups <- sum(lengths(data$annotation_colors))
  height <- min(24, max(7, n_taxa * 0.19 + 3, (150 + (2 * n_groups + 5) * 10) / 72 + 2))
  filename <- paste0("heatmap_top", top_n, "_", tolower(taxrank))
  # Reuse the existing device helper: PNG is 300 DPI; PDF is vector output.
  save_base_plot(function() {
    heatmap <- pheatmap::pheatmap(
      data$z_scores, scale = "none",
      color = grDevices::colorRampPalette(c("#4575B4", "#FFFFBF", "#D73027"))(101),
      breaks = seq(-limit, limit, length.out = 102),
      cluster_rows = n_taxa > 1L, cluster_cols = n_samples > 1L,
      annotation_col = data$annotation, annotation_colors = data$annotation_colors,
      show_rownames = TRUE, show_colnames = n_samples <= 60L,
      fontsize = 10, fontsize_row = 9, fontsize_col = 8, angle_col = "45",
      border_color = NA, main = paste0(taxrank, " relative abundance (row Z-scores)"),
      legend_breaks = c(-limit, 0, limit),
      legend_labels = c(sprintf("%.1f", -limit), "0 (row Z-score)", sprintf("%.1f", limit)),
      silent = TRUE
    )
    grid::grid.newpage()
    grid::grid.draw(stack_heatmap_legends(heatmap$gtable))
  }, filename, width = width, height = height)
  invisible(filename)
}

# Older configurations and direct stage callers retain the barplot workflow.
if (selected_plot_type %in% c("heatmap", "both")) {
  selected_rank <- if (exists("heatmap_rank")) heatmap_rank else relative_abundance_ranks[2]
  selected_top_n <- if (exists("heatmap_top_n")) heatmap_top_n else 30L
  heatmap_status <- tryCatch({
    heatmap_data <- prepare_abundance_heatmap(physeq, selected_rank, selected_top_n, group_column)
    save_abundance_heatmap(heatmap_data, selected_rank, selected_top_n)
    data.frame(rank = selected_rank, top_n = selected_top_n, status = "run", reason = "")
  }, error = function(error) {
    reason <- conditionMessage(error)
    message("ERROR: Abundance heatmap failed (", selected_rank, "): ", reason,
            ". Continuing remaining pipeline steps.")
    data.frame(rank = selected_rank, top_n = selected_top_n, status = "failed", reason = reason)
  })
  write.csv(heatmap_status, file.path(dir_tables, "heatmap_status.csv"), row.names = FALSE)
}
if (selected_plot_type %in% c("barplot", "both")) {
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
}
