# 16. ANCOM-BC

suppressPackageStartupMessages({
  library(phyloseq)
  library(tidyverse)
})

physeq <- readRDS(file.path(dir_rds, "physeq.rds"))

# Check metadata ---------------------------------------------------------------
meta <- as(sample_data(physeq), "data.frame")

if (!(group_column %in% colnames(meta))) {
  stop("Column '", group_column, "' not found in the selected metadata file: ", metadata_file)
}

if (!(reference_group %in% meta[[group_column]])) {
  stop("Reference group '", reference_group,
       "' not found in metadata column '", group_column, "'.")
}

sample_data(physeq)[[group_column]] <-
  relevel(factor(sample_data(physeq)[[group_column]]),
          ref = reference_group)

# Match ANCOM-BC2's aggregation and filtering before fitting. In ANCOMBC
# 2.14.0, a single modeled taxon collapses q_hat to a vector, then the
# pseudo-count sensitivity analysis fails at apply(ss_3d_prim, c(1, 2), ...).
inspect_ancombc_taxa <- function(data, rank, group, prv_cut = 0.1) {
  checked <- ANCOMBC::data_sanity_check(
    data = data, tax_level = rank, fix_formula = group, group = group,
    struc_zero = TRUE, global = FALSE, verbose = FALSE
  )
  aggregated <- checked$feature_table_aggregate
  # Reuse the package's structural-zero detector so this guard agrees with
  # the model's struc_zero = TRUE and neg_lb = TRUE settings.
  zero_ind <- getFromNamespace(".get_struc_zero", "ANCOMBC")(
    data = aggregated, meta_data = checked$meta_data, group = group, neg_lb = TRUE
  )
  structural <- rowSums(as.matrix(zero_ind[, -1, drop = FALSE])) > 0L
  prevalence <- rowMeans(aggregated != 0, na.rm = TRUE)
  eligible <- !structural & is.finite(prevalence) & prevalence >= prv_cut
  list(
    zero_ind = zero_ind,
    filtering = data.frame(
      taxon = rownames(aggregated), prevalence = prevalence,
      structural_zero_any_group = structural, eligible_for_model = eligible,
      row.names = NULL
    ),
    modeled_taxa = sum(eligible)
  )
}

if (differential_method == "ANCOM-BC") {
  differential_rank <- relative_abundance_ranks[2]
  assignments <- as.character(tax_table(physeq)[, differential_rank])
  if (!any(!is.na(assignments) & nzchar(trimws(assignments)))) {
    reason <- paste("No assigned taxa at rank", differential_rank)
    warning("Skipping ANCOM-BC: ", reason)
    write.csv(data.frame(method = differential_method, rank = differential_rank,
                         status = "skipped", reason = reason),
              file.path(dir_tables, "differential_abundance_status.csv"),
              row.names = FALSE)
  } else {
  suppressPackageStartupMessages(library(ANCOMBC))
  ancombc_input <- inspect_ancombc_taxa(physeq, differential_rank, group_column)
  write.csv(ancombc_input$zero_ind,
            file.path(dir_tables, "ancombc_structural_zeros.csv"), row.names = FALSE)
  write.csv(ancombc_input$filtering,
            file.path(dir_tables, "ancombc_taxa_filtering.csv"), row.names = FALSE)
  message(sprintf("ANCOM-BC2: %d aggregated taxa; %d eligible after structural-zero and prevalence filtering.",
                  nrow(ancombc_input$filtering), ancombc_input$modeled_taxa))
  if (ancombc_input$modeled_taxa < 2L) {
    reason <- sprintf(
      "Only %d taxa remain at rank %s after structural-zero and prevalence filtering; at least 2 are required by BIOME's ANCOM-BC2 guard. See ancombc_taxa_filtering.csv and ancombc_structural_zeros.csv.",
      ancombc_input$modeled_taxa, differential_rank
    )
    warning("Skipping ANCOM-BC2: ", reason)
    write.csv(data.frame(method = differential_method, rank = differential_rank,
                         status = "skipped", reason = reason,
                         modeled_taxa = ancombc_input$modeled_taxa,
                         package_version = as.character(packageVersion("ANCOMBC"))),
              file.path(dir_tables, "differential_abundance_status.csv"), row.names = FALSE)
  } else {
  ancombc_fit <- compute_and_save(
  file.path(dir_rds, "ancombc_fit.rds"),
  function() ancombc2(
      data = physeq, tax_level = relative_abundance_ranks[2],
      fix_formula = group_column,
      rand_formula = NULL, p_adj_method = "holm", group = group_column,
      struc_zero = TRUE, neg_lb = TRUE, prv_cut = 0.1, lib_cut = 0,
      pseudo_sens = TRUE, alpha = 0.05, global = FALSE
    )
)
  results <- ancombc_fit$res
  write.csv(data.frame(method = differential_method, rank = differential_rank,
                       status = "run", reason = "",
                       modeled_taxa = ancombc_input$modeled_taxa,
                       package_version = as.character(packageVersion("ANCOMBC"))),
            file.path(dir_tables, "differential_abundance_status.csv"), row.names = FALSE)
  write.csv(results, file.path(dir_tables, "ancombc_results.csv"), row.names = FALSE)
  lfc_cols <- grep("^lfc_", names(results), value = TRUE)
  lfc_cols <- lfc_cols[!grepl("intercept", lfc_cols, ignore.case = TRUE)]
  for (contrast_index in seq_along(lfc_cols)) {
    lfc_col <- lfc_cols[contrast_index]
    coefficient <- sub("^lfc_", "", lfc_col)
    q_col <- paste0("q_", coefficient)
    if (!(q_col %in% names(results))) {
      warning("Skipping ANCOM-BC coefficient without a matching q-value: ", coefficient)
      next
    }
    volcano_data <- results %>% transmute(
      taxon, lfc = .data[[lfc_col]], q = .data[[q_col]],
      Significance = factor(significance_label(q),
                            levels = c("ns", "*", "**", "***", "****", "not tested"))
    )
    save_plot(
      ggplot(volcano_data, aes(x = lfc, y = -log10(pmax(q, .Machine$double.xmin)), color = Significance)) +
        geom_point(alpha = 0.7) +
        scale_color_manual(values = c("ns" = "grey65", "*" = "#e69f00",
                                      "**" = "#d55e00", "***" = "#cc3311",
                                      "****" = "#882255", "not tested" = "grey85")) +
        labs(title = paste0(coefficient, " (reference: ", reference_group, ")"),
             x = "Log fold change (natural log)",
             y = expression(-log[10](q)), color = "Significance (adjusted p)",
             caption = "**** <= 0.0001; *** <= 0.001; ** <= 0.01; * <= 0.05; ns > 0.05") +
        theme_bw(base_size = 13),
      paste0("ancombc_volcano_", contrast_index, "_",
             gsub("[^A-Za-z0-9]+", "_", coefficient))
    )
  }
  }
  }
} else {
  differential_results <- compute_and_save(
  file.path(dir_rds, paste0(
      "differential_", tolower(differential_method), "_results.rds"
    )),
  function() {
      meta <- as(sample_data(physeq), "data.frame")
      groups <- droplevels(factor(meta[[group_column]]))
      comparisons <- setdiff(levels(groups), reference_group)
      if (differential_method == "DESeq2") {
        suppressPackageStartupMessages(library(DESeq2))
        dds <- phyloseq_to_deseq2(physeq, as.formula(paste("~", group_column)))
        dds <- DESeq(dds, quiet = TRUE)
        results <- lapply(comparisons, function(comparison) {
          result <- as.data.frame(DESeq2::results(
            dds, contrast = c(group_column, comparison, reference_group)
          ))
          result$taxon <- rownames(result)
          result$contrast <- paste(comparison, "vs", reference_group)
          result
        })
      } else if (differential_method == "ALDEx2") {
        suppressPackageStartupMessages(library(ALDEx2))
        count_matrix <- as(otu_table(physeq), "matrix")
        if (!taxa_are_rows(physeq)) count_matrix <- t(count_matrix)
        groups <- factor(meta[sample_names(physeq), group_column], levels = levels(groups))
        results <- lapply(comparisons, function(comparison) {
          keep <- groups %in% c(reference_group, comparison)
          pair_groups <- droplevels(groups[keep])
          pair_counts <- count_matrix[, sample_names(physeq)[keep], drop = FALSE]
          clr <- aldex.clr(
            pair_counts, pair_groups, mc.samples = 128, denom = "all", verbose = FALSE
          )
          result <- as.data.frame(aldex.tt(clr, paired.test = FALSE, verbose = FALSE))
          result$taxon <- rownames(result)
          result$contrast <- paste(comparison, "vs", reference_group)
          result
        })
      } else {
        stop("Unknown differential abundance method: ", differential_method)
      }
      bind_rows(results)
    }
)
  write.csv(differential_results,
            file.path(dir_tables, paste0(tolower(gsub("[^A-Za-z0-9]+", "_", differential_method)),
                                         "_results.csv")),
            row.names = FALSE)
}
