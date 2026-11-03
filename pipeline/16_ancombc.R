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
  ancombc_fit <- compute_and_save(
  file.path(dir_rds, "ancombc_fit.rds"),
  function() ancombc2(
      data = physeq, tax_level = relative_abundance_ranks[2],
      fix_formula = group_column,
      rand_formula = NULL, p_adj_method = "holm", group = group_column,
      struc_zero = TRUE, neg_lb = TRUE, alpha = 0.05, global = FALSE
    )
)
  results <- ancombc_fit$res
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
      taxon, lfc = .data[[lfc_col]], q = .data[[q_col]], significant = q < 0.05
    )
    save_plot(
      ggplot(volcano_data, aes(x = lfc, y = -log10(pmax(q, .Machine$double.xmin)), color = significant)) +
        geom_point(alpha = 0.7) +
        labs(title = paste0(coefficient, " (reference: ", reference_group, ")"),
             x = "Log fold change (natural log)",
             y = expression(-log[10](q)), color = "q < 0.05") +
        theme_bw(base_size = 13),
      paste0("ancombc_volcano_", contrast_index, "_",
             gsub("[^A-Za-z0-9]+", "_", coefficient))
    )
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
