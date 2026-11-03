# 13. ALPHA DIVERSITY

suppressPackageStartupMessages({
  library(phyloseq)
  library(tidyverse)
})

# Estimate Richness & Merge Metadata -------------------------------------------

physeq_rarefied <- readRDS(file.path(dir_rds, "physeq_rarefied.rds"))

# Estimate alpha diversity metrics
alpha_div <- estimate_richness(physeq_rarefied, measures = alpha_measures)
alpha_div$SampleID <- rownames(alpha_div)

# Merge sample metadata
sample_df <- as(sample_data(physeq_rarefied), "data.frame")
sample_df$SampleID <- rownames(sample_df)

alpha_div <- alpha_div %>%
  left_join(sample_df, by = "SampleID") %>%
  relocate(SampleID, all_of(group_column))

# Remove samples with missing group information
if (any(is.na(alpha_div[[group_column]]))) {
  n_na <- sum(is.na(alpha_div[[group_column]]))
  warning(sprintf("Removing %d sample(s) with NA in '%s'.", n_na, group_column))
  alpha_div <- alpha_div[!is.na(alpha_div[[group_column]]), ]
}

alpha_div[[group_column]] <- droplevels(factor(alpha_div[[group_column]]))
present_group_levels <- levels(alpha_div[[group_column]])
n_groups <- length(present_group_levels)

# Save sample-level alpha diversity table
write.csv(alpha_div, file.path(dir_tables, "alpha_diversity.csv"), row.names = FALSE)
message(sprintf("Alpha diversity calculated for %d samples across %d groups.", nrow(alpha_div), n_groups))

# Statistical Testing ----------------------------------------------------------

p_labels <- list() # Store p-values for plot facet headers

if (n_groups < 2) {
  
  message("Only 1 group present. Skipping hypothesis testing.")
  
} else if (n_groups == 2) {
  
  # --- Case A: Exactly 2 Groups -> Wilcoxon Rank-Sum Test --------------------
  message(sprintf("Running Wilcoxon rank-sum tests (%s vs %s)...",
                  present_group_levels[1], present_group_levels[2]))
  
  wilcox_results <- lapply(alpha_measures, function(metric) {
    test <- suppressWarnings(
      wilcox.test(reformulate(group_column, response = metric), data = alpha_div)
    )
    data.frame(
      Metric     = metric,
      Group1     = present_group_levels[1],
      Group2     = present_group_levels[2],
      Statistic  = unname(test$statistic),
      p_value    = test$p.value
    )
  })
  
  wilcox_results <- dplyr::bind_rows(wilcox_results) %>%
    mutate(
      p_adj = p.adjust(p_value, method = "BH"),
      Significance = case_when(
        p_adj < 0.001 ~ "***",
        p_adj < 0.01  ~ "**",
        p_adj < 0.05  ~ "*",
        TRUE          ~ "ns"
      )
    )
  
  write.csv(wilcox_results, file.path(dir_tables, "alpha_diversity_wilcoxon.csv"), row.names = FALSE)
  print(wilcox_results)
  
  # Format facet labels with p-value
  p_labels <- setNames(
    sprintf("%s\n(Wilcoxon p = %.3g)", wilcox_results$Metric, wilcox_results$p_value),
    wilcox_results$Metric
  )
  
} else {
  
  # --- Case B: > 2 Groups -> Kruskal-Wallis + Pairwise Wilcoxon Post-Hoc ------
  message(sprintf("Running Kruskal-Wallis omnibus tests across %d groups...", n_groups))
  
  # 1. Omnibus Test: Kruskal-Wallis
  kw_results <- lapply(alpha_measures, function(metric) {
    test <- suppressWarnings(
      kruskal.test(reformulate(group_column, response = metric), data = alpha_div)
    )
    data.frame(
      Metric    = metric,
      ChiSq     = unname(test$statistic),
      df        = unname(test$parameter),
      p_value   = test$p.value
    )
  })
  
  kw_results <- dplyr::bind_rows(kw_results) %>%
    mutate(
      p_adj = p.adjust(p_value, method = "BH"),
      Significance = case_when(
        p_adj < 0.001 ~ "***",
        p_adj < 0.01  ~ "**",
        p_adj < 0.05  ~ "*",
        TRUE          ~ "ns"
      )
    )
  
  write.csv(kw_results, file.path(dir_tables, "alpha_diversity_kruskal_omnibus.csv"), row.names = FALSE)
  message("--- Kruskal-Wallis Omnibus Results ---")
  print(kw_results)
  
  # 2. Post-Hoc Test: Pairwise Wilcoxon with FDR Correction
  message("Running pairwise Wilcoxon post-hoc comparisons (FDR-adjusted)...")
  pairs <- combn(present_group_levels, 2, simplify = FALSE)
  
  pairwise_results <- lapply(alpha_measures, function(metric) {
    pair_rows <- lapply(pairs, function(p) {
      sub_df <- alpha_div %>% filter(.data[[group_column]] %in% p)
      test <- suppressWarnings(
        wilcox.test(reformulate(group_column, response = metric), data = sub_df)
      )
      data.frame(
        Metric    = metric,
        Group1    = p[1],
        Group2    = p[2],
        W         = unname(test$statistic),
        p_value   = test$p.value
      )
    })
    
    # Adjust p-values across all pairs within this metric
    dplyr::bind_rows(pair_rows) %>%
      mutate(
        p_adj = p.adjust(p_value, method = "BH"),
        Significance = case_when(
          p_adj < 0.001 ~ "***",
          p_adj < 0.01  ~ "**",
          p_adj < 0.05  ~ "*",
          TRUE          ~ "ns"
        )
      )
  })
  
  pairwise_results <- dplyr::bind_rows(pairwise_results)
  write.csv(pairwise_results, file.path(dir_tables, "alpha_diversity_pairwise_wilcoxon.csv"), row.names = FALSE)
  message("--- Pairwise Post-Hoc Results (Top 10) ---")
  print(head(pairwise_results, 10))
  
  # Format facet labels with omnibus Kruskal-Wallis p-value
  p_labels <- setNames(
    sprintf("%s\n(Kruskal-Wallis p = %.3g)", kw_results$Metric, kw_results$p_value),
    kw_results$Metric
  )
}

# Alpha Diversity Boxplot ------------------------------------------------------

alpha_long <- alpha_div %>%
  pivot_longer(cols = all_of(alpha_measures), names_to = "Metric", values_to = "Value")

# Apply updated labels with p-values to facets (if statistical test was performed)
if (length(p_labels) > 0) {
  alpha_long$FacetLabel <- factor(p_labels[alpha_long$Metric], levels = unname(p_labels))
} else {
  alpha_long$FacetLabel <- factor(alpha_long$Metric)
}

alpha_plot <- ggplot(alpha_long, aes(x = .data[[group_column]], y = Value, fill = .data[[group_column]])) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1.8, shape = 21, color = "black") +
  facet_wrap(~FacetLabel, scales = "free_y") +
  labs(
    x = group_column,
    y = "Diversity Metric Value",
    fill = group_column
  ) +
  theme_bw(base_size = 13) +
  theme(
    strip.background = element_rect(fill = "#f8f9fa", color = "#dcdde1"),
    strip.text       = element_text(face = "bold", size = 11),
    axis.text.x      = element_text(angle = if (n_groups > 3) 45 else 0, hjust = if (n_groups > 3) 1 else 0.5),
    legend.position  = "none" # Redundant with x-axis labels
  )

save_plot(alpha_plot, "alpha_diversity_boxplots")
