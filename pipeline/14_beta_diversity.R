# 14. BETA DIVERSITY

suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(tidyverse)
})

physeq_rarefied <- readRDS(file.path(dir_rds, "physeq_rarefied.rds"))

# Extract metadata from phyloseq -----------------------------------------------
metadata_df <- as(sample_data(physeq_rarefied), "data.frame")
metadata_df$SampleID <- rownames(metadata_df)

# Remove samples with missing group information
metadata_df <- metadata_df %>%
  filter(!is.na(.data[[group_column]]))

# Keep phyloseq in the same sample order
physeq_rarefied <- prune_samples(metadata_df$SampleID, physeq_rarefied)
metadata_df <- metadata_df[match(sample_names(physeq_rarefied), metadata_df$SampleID), ]

metadata_df[[group_column]] <- droplevels(factor(metadata_df[[group_column]]))

distance_labels <- c(
  bray      = "Bray-Curtis",
  wunifrac  = "Weighted UniFrac",
  unifrac   = "Unweighted UniFrac",
  jaccard   = "Jaccard"
)

permanova_results  <- list()
betadisper_results <- list()

# Rarefaction may remove a whole group or leave too few replicates for testing.
group_counts <- table(metadata_df[[group_column]])
write.csv(
  data.frame(group = names(group_counts), samples = as.integer(group_counts)),
  file.path(dir_tables, "beta_group_counts.csv"), row.names = FALSE
)
statistics_reason <- if (length(group_counts) < 2L) {
  "Fewer than two groups remain after rarefaction."
} else if (any(group_counts < 2L)) {
  "At least one group has fewer than two samples after rarefaction."
} else {
  ""
}
run_beta_statistics <- run_permanova && !nzchar(statistics_reason)
if (run_permanova && nzchar(statistics_reason)) {
  warning("Skipping PERMANOVA and betadisper: ", statistics_reason)
}
if (run_permanova) {
  write.csv(
    data.frame(status = if (run_beta_statistics) "run" else "skipped",
               reason = statistics_reason),
    file.path(dir_tables, "beta_statistics_status.csv"), row.names = FALSE
  )
}

# Sparse retained data may not support a two-dimensional ordination plot.
save_beta_ordination <- function(method, filename, title, distance = NULL) {
  if (nsamples(physeq_rarefied) < 3L) {
    warning("Skipping ", title, ": fewer than three samples remain after rarefaction.")
    return(invisible(NULL))
  }
  tryCatch({
    ordination <- if (method == "PCA") {
      ordinate(physeq_rarefied, method = "RDA", formula = ~1)
    } else {
      ordinate(physeq_rarefied, method = method, distance = distance)
    }
    plot <- plot_ordination(physeq_rarefied, ordination, color = group_column) +
      stat_ellipse(type = "norm", linetype = 2) +
      labs(title = title) + theme_bw(base_size = 13)
    save_plot(plot, filename)
  }, error = function(error) {
    warning("Skipping ", title, ": ", conditionMessage(error))
  })
}

# PCA uses the abundance table, independently of the selected distances.
if (ordination_method == "PCA") {
  save_beta_ordination("PCA", "pca_asv_counts", "PCA of rarefied ASV counts")
}

for (dist_method in beta_distances) {
  
  dist_matrix <- compute_and_save(
  file.path(dir_rds, paste0("dist_", dist_method, ".rds")),
  function() {
      phyloseq::distance(physeq_rarefied, method = dist_method)
    }
)
  
  # Distance-based ordination --------------------------------------------------
  if (ordination_method != "PCA") {
    save_beta_ordination(
      ordination_method, paste0(tolower(ordination_method), "_", dist_method),
      paste(ordination_method, "-", distance_labels[[dist_method]]), dist_matrix
    )
  }
  
  if (run_beta_statistics) {
    # Test homogeneity of dispersion before PERMANOVA.
    if (run_betadisper) {
      disper_fit <- betadisper(dist_matrix, metadata_df[[group_column]])
      disper_test <- permutest(disper_fit, permutations = permanova_permutations)
      betadisper_results[[dist_method]] <-
        as.data.frame(disper_test$tab) %>%
        rownames_to_column("term") %>%
        mutate(distance = distance_labels[[dist_method]], .before = 1)
    }

    permanova_fit <- adonis2(
      as.formula(paste("dist_matrix ~", group_column)),
      data = metadata_df,
      permutations = permanova_permutations
    )

    permanova_results[[dist_method]] <-
      as.data.frame(permanova_fit) %>%
      rownames_to_column("term") %>%
      mutate(distance = distance_labels[[dist_method]], .before = 1)
  }
}

if (run_beta_statistics) {
  write.csv(
    bind_rows(permanova_results),
    file.path(dir_tables, "permanova_results.csv"),
    row.names = FALSE
  )

  if (run_betadisper) {
    write.csv(bind_rows(betadisper_results),
              file.path(dir_tables, "betadisper_results.csv"), row.names = FALSE)
  }
}
