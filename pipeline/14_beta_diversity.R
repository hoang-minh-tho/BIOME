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
geometry_results <- list()

inspect_distance_geometry <- function(distance, relative_tolerance = 1e-8) {
  squared <- as.matrix(distance)^2
  centered <- sweep(sweep(squared, 1, rowMeans(squared), "-"), 2, colMeans(squared), "-") + mean(squared)
  values <- eigen(-0.5 * centered, symmetric = TRUE, only.values = TRUE)$values
  tolerance <- max(abs(values)) * relative_tolerance
  list(eigenvalues = values, tolerance = tolerance, negative = sum(values < -tolerance))
}

prepare_beta_distance <- function(distance, label) {
  raw <- inspect_distance_geometry(distance)
  use_sqrt <- raw$negative > 0L
  used <- if (use_sqrt) sqrt(distance) else distance
  checked <- if (use_sqrt) inspect_distance_geometry(used) else raw
  message(label, ": ", if (use_sqrt) "applying square-root distances" else "preserving original distances",
          " (", raw$negative, " negative eigenvalues beyond numerical tolerance).")
  if (use_sqrt && checked$negative > 0L) {
    warning(label, ": square-root distances still have ", checked$negative,
            " negative eigenvalues beyond tolerance; geometry is not fully Euclidean. ",
            "No additional correction was applied.")
  }
  list(distance = used, sqrt_dist = use_sqrt,
       diagnostics = data.frame(
         distance = label, sqrt_dist = use_sqrt,
         raw_min_eigenvalue = min(raw$eigenvalues), raw_tolerance = raw$tolerance,
         raw_negative_eigenvalues = raw$negative,
         analysis_min_eigenvalue = min(checked$eigenvalues), analysis_tolerance = checked$tolerance,
         analysis_negative_eigenvalues = checked$negative
       ))
}

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
  
  geometry <- prepare_beta_distance(dist_matrix, distance_labels[[dist_method]])
  geometry_results[[dist_method]] <- geometry$diagnostics
  saveRDS(geometry$distance, file.path(dir_rds, paste0("dist_", dist_method, "_analysis.rds")))
  write.csv(bind_rows(geometry_results), file.path(dir_tables, "beta_distance_geometry.csv"), row.names = FALSE)

  # Distance-based ordination --------------------------------------------------
  if (ordination_method != "PCA") {
    save_beta_ordination(
      ordination_method, paste0(tolower(ordination_method), "_", dist_method),
      paste(ordination_method, "-", distance_labels[[dist_method]],
            if (geometry$sqrt_dist) "(square-root distances)" else ""), geometry$distance
    )
  }
  
  if (run_beta_statistics) {
    # Test homogeneity of dispersion before PERMANOVA.
    if (run_betadisper) {
      disper_fit <- betadisper(dist_matrix, metadata_df[[group_column]], sqrt.dist = geometry$sqrt_dist)
      disper_test <- permutest(disper_fit, permutations = permanova_permutations)
      betadisper_results[[dist_method]] <-
        as.data.frame(disper_test$tab) %>%
        rownames_to_column("term") %>%
        mutate(distance = distance_labels[[dist_method]], sqrt_dist = geometry$sqrt_dist, .before = 1)
    }

    permanova_fit <- adonis2(
      as.formula(paste("dist_matrix ~", group_column)),
      data = metadata_df,
      permutations = permanova_permutations,
      sqrt.dist = geometry$sqrt_dist
    )

    permanova_results[[dist_method]] <-
      as.data.frame(permanova_fit) %>%
      rownames_to_column("term") %>%
      mutate(distance = distance_labels[[dist_method]], sqrt_dist = geometry$sqrt_dist, .before = 1)
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
