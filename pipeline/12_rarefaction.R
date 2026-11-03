# 12. RAREFACTION

suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
})

# Load raw phyloseq object -----------------------------------------------------
physeq <- readRDS(file.path(dir_rds, "physeq.rds"))

# Extract sequencing depths and OTU table --------------------------------------
depths  <- sample_sums(physeq)
n_total <- length(depths)

asv_mat <- as(otu_table(physeq), "matrix")
if (taxa_are_rows(physeq)) {
  asv_mat <- t(asv_mat)
}

message(sprintf(
  "Loaded phyloseq: %d samples (depths: %d to %d reads).",
  n_total, min(depths), max(depths)
))

# Depth Selection --------------------------------------------------------------

depth_selection <- compute_and_save(
  file.path(dir_rds, "rarefaction_depth.rds"),
  function() {
    opt_df <- NULL
    if (rarefaction_mode == "min") {
      target_depth <- min(depths)
      message(sprintf(
        "[Mode: min] Selected depth = %d reads (100%% samples kept).",
        target_depth
      ))
    } else if (rarefaction_mode == "custom") {
      target_depth <- custom_rarefaction_depth
      n_kept <- sum(depths >= target_depth)
      message(sprintf(
        "[Mode: custom] Selected depth = %d reads (retains %d/%d samples, %.1f%%).",
        target_depth, n_kept, n_total,
        100 * n_kept / n_total
      ))
    } else if (rarefaction_mode == "auto_rarefy") {
      message("[Mode: auto_rarefy] Evaluating candidate depths...")
      cand_depths <- sort(unique(depths))
      min_samples_to_keep <- ceiling(n_total * (1 - max_sample_loss))
      opt_records <- lapply(cand_depths, function(D) {
        retained_mask <- depths >= D
        n_retained <- sum(retained_mask)
        if (n_retained < min_samples_to_keep) return(NULL)

        retained_idx <- which(retained_mask)
        sub_mat <- asv_mat[retained_idx, , drop = FALSE]
        sub_depths <- depths[retained_idx]
        s_at_D <- vegan::rarefy(sub_mat, sample = D)
        step_targets <- pmin(sub_depths, D + read_step)
        rel_gains <- sapply(seq_along(retained_idx), function(j) {
          if (step_targets[j] <= D || s_at_D[j] <= 0) return(NA_real_)
          s_step <- vegan::rarefy(
            sub_mat[j, , drop = FALSE],
            sample = step_targets[j]
          )
          eff_step <- step_targets[j] - D
          ((s_step - s_at_D[j]) / s_at_D[j]) * (read_step / eff_step)
        })
        assessable <- is.finite(rel_gains)
        n_assessed <- sum(assessable)
        # Require evidence of plateau in the configured fraction of ALL retained
        # samples. Unassessable samples neither demonstrate nor disprove plateau.
        n_plateaued <- sum(rel_gains[assessable] <= plateau_gain)
        pct_plateaued <- if (n_assessed) n_plateaued / n_assessed else NA_real_
        data.frame(
          depth = D,
          samples_kept = n_retained,
          pct_kept = 100 * n_retained / n_total,
          samples_assessed = n_assessed,
          pct_assessed = 100 * n_assessed / n_retained,
          pct_plateaued = 100 * pct_plateaued,
          is_plateaued = n_assessed > 0 &&
            n_plateaued / n_retained >= plateau_samples,
          yield = D * n_retained
        )
      })
      opt_df <- do.call(rbind, opt_records)
      feasible <- opt_df[opt_df$is_plateaued, ]
      if (nrow(feasible) > 0) {
        best_row <- feasible[which.max(feasible$yield), ]
        target_depth <- best_row$depth
        message(sprintf(
          "Auto-yield optimal: %d reads | Retains %d/%d samples (%.1f%%) | %.1f%% at plateau.",
          best_row$depth,
          best_row$samples_kept,
          n_total,
          best_row$pct_kept,
          best_row$pct_plateaued
        ))
      } else {
        warning(
          "Plateau criterion not met. Falling back to maximum yield under the sample-retention constraint."
        )
        best_row <- opt_df[which.max(opt_df$yield), ]
        target_depth <- best_row$depth
        message(sprintf(
          "Fallback depth: %d reads | Retains %d/%d samples (%.1f%%).",
          best_row$depth,
          best_row$samples_kept,
          n_total,
          best_row$pct_kept
        ))
      }
      opt_df$selected <- opt_df$depth == target_depth
    } else {
      stop("Unknown rarefaction mode: ", rarefaction_mode)
    }
    list(target_depth = target_depth, optimization = opt_df)
  }
)
target_depth <- depth_selection$target_depth
opt_df <- depth_selection$optimization
if (!is.null(opt_df)) {
  write.csv(
    opt_df,
    file.path(dir_tables, "rarefaction_optimization.csv"),
    row.names = FALSE
  )
}

# Rarefaction Curve ------------------------------------------------------------

n_kept <- sum(depths >= target_depth)

save_base_plot(
  function() {
    vegan::rarecurve(
      asv_mat,
      step = max(50, round(max(depths) / 100)),
      col = "#34495e",
      lwd = 1,
      lty = 1,
      label = FALSE,
      xlab = "Sequencing Depth (Reads)",
      ylab = "Expected ASVs",
      main = "Rarefaction Curves"
    )

    abline(v = target_depth, col = "#e74c3c", lwd = 2, lty = 2)
    legend(
      "bottomright",
      legend = c(
        sprintf("Selected depth: %s reads", format(target_depth, big.mark = ",")),
        sprintf("Retained: %d/%d samples (%.1f%%)", n_kept, n_total,
                100 * n_kept / n_total)
      ),
      col = c("#e74c3c", NA),
      lty = c(2, NA),
      lwd = c(2, NA),
      bty = "n",
      cex = 0.85
    )
  },
  "rarefaction_curve",
  width = 6.5,
  height = 5.5
)

message("Rarefaction curve exported.")

# Rarefy Phyloseq Object -------------------------------------------------------

physeq_rarefied <- compute_and_save(
  file.path(dir_rds, "physeq_rarefied.rds"),
  function() {
    rarefy_even_depth(
      physeq,
      sample.size = target_depth,
      rngseed = SEED,
      replace = FALSE,
      verbose = TRUE
    )
  }
)

message(sprintf(
  "Rarefaction complete: %d samples retained at an even depth of %d reads/sample.",
  nsamples(physeq_rarefied),
  target_depth
))
