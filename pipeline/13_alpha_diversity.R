# 13. ALPHA DIVERSITY

suppressPackageStartupMessages({
  library(phyloseq)
  library(tidyverse)
})

# Estimate Richness & Merge Metadata -------------------------------------------

physeq_rarefied <- readRDS(file.path(dir_rds, "physeq_rarefied.rds"))

# Estimate alpha diversity metrics
base_measures <- setdiff(alpha_measures, "Pielou")
if ("Pielou" %in% alpha_measures) base_measures <- union(base_measures, c("Observed", "Shannon"))
alpha_div <- estimate_richness(physeq_rarefied, measures = base_measures)
if ("Pielou" %in% alpha_measures) {
  alpha_div$Pielou <- NA_real_
  valid_richness <- alpha_div$Observed > 1
  alpha_div$Pielou[valid_richness] <-
    alpha_div$Shannon[valid_richness] / log(alpha_div$Observed[valid_richness])
  # Prerequisites are calculated even if not selected, but stay out of exports.
  alpha_div <- alpha_div[, !names(alpha_div) %in% setdiff(c("Observed", "Shannon"), alpha_measures),
                         drop = FALSE]
}
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

# Decimal p-values where readable; scientific notation preserves tiny values.
format_alpha_p <- function(p) {
  vapply(p, function(value) {
    if (!is.finite(value)) return("not tested")
    paste0("p = ", if (value > 0 && value < 0.0001) {
      formatC(value, format = "e", digits = 3)
    } else sprintf("%.4f", value))
  }, character(1))
}

alpha_facet_labels <- function(metrics, p, tests) {
  paste0(metrics, "\n", ifelse(is.finite(p),
    paste(tests, format_alpha_p(p)), "not tested"))
}

analyse_alpha_metric <- function(values, groups, metric) {
  groups <- factor(groups)
  group_names <- levels(groups)
  k <- length(group_names)
  valid <- is.finite(values) & !is.na(groups)
  data <- data.frame(Value = values[valid],
                     Group = factor(as.integer(groups[valid]), levels = seq_len(k)))
  # Numeric internal labels make pair mapping safe even for names with hyphens.
  assumptions <- bind_rows(lapply(seq_len(k), function(i) {
    x <- data$Value[data$Group == i]
    reason <- if (length(x) < 3L) "Fewer than 3 finite values" else if (length(x) > 5000L) {
      "More than 5000 finite values (Shapiro-Wilk limit)"
    } else if (length(unique(x)) < 2L) "Constant values" else ""
    sw <- if (!nzchar(reason)) shapiro.test(x) else NULL
    data.frame(Metric = metric, Group = group_names[i], N = length(x),
               Test = "Shapiro-Wilk", Statistic = if (is.null(sw)) NA_real_ else unname(sw$statistic),
               p_value = if (is.null(sw)) NA_real_ else sw$p.value,
               Result = if (is.null(sw)) "not testable" else if (sw$p.value > 0.05) "passed" else "failed",
               Reason = reason)
  }))
  normality <- if (all(assumptions$Result == "passed")) "passed" else if (
    any(assumptions$Result == "failed")) "failed" else "not testable"
  variance <- "not applicable"
  variance_p <- variance_stat <- NA_real_
  if (normality == "passed" && k >= 2L) {
    deviations <- data
    deviations$Value <- abs(data$Value - ave(data$Value, data$Group, FUN = median))
    lev <- summary(aov(Value ~ Group, data = deviations))[[1]]
    variance_stat <- unname(lev$`F value`[1])
    variance_p <- unname(lev$`Pr(>F)`[1])
    variance <- if (!is.finite(variance_p)) "not testable" else if (variance_p > 0.05) "passed" else "failed"
  }
  selected <- if (normality == "passed" && variance == "passed") "ANOVA" else "Kruskal-Wallis"
  statistic <- p_value <- df1 <- df2 <- NA_real_
  status <- "not testable"
  reason <- "Fewer than two groups have finite values"
  present <- which(table(data$Group) > 0)
  fit <- NULL
  if (length(present) >= 2L) {
    if (selected == "ANOVA") {
      fit <- aov(Value ~ Group, data = data)
      result <- summary(fit)[[1]]
      statistic <- unname(result$`F value`[1])
      p_value <- unname(result$`Pr(>F)`[1])
      df1 <- result$Df[1]
      df2 <- result$Df[2]
    } else {
      result <- kruskal.test(Value ~ Group, data = data)
      statistic <- unname(result$statistic)
      p_value <- result$p.value
      df1 <- unname(result$parameter)
    }
    status <- if (is.finite(p_value)) "tested" else "not testable"
    reason <- if (is.finite(p_value)) "" else "No rank variation across finite values"
  }
  if (!is.finite(p_value)) p_value <- NA_real_
  if (!is.finite(statistic)) statistic <- NA_real_
  overall <- data.frame(
    Metric = metric, Normality = normality, Equal_variance = variance,
    Variance_test = if (normality == "passed" && k >= 2L) "Levene (median-centered)" else "not applicable",
    Variance_statistic = variance_stat, Variance_p_value = variance_p,
    Test = selected, Statistic = statistic, df = df1, Residual_df = df2,
    p_value = p_value, N = nrow(data), Groups_tested = length(present),
    Status = status, Reason = reason
  )
  pairs <- if (k >= 2L) combn(seq_len(k), 2, simplify = FALSE) else list()
  pairwise <- bind_rows(lapply(pairs, function(pair) data.frame(
    Metric = metric, Group1 = group_names[pair[1]], Group2 = group_names[pair[2]],
    Test = if (selected == "ANOVA") "Tukey HSD" else "Dunn",
    Statistic_type = if (selected == "ANOVA") "q (studentized range)" else "Z",
    Raw_p_method = if (selected == "ANOVA") "Two-sided pooled residual t" else "Two-sided normal (Dunn Z)",
    Statistic = NA_real_, p_value = NA_real_, p_adj = NA_real_,
    Adjustment = if (selected == "ANOVA") "Tukey" else "BH",
    Significance = "not tested", Status = "not testable"
  )))
  if (status == "tested" && nrow(pairwise)) {
    if (selected == "ANOVA") {
      tukey <- TukeyHSD(fit, "Group")$Group
      means <- tapply(data$Value, data$Group, mean)
      sizes <- table(data$Group)
      mse <- deviance(fit) / df.residual(fit)
      for (j in seq_along(pairs)) {
        pair <- pairs[[j]]
        t <- (means[pair[1]] - means[pair[2]]) /
          sqrt(mse * (1 / sizes[pair[1]] + 1 / sizes[pair[2]]))
        pairwise$Statistic[j] <- unname(abs(t) * sqrt(2))
        # Tukey does not supply an unadjusted p; use the equivalent pooled t.
        pairwise$p_value[j] <- 2 * pt(-abs(t), df.residual(fit))
        pairwise$p_adj[j] <- tukey[paste(rev(pair), collapse = "-"), "p adj"]
      }
    } else {
      if (!requireNamespace("dunn.test", quietly = TRUE)) {
        stop("Missing R package dunn.test. Run BIOME dependency setup before alpha statistical testing.")
      }
      # Use two-sided Dunn p-values and stats::p.adjust for monotone BH values.
      invisible(capture.output(dunn <- dunn.test::dunn.test(
        data$Value, droplevels(data$Group), method = "none", altp = TRUE,
        kw = FALSE, table = FALSE, list = FALSE, interpret = FALSE
      )))
      ids <- gsub(" ", "", dunn$comparisons, fixed = TRUE)
      raw <- dunn$altP
      adjusted <- p.adjust(raw, method = "BH")
      for (j in seq_along(pairs)) {
        pair <- pairs[[j]]
        index <- match(paste(pair, collapse = "-"), ids)
        if (is.na(index)) next
        pairwise$Statistic[j] <- dunn$Z[index]
        pairwise$p_value[j] <- raw[index]
        pairwise$p_adj[j] <- adjusted[index]
      }
    }
    pairwise$Significance <- significance_label(pairwise$p_adj)
    pairwise$Status <- ifelse(is.finite(pairwise$p_adj), "tested", "not testable")
  }
  list(overall = overall, assumptions = assumptions, pairwise = pairwise)
}

p_labels <- list()
alpha_test_caption <- NULL
alpha_comparisons <- data.frame()
alpha_overall <- alpha_assumptions <- data.frame()
if (run_alpha_statistics) {
  analyses <- lapply(alpha_measures, function(metric) {
    analyse_alpha_metric(alpha_div[[metric]], alpha_div[[group_column]], metric)
  })
  alpha_overall <- bind_rows(lapply(analyses, `[[`, "overall"))
  alpha_assumptions <- bind_rows(lapply(analyses, `[[`, "assumptions"))
  alpha_comparisons <- bind_rows(lapply(analyses, `[[`, "pairwise"))
  # Preserve the pairwise schema even when no group pair exists.
  if (!ncol(alpha_comparisons)) alpha_comparisons <- data.frame(
    Metric = character(), Group1 = character(), Group2 = character(), Test = character(),
    Statistic_type = character(), Raw_p_method = character(), Statistic = numeric(), p_value = numeric(),
    p_adj = numeric(), Adjustment = character(), Significance = character(), Status = character())
  write.csv(alpha_overall, file.path(dir_tables, "alpha_diversity_overall.csv"), row.names = FALSE)
  write.csv(alpha_assumptions, file.path(dir_tables, "alpha_diversity_assumptions.csv"), row.names = FALSE)
  write.csv(alpha_comparisons, file.path(dir_tables, "alpha_diversity_pairwise.csv"), row.names = FALSE)
  p_labels <- setNames(alpha_facet_labels(alpha_overall$Metric, alpha_overall$p_value,
                                         alpha_overall$Test), alpha_overall$Metric)
  if (n_groups %in% 2:3) {
    alpha_test_caption <- "Pairwise Tukey HSD or Dunn tests; labels use Tukey- or BH-adjusted p-values."
  }
  print(alpha_overall)
} else {
  message("Alpha statistical comparisons disabled; saving metrics and plots only.")
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

# Stack pair brackets separately in each facet, using its own data range.
alpha_annotations <- data.frame()
if (n_groups %in% 2:3 && nrow(alpha_comparisons)) {
  alpha_annotations <- alpha_comparisons %>% filter(is.finite(p_adj))
  if (nrow(alpha_annotations)) {
    alpha_annotations <- alpha_annotations %>%
      group_by(Metric) %>%
      group_modify(function(comparisons, key) {
        values <- alpha_long$Value[alpha_long$Metric == key$Metric]
        bounds <- range(values[is.finite(values)])
        span <- max(diff(bounds), abs(bounds) * 0.1, 0.1)
        comparisons %>% mutate(
          xmin = match(Group1, present_group_levels),
          xmax = match(Group2, present_group_levels)
        ) %>% arrange(xmax - xmin, xmin) %>% mutate(
          y = bounds[2] + span * (0.12 + 0.16 * (row_number() - 1)),
          tip = span * 0.035,
          label_y = y + span * 0.025
        )
      }) %>% ungroup()
    alpha_annotations$FacetLabel <- factor(
      p_labels[alpha_annotations$Metric], levels = levels(alpha_long$FacetLabel))
  }
}

# Many groups need the full figure width, rather than three narrow panels.
alpha_ncol <- if (n_groups > 4L) 1L else min(3L, length(alpha_measures))
alpha_nrow <- ceiling(length(alpha_measures) / alpha_ncol)
alpha_plot_width <- max(8, 3.5 * alpha_ncol)
alpha_plot_height <- max(6, 4.2 * alpha_nrow)
if (nrow(alpha_annotations)) {
  max_brackets <- max(table(alpha_annotations$Metric))
  alpha_plot_height <- alpha_plot_height + max(0, max_brackets - 3) * 0.18 * alpha_nrow
}

alpha_plot <- ggplot(alpha_long, aes(x = .data[[group_column]], y = Value, fill = .data[[group_column]])) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.6) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.6, size = 1.8, shape = 21, color = "black") +
  facet_wrap(~FacetLabel, scales = "free_y", ncol = alpha_ncol) +
  labs(
    x = group_column,
    y = "Diversity Metric Value",
    fill = group_column,
    caption = if (!is.null(alpha_test_caption)) {
      paste(strwrap(alpha_test_caption, width = 85), collapse = "\n")
    }
  ) +
  theme_bw(base_size = 13) +
  theme(
    strip.background = element_rect(fill = "#f8f9fa", color = "#dcdde1"),
    strip.text       = element_text(face = "bold", size = 11),
    panel.spacing    = grid::unit(1.2, "lines"),
    axis.text.x      = element_text(angle = if (n_groups > 3) 45 else 0, hjust = if (n_groups > 3) 1 else 0.5),
    plot.caption     = element_text(size = 9, hjust = 0, margin = margin(t = 12)),
    plot.caption.position = "plot",
    legend.position  = "none" # Redundant with x-axis labels
  )

if (nrow(alpha_annotations)) {
  alpha_plot <- alpha_plot +
    geom_segment(data = alpha_annotations,
                 aes(x = xmin, xend = xmax, y = y, yend = y), inherit.aes = FALSE) +
    geom_segment(data = alpha_annotations,
                 aes(x = xmin, xend = xmin, y = y, yend = y - tip), inherit.aes = FALSE) +
    geom_segment(data = alpha_annotations,
                 aes(x = xmax, xend = xmax, y = y, yend = y - tip), inherit.aes = FALSE) +
    geom_text(data = alpha_annotations,
              aes(x = (xmin + xmax) / 2, y = label_y, label = Significance),
              inherit.aes = FALSE, vjust = 0, size = 3.5) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.12)))
}

save_plot(alpha_plot, "alpha_diversity_boxplots",
          width = alpha_plot_width, height = alpha_plot_height)
