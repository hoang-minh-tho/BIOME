# Run from the repository root: Rscript tests/test_alpha_posthoc.R
suppressPackageStartupMessages({library(phyloseq); library(tidyverse)})
for (expr in as.list(parse("pipeline/01_helpers.R"))) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      identical(expr[[2]], as.name("significance_label"))) eval(expr)
}
for (expr in as.list(parse("pipeline/13_alpha_diversity.R"))) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) && is.symbol(expr[[2]]) &&
      as.character(expr[[2]]) %in% c("analyse_alpha_metric", "format_alpha_p", "alpha_facet_labels")) eval(expr)
}
stopifnot(identical(significance_label(c(0.0001, 0.001, 0.01, 0.05, 0.05001, NA)),
                    c("****", "***", "**", "*", "ns", "not tested")))
# Deterministic normal groups with identical variance select ANOVA.
z <- qnorm(ppoints(15))
groups <- factor(rep(c("A-B", "C", "D-E"), each = 15))
values <- c(z, z + 1, z + 2)
normal <- analyse_alpha_metric(values, groups, "Shannon")
stopifnot(normal$overall$Test == "ANOVA", normal$overall$Normality == "passed",
          normal$overall$Equal_variance == "passed", nrow(normal$assumptions) == 3,
          all(normal$pairwise$Test == "Tukey HSD"), all(normal$pairwise$Adjustment == "Tukey"))
reference <- aov(values ~ groups)
stopifnot(isTRUE(all.equal(normal$overall$p_value, summary(reference)[[1]]$`Pr(>F)`[1])))
tukey <- TukeyHSD(reference)$groups
stopifnot(all(normal$pairwise$Statistic_type == "q (studentized range)"),
          isTRUE(all.equal(normal$pairwise$p_adj,
                            ptukey(normal$pairwise$Statistic, 3, df.residual(reference), lower.tail = FALSE))))
for (j in seq_len(nrow(normal$pairwise))) {
  row <- normal$pairwise[j, ]
  stopifnot(isTRUE(all.equal(row$p_adj, unname(tukey[paste(row$Group2, row$Group1, sep = "-"), "p adj"]))))
}
raw <- pairwise.t.test(values, groups, pool.sd = TRUE, p.adjust.method = "none")$p.value
for (j in seq_len(nrow(normal$pairwise))) {
  row <- normal$pairwise[j, ]
  stopifnot(isTRUE(all.equal(row$p_value, unname(raw[row$Group2, row$Group1]))))
}
# Normal but heteroscedastic groups must not use ordinary ANOVA/Tukey.
unequal <- analyse_alpha_metric(c(z, 20*z, 100*z), groups, "Shannon")
stopifnot(unequal$overall$Normality == "passed", unequal$overall$Equal_variance == "failed",
          unequal$overall$Test == "Kruskal-Wallis")
# Check two-sided Dunn Z and BH against the independent tied-rank formula.
x <- c(rep(0, 10), 1:5, rep(0, 7), 1:8, rep(1, 10), 2:6)
ranked <- analyse_alpha_metric(x, groups, "Observed")
stopifnot(ranked$overall$Test == "Kruskal-Wallis", ranked$overall$Normality == "failed",
          ranked$overall$Equal_variance == "not applicable")
ranks <- rank(x)
n <- length(x)
ties <- table(x)
variance <- n*(n+1)/12 - sum(ties^3 - ties)/(12*(n-1))
means <- tapply(ranks, groups, mean)
sizes <- table(groups)
pairs <- combn(1:3, 2, simplify = FALSE)
expected_z <- vapply(pairs, function(p) (means[p[1]] - means[p[2]]) /
                       sqrt(variance * (1/sizes[p[1]] + 1/sizes[p[2]])), numeric(1))
expected_p <- 2*pnorm(-abs(expected_z))
stopifnot(isTRUE(all.equal(ranked$pairwise$Statistic, unname(expected_z))),
          isTRUE(all.equal(ranked$pairwise$p_value, unname(expected_p))),
          isTRUE(all.equal(ranked$pairwise$p_adj, unname(p.adjust(expected_p, "BH")))))
# A nonsignificant overall result still exports all pairs and their ns labels.
null <- analyse_alpha_metric(rep(z, 3), groups, "Shannon")
stopifnot(null$overall$p_value > 0.05, all(null$pairwise$Significance == "ns"))
small <- analyse_alpha_metric(c(1, 2, 3, 7), factor(c("A", "A", "B", "B")), "Observed")
stopifnot(small$overall$Test == "Kruskal-Wallis", small$overall$Normality == "not testable")
missing <- analyse_alpha_metric(c(NA, Inf, 1, 2), factor(c("A", "A", "B", "B")), "Pielou")
stopifnot(missing$overall$Status == "not testable", all(is.na(missing$pairwise$p_adj)))
partial <- analyse_alpha_metric(c(NA, Inf, 1, 2, 4, 5),
                                factor(rep(c("A", "B", "C"), each = 2)), "Pielou")
stopifnot(partial$overall$Groups_tested == 2,
          is.finite(partial$pairwise$p_adj[partial$pairwise$Group1 == "B"]),
          all(is.na(partial$pairwise$p_adj[partial$pairwise$Group1 == "A"])))
constant <- analyse_alpha_metric(rep(1, 45), groups, "Observed")
stopifnot(is.na(constant$overall$p_value), all(constant$pairwise$Significance == "not tested"))
# Non-testable Shapiro sample sizes are explicit and do not subsample randomly.
large <- analyse_alpha_metric(c(seq_len(5001), 1:3), factor(c(rep("A", 5001), rep("B", 3))), "Observed")
stopifnot(large$assumptions$Result[1] == "not testable", large$overall$Test == "Kruskal-Wallis")

root <- tempfile("biome-alpha-posthoc-")
dir.create(root)
tryCatch({
  context <- new.env(parent = globalenv())
  for (folder in c("rds", "tables", "figures")) {
    path <- file.path(root, folder)
    dir.create(path)
    assign(paste0("dir_", folder), path, context)
  }
  context$group_column <- "Group"
  context$alpha_measures <- c("Shannon", "Observed")
  context$run_alpha_statistics <- TRUE
  context$save_plot <- function(plot, name, width, height) {
    ggplot_build(plot)
    preview <- Sys.getenv("BIOME_ALPHA_PREVIEW_DIR")
    if (nzchar(preview)) {
      dir.create(preview, recursive = TRUE, showWarnings = FALSE)
      ggsave(file.path(preview, paste0("alpha_", context$n_groups, "_groups.png")),
             plot, width = width, height = height, dpi = 120)
    }
  }
  baseline <- rbind(c(90, 5, 3, 2), c(85, 7, 5, 3), c(80, 10, 6, 4))
  for (k in c(2L, 3L, 4L, 7L, 12L)) {
    counts <- baseline[rep(1:3, k * 3), ]
    dimnames(counts) <- list(paste0("S", seq_len(nrow(counts))), paste0("ASV", 1:4))
    metadata <- data.frame(Group = rep(paste0("G", seq_len(k)), each = 9), row.names = rownames(counts))
    saveRDS(phyloseq(otu_table(counts, taxa_are_rows = FALSE), sample_data(metadata)),
            file.path(context$dir_rds, "physeq_rarefied.rds"))
    suppressWarnings(sys.source("pipeline/13_alpha_diversity.R", context))
    stopifnot(nrow(context$alpha_comparisons) == 2*choose(k, 2),
              all(is.finite(context$alpha_comparisons$p_adj[context$alpha_comparisons$Metric == "Shannon"])),
              !any(grepl("Shapiro|Levene|\\*|\\bns\\b", unlist(context$p_labels))))
    if (k <= 3) stopifnot(nrow(context$alpha_annotations) == choose(k, 2),
                          all(context$alpha_annotations$Significance == "ns"))
    if (k >= 4) stopifnot(nrow(context$alpha_annotations) == 0,
                          is.null(context$alpha_plot$labels$caption),
                          length(context$alpha_plot$layers) == 2)
    for (name in c("overall", "assumptions", "pairwise")) {
      exported <- read.csv(file.path(context$dir_tables, paste0("alpha_diversity_", name, ".csv")))
      stopifnot(nrow(exported) > 0)
    }
  }
  context$run_alpha_statistics <- FALSE
  suppressWarnings(sys.source("pipeline/13_alpha_diversity.R", context))
  stopifnot(nrow(context$alpha_annotations) == 0, length(context$p_labels) == 0,
            is.null(context$alpha_test_caption))
}, finally = unlink(root, recursive = TRUE))
cat("Automatic assumptions, ANOVA/Tukey, tied Dunn/BH, full exports and 2/3/4/7/12-group plots passed.\n")
