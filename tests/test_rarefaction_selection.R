# Run from the repository root: Rscript tests/test_rarefaction_selection.R
suppressPackageStartupMessages(library(vegan))

# Exercise the production selection and CSV export without plotting or rarefying.
expressions <- as.list(parse("pipeline/12_rarefaction.R"))
assignment_to <- function(expr, name) {
  is.call(expr) && identical(expr[[1]], as.name("<-")) &&
    identical(expr[[2]], as.name(name))
}
find_assignment <- function(name) {
  expressions[[which(vapply(expressions, assignment_to, logical(1),
                            name = name))]]
}
eval(find_assignment("rank_rarefaction_candidates"))

# The higher-yield candidate wins even when it retains fewer samples.
candidates <- data.frame(depth = c(19895, 16023), samples_kept = c(17, 18),
                         yield = c(338215, 288414))
stopifnot(rank_rarefaction_candidates(candidates)$depth[1L] == 19895)
candidates$samples_kept <- c(18, 18)
candidates$yield <- candidates$depth * candidates$samples_kept
stopifnot(rank_rarefaction_candidates(candidates)$depth[1L] == 19895)
stopifnot(rank_rarefaction_candidates(candidates[2:1, ])$depth[1L] == 19895)

run_selection <- function(counts, plateau_fraction, max_loss, expected_depth,
                          expect_fallback = FALSE) {
  context <- new.env(parent = globalenv())
  context$depths <- rowSums(counts)
  context$n_total <- nrow(counts)
  context$asv_mat <- counts
  context$rarefaction_mode <- "auto_rarefy"
  context$max_sample_loss <- max_loss
  context$plateau_samples <- plateau_fraction
  context$plateau_gain <- 0.01
  context$read_step <- 20
  context$dir_rds <- context$dir_tables <- tempfile("rarefaction-selection-")
  dir.create(context$dir_rds)
  on.exit(unlink(context$dir_rds, recursive = TRUE))
  context$compute_and_save <- function(path, compute_fun) {
    result <- compute_fun()
    saveRDS(result, path)
    result
  }
  warnings <- character()
  withCallingHandlers({
    eval(find_assignment("depth_selection"), context)
    eval(find_assignment("target_depth"), context)
    eval(find_assignment("opt_df"), context)
    export <- Filter(function(expr) {
      is.call(expr) && identical(expr[[1]], as.name("if")) &&
        identical(expr[[2]], quote(!is.null(opt_df)))
    }, expressions)
    stopifnot(length(export) == 1L)
    eval(export[[1L]], context)
  }, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  summary <- read.csv(file.path(context$dir_tables, "rarefaction_optimization.csv"))
  saved <- readRDS(file.path(context$dir_rds, "rarefaction_depth.rds"))
  stopifnot(context$target_depth == expected_depth,
            saved$target_depth == expected_depth,
            sum(summary$selected) == 1L,
            summary$depth[summary$selected] == expected_depth,
            all(summary$samples_kept >= ceiling(nrow(counts) * (1 - max_loss))),
            any(grepl("Plateau criterion not met", warnings)) == expect_fallback)
  if (!expect_fallback) stopifnot(summary$is_plateaued[summary$selected])
}

# With two abundant ASVs, all assessable curves have reached a plateau.
# Depth 150 has higher yield (450) than depth 100 (400), and is selected despite losing a sample.
counts <- cbind(c(50, 75, 100, 125), c(50, 75, 100, 125))
run_selection(counts, 0.5, 0.5, 150)
run_selection(counts, 1, 0.5, 150, expect_fallback = TRUE)
run_selection(counts, 0.5, 0, 100)

# A diverse sample prevents depth 100 from qualifying. Depth 150 qualifies;
# when the loss limit excludes it, fallback must still retain every sample.
counts <- matrix(0, nrow = 4L, ncol = 150L)
counts[c(1, 3, 4), 1] <- c(100, 200, 250)
counts[2, ] <- 1
run_selection(counts, 0.6, 0.5, 150)
run_selection(counts, 0.6, 0, 100, expect_fallback = TRUE)
cat("Rarefaction yield ranking, plateau filtering, fallback, and exports passed.\n")
