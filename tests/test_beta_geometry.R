suppressPackageStartupMessages({library(phyloseq); library(vegan); library(tidyverse)})
for (expr in as.list(parse("pipeline/14_beta_diversity.R"))) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) && is.symbol(expr[[2]]) &&
      as.character(expr[[2]]) %in% c("inspect_distance_geometry", "prepare_beta_distance")) eval(expr)
}
euclidean <- dist(matrix(c(0, 0, 1, 0, 0, 1, 1, 1), ncol = 2, byrow = TRUE))
unchanged <- prepare_beta_distance(euclidean, "Euclidean")
stopifnot(!unchanged$sqrt_dist, identical(unchanged$distance, euclidean))
# A square with Manhattan distances has a negative PCoA eigenvalue;
# its square root is Euclidean. Include labels and non-unit scaling.
non_euclidean <- dist(matrix(c(0, 0, 1, 0, 0, 1, 1, 1), ncol = 2, byrow = TRUE), method = "manhattan")
attr(non_euclidean, "Labels") <- paste0("S", 1:4)
for (scale in c(1e-6, 1, 1e6)) {
  distance <- non_euclidean * scale
  transformed <- prepare_beta_distance(distance, "Manhattan")
  stopifnot(transformed$sqrt_dist,
            identical(as.numeric(transformed$distance), sqrt(as.numeric(distance))),
            identical(attr(transformed$distance, "Labels"), attr(distance, "Labels")),
            transformed$diagnostics$raw_negative_eigenvalues > 0,
            transformed$diagnostics$analysis_negative_eigenvalues == 0)
}
stopifnot(inspect_distance_geometry(dist(matrix(0, nrow = 4, ncol = 2)))$negative == 0)
stopifnot(inspect_distance_geometry(euclidean)$negative == 0)
rounding_noise <- as.dist(matrix(c(0, 1, 2 + 1e-10, 1, 0, 1, 2 + 1e-10, 1, 0), nrow = 3))
noise_check <- inspect_distance_geometry(rounding_noise)
stopifnot(min(noise_check$eigenvalues) < 0, noise_check$negative == 0)
# Square root is not guaranteed to repair arbitrary dissimilarities.
residual <- as.dist(matrix(c(0, 1, 1, 1, 0, 10, 1, 10, 0), nrow = 3))
warnings <- character()
withCallingHandlers(remaining <- prepare_beta_distance(residual, "Residual"), warning = function(w) {
  warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
})
stopifnot(remaining$sqrt_dist, remaining$diagnostics$analysis_negative_eigenvalues > 0,
          any(grepl("not fully Euclidean", warnings)))
# Check the actual pipeline consumers use the same distance geometry once.
root <- tempfile("biome-geometry-")
dir.create(root)
tryCatch({
  for (use_sqrt in c(FALSE, TRUE)) {
    context <- new.env(parent = globalenv())
    for (folder in c("rds", "tables", "figures")) {
      path <- file.path(root, as.character(use_sqrt), folder)
      dir.create(path, recursive = TRUE)
      assign(paste0("dir_", folder), path, context)
    }
    coordinates <- matrix(c(0,0, 1,0, 0,1, 1,1, 2,0, 0,2), ncol = 2, byrow = TRUE)
    rownames(coordinates) <- paste0("S", 1:6)
    distance <- dist(coordinates, method = if (use_sqrt) "manhattan" else "euclidean")
    counts <- matrix(10, nrow = 6, ncol = 3,
                     dimnames = list(rownames(coordinates), paste0("ASV", 1:3)))
    metadata <- data.frame(Group = rep(c("A", "B"), each = 3), row.names = rownames(coordinates))
    physeq <- phyloseq(otu_table(counts, taxa_are_rows = FALSE), sample_data(metadata))
    saveRDS(physeq, file.path(context$dir_rds, "physeq_rarefied.rds"))
    context$group_column <- "Group"
    context$beta_distances <- "bray"
    context$ordination_method <- "PCoA"
    context$run_permanova <- TRUE
    context$run_betadisper <- TRUE
    context$permanova_permutations <- 9L
    context$compute_and_save <- function(path, fun) { saveRDS(distance, path); distance }
    context$ordinate <- function(object, method, distance, ...) {
      context$ordination_distance <- distance
      phyloseq::ordinate(object, method = method, distance = distance, ...)
    }
    context$save_plot <- function(...) invisible(NULL)
    context$betadisper <- function(d, group, sqrt.dist) {
      context$dispersion_flag <- sqrt.dist
      vegan::betadisper(d, group, sqrt.dist = sqrt.dist)
    }
    context$adonis2 <- function(formula, data, permutations, sqrt.dist) {
      context$permanova_flag <- sqrt.dist
      dist_matrix <- context$dist_matrix
      vegan::adonis2(formula, data = data, permutations = permutations, sqrt.dist = sqrt.dist)
    }
    sys.source("pipeline/14_beta_diversity.R", context)
    diagnostics <- read.csv(file.path(context$dir_tables, "beta_distance_geometry.csv"))
    stopifnot(diagnostics$sqrt_dist == use_sqrt,
              identical(context$dispersion_flag, use_sqrt),
              identical(context$permanova_flag, use_sqrt),
              identical(as.numeric(context$ordination_distance),
                        as.numeric(if (use_sqrt) sqrt(distance) else distance)),
              identical(readRDS(file.path(context$dir_rds, "dist_bray.rds")), distance),
              identical(readRDS(file.path(context$dir_rds, "dist_bray_analysis.rds")),
                        context$ordination_distance))
    tests <- read.csv(file.path(context$dir_tables, "permanova_results.csv"))
    stopifnot(all(tests$sqrt_dist == use_sqrt))
  }
}, finally = unlink(root, recursive = TRUE))
cat("Geometry detection, numerical tolerance, residual warnings, and consistent single transformation passed.\n")
