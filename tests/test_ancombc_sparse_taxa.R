# Run from the repository root: Rscript tests/test_ancombc_sparse_taxa.R
suppressPackageStartupMessages({library(phyloseq); library(ANCOMBC)})
expressions <- as.list(parse("pipeline/16_ancombc.R"))
for (expr in expressions) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      identical(expr[[2]], as.name("inspect_ancombc_taxa"))) eval(expr)
}

# Reproduce the package's dimension collapse without a lengthy model fit.
p_hat <- matrix(seq(0.1, 0.7, 0.1), nrow = 1,
                dimnames = list("Unknown", paste0("Group", 1:7)))
q_hat <- apply(p_hat, 2, function(x) p.adjust(x, method = "holm"))
q_prim <- data.frame(q_hat, check.names = FALSE)
colnames(q_prim) <- paste0("q_", colnames(q_hat))
ss_list <- rep(list(q_prim), 4)
ss_list_prim <- lapply(ss_list, function(x) x[, grepl("q_", colnames(x))])
ss_3d_prim <- array(unlist(ss_list_prim), c(dim(ss_list_prim[[1]]), length(ss_list_prim)))
error <- tryCatch(apply(ss_3d_prim, c(1, 2), mean), error = identity)
stopifnot(inherits(error, "error"), grepl("MARGIN", conditionMessage(error)))

groups <- rep(LETTERS[1:7], c(3, 2, 3, 3, 2, 2, 3))
counts <- rbind(rep(100, 18), rep(80, 18),
                t(vapply(LETTERS[1:7], function(g) as.numeric(groups == g) * 50,
                         numeric(18))))
dimnames(counts) <- list(paste0("ASV", 1:9), paste0("S", 1:18))
taxonomy <- cbind(Kingdom = "Bacteria", Genus = c("Unknown", "Unknown", LETTERS[1:7]))
rownames(taxonomy) <- rownames(counts)
phy <- phyloseq(otu_table(counts, taxa_are_rows = TRUE), tax_table(taxonomy),
                 sample_data(data.frame(Group = groups, row.names = colnames(counts))))
check <- function(data) suppressWarnings(inspect_ancombc_taxa(data, "Genus", "Group"))
one <- check(phy)
stopifnot(one$modeled_taxa == 1L,
          one$filtering$taxon[one$filtering$eligible_for_model] == "Unknown")
none <- check(prune_taxa(taxa_names(phy)[-(1:2)], phy))
stopifnot(none$modeled_taxa == 0L)
two_phy <- phy
tax_table(two_phy)[2, "Genus"] <- "Shared2"
two <- check(two_phy)
stopifnot(two$modeled_taxa == 2L)

# Check the guard's counts against the installed package's own preprocessing.
for (data in list(phy, two_phy)) {
  checked <- suppressWarnings(ANCOMBC::data_sanity_check(
    data = data, tax_level = "Genus", fix_formula = "Group", group = "Group",
    struc_zero = TRUE, verbose = FALSE))
  inspected <- check(data)
  keep <- which(rowSums(inspected$zero_ind[, -1, drop = FALSE]) == 0)
  core <- getFromNamespace(".data_core", "ANCOMBC")(
    checked$feature_table_aggregate, checked$meta_data, prv_cut = 0.1, lib_cut = 0,
    tax_keep = keep)
  stopifnot(nrow(core$feature_table) == inspected$modeled_taxa)
}

root <- tempfile("biome-ancombc-sparse-")
dir.create(root)
tryCatch({
  for (data in list(phy, prune_taxa(taxa_names(phy)[-(1:2)], phy))) {
    context <- new.env(parent = globalenv())
    context$dir_rds <- context$dir_tables <- root
    context$group_column <- "Group"
    context$reference_group <- "A"
    context$differential_method <- "ANCOM-BC"
    context$relative_abundance_ranks <- c("Phylum", "Genus")
    context$compute_and_save <- function(...) stop("Guard failed: attempted model fitting")
    context$save_plot <- function(...) stop("Guard failed: attempted statistical plotting")
    saveRDS(data, file.path(root, "physeq.rds"))
    suppressWarnings(sys.source("pipeline/16_ancombc.R", context))
    status <- read.csv(file.path(root, "differential_abundance_status.csv"))
    stopifnot(status$status == "skipped", status$modeled_taxa == check(data)$modeled_taxa,
              !file.exists(file.path(root, "ancombc_fit.rds")),
              !file.exists(file.path(root, "ancombc_results.csv")),
              file.exists(file.path(root, "ancombc_structural_zeros.csv")),
              file.exists(file.path(root, "ancombc_taxa_filtering.csv")))
  }
}, finally = unlink(root, recursive = TRUE))
cat("ANCOMBC dimension bug reproduced; zero/one-taxon guard and filtering diagnostics passed.\n")
