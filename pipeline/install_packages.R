args <- commandArgs(trailingOnly = TRUE)
check_only <- "--check" %in% args
options(repos = c(CRAN = "https://cloud.r-project.org"))

user_library <- Sys.getenv("R_LIBS_USER")
if (!nzchar(user_library)) {
  user_library <- file.path(
    path.expand("~"),
    "R",
    "win-library",
    paste(R.version$major, strsplit(R.version$minor, ".", fixed = TRUE)[[1]][1], sep = ".")
  )
}
dir.create(user_library, recursive = TRUE, showWarnings = FALSE)
.libPaths(unique(c(user_library, .libPaths())))

cran_pkgs <- c("tidyverse", "vegan", "Polychrome", "phangorn", "jsonlite", "rlang", "pheatmap", "dunn.test")
bioc_pkgs <- c("dada2", "phyloseq", "ShortRead", "Biostrings", "DECIPHER", "ANCOMBC",
               "DESeq2", "ALDEx2", "decontam")
missing_cran <- cran_pkgs[!vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)]
missing_bioc <- bioc_pkgs[!vapply(bioc_pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (check_only) {
  missing <- c(missing_cran, missing_bioc)
  if (length(missing)) {
    stop("Missing R packages: ", paste(missing, collapse = ", "),
         ". Choose setup in the application to install them.")
  }
  message("R_PACKAGES_OK")
  quit(save = "no", status = 0)
}

if (length(missing_cran)) {
  message("Installing CRAN packages: ", paste(missing_cran, collapse = ", "))
  install.packages(missing_cran)
}

if (length(missing_bioc)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  message("Installing Bioconductor packages: ", paste(missing_bioc, collapse = ", "))
  BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)
}

still_missing <- c(cran_pkgs, bioc_pkgs)[
  !vapply(c(cran_pkgs, bioc_pkgs), requireNamespace, logical(1), quietly = TRUE)
]
if (length(still_missing)) {
  stop("Package installation did not complete: ", paste(still_missing, collapse = ", "))
}
message("R_PACKAGES_OK")
