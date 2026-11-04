# Exercise the computation callback without requiring DADA2 or sequencing data.
# Run from the repository root: Rscript tests/test_chimera_threads.R
expressions <- parse("pipeline/07_chimera_length_filter.R")
assignment <- Filter(function(expr) {
  is.call(expr) && identical(expr[[1]], as.name("<-")) &&
    identical(expr[[2]], as.name("seqtab_nochim"))
}, as.list(expressions))[[1]]
callback_expression <- assignment[[3]][[3]]

for (os_type in c("windows", "unix")) {
  for (method in c("consensus", "pooled", "none")) {
    context <- new.env(parent = baseenv())
    context$.Platform <- list(OS.type = os_type)
    context$N_THREADS <- 7L
    context$chimera_method <- method
    context$seqtab <- matrix(1L, nrow = 1, ncol = 1)
    context$removeBimeraDenovo <- function(seqtab, method, multithread, verbose) {
      multithread
    }
    callback <- eval(callback_expression, envir = context)
    expected <- if (method == "none") {
      context$seqtab
    } else if (os_type == "windows" && method == "pooled") {
      FALSE
    } else {
      7L
    }
    stopifnot(identical(callback(), expected))
  }
}
cat("Chimera worker selection passed for Windows and Unix methods.\n")
