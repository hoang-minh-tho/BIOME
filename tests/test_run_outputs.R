# Run from the repository root: Rscript tests/test_run_outputs.R
script <- readLines("pipeline/01_helpers.R", warn = FALSE)
prefix <- script[seq_len(grep("^# Set plot format", script)[1] - 1L)]
expressions <- parse(text = prefix)
run_root <- tempfile("biome-output-test-")
dir.create(run_root)
previous_directory <- getwd()
tryCatch({
  setwd(run_root)
  context <- new.env(parent = baseenv())
  # The first expression loads plotting packages, unnecessary for this check.
  for (expression in as.list(expressions)[-1]) eval(expression, context)
  stopifnot(identical(sort(list.files(".")),
                      sort(c("rds", "tables", "figures", "logs"))))
  calls <- 0L
  compute <- function() { calls <<- calls + 1L; calls }
  output <- file.path("rds", "example.rds")
  stopifnot(identical(context$compute_and_save(output, compute), 1L))
  stopifnot(identical(context$compute_and_save(output, compute), 2L))
  stopifnot(identical(readRDS(output), 2L))
  stopifnot(!any(grepl("[.]hash", list.files("rds"))))
}, finally = {
  setwd(previous_directory)
  unlink(run_root, recursive = TRUE)
})
cat("Four run folders and fresh RDS computations verified.\n")
