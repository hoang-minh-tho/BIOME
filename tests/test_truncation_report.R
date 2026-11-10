# Verify report coordinates without changing DADA2 filtering parameters.
expressions <- as.list(parse("pipeline/04_filter_trim.R"))
definition <- Filter(function(expr) {
  is.call(expr) && identical(expr[[1]], as.name("<-")) &&
    identical(expr[[2]], as.name("write_truncation_report"))
}, expressions)[[1]]
eval(definition)
path <- tempfile(fileext = ".txt")
tryCatch({
  write_truncation_report(path, "cutadapt", c(236L, 220L), c(0L, 0L))
  report <- readLines(path)
  stopifnot(any(grepl("AFTER Cutadapt", report)),
            "trunc_len_fwd: 236" %in% report,
            "Forward: 236" %in% report,
            "Reverse: 220" %in% report,
            any(grepl("undetected primer bases may remain", report)))
  write_truncation_report(path, "trim_left", c(252L, 240L), c(16L, 20L))
  report <- readLines(path)
  stopifnot(any(grepl("BEFORE leading-base removal", report)),
            "trunc_len_fwd: 252" %in% report,
            "trim_left_fwd: 16" %in% report,
            "Forward: 236" %in% report,
            "Reverse: 220" %in% report,
            any(grepl("not measured final lengths", report)))
  write_truncation_report(path, "trim_left", c(0L, 0L), c(16L, 20L))
  report <- readLines(path)
  stopifnot("Forward: variable (truncation disabled)" %in% report,
            "Reverse: variable (truncation disabled)" %in% report)
}, finally = unlink(path))
cat("Truncation reports explain primer coordinates and consistent retained-length targets.\n")
