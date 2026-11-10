# Confirm real DADA2 output matches manual post-primer targets for both methods.
suppressPackageStartupMessages({library(dada2); library(ShortRead)})
root <- tempfile("biome-manual-length-")
dir.create(root)
tryCatch({
  target <- paste(rep("ACGT", 25), collapse = "")
  primer <- paste(rep("TGCA", 5), collapse = "")
  for (method in c("cutadapt", "trim_left", "legacy")) {
    sequence <- if (method == "cutadapt") target else paste0(primer, target)
    input <- file.path(root, paste0(method, ".fastq"))
    output <- file.path(root, paste0(method, "_filtered.fastq.gz"))
    writeLines(unlist(lapply(1:4, function(index) {
      c(paste0("@read", index), sequence, "+", paste(rep("I", nchar(sequence)), collapse = ""))
    })), input)
    # New GUI adds primer length; old runtime configurations retain DADA2 coordinates.
    position <- if (method == "trim_left") 100L else 80L
    leading <- if (method == "cutadapt") 0L else 20L
    filterAndTrim(input, output, truncLen = position, trimLeft = leading,
                  minLen = 1, maxN = 0, maxEE = Inf, truncQ = 0,
                  rm.phix = FALSE, multithread = FALSE, verbose = FALSE)
    reads <- readFastq(output)
    expected_length <- if (method == "legacy") 60L else 80L
    stopifnot(length(reads) == 4L, all(width(sread(reads)) == expected_length),
              all(as.character(sread(reads)) == substr(target, 1L, expected_length)))
  }
}, finally = unlink(root, recursive = TRUE))
cat("Manual lengths produce matching post-primer reads; legacy DADA2 coordinates remain unchanged.\n")
