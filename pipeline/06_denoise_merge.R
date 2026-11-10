# 6. DENOISE + MERGE PAIRS

suppressPackageStartupMessages(library(dada2))

fastq_info <- readRDS(file.path(dir_rds, "fastq_info_filtered.rds"))
err_fwd    <- readRDS(file.path(dir_rds, "err_fwd.rds"))
err_rev    <- readRDS(file.path(dir_rds, "err_rev.rds"))

# Share dereplicated reads within this run ----------------------

get_derep_fwd <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      cached <- derepFastq(fastq_info$fnF.filt, verbose = FALSE)
      names(cached) <- fastq_info$SampleID
    }
    cached
  }
})

get_derep_rev <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      cached <- derepFastq(fastq_info$fnR.filt, verbose = FALSE)
      names(cached) <- fastq_info$SampleID
    }
    cached
  }
})

# Denoise ----------------------------------------------------------------------

dada_fwd <- compute_and_save(
  file.path(dir_rds, "dada_fwd.rds"),
  function() {
    dada(
      get_derep_fwd(),
      err = err_fwd,
      pool = pooling,
      multithread = N_THREADS
    )
  }
)

dada_rev <- compute_and_save(
  file.path(dir_rds, "dada_rev.rds"),
  function() {
    dada(
      get_derep_rev(),
      err = err_rev,
      pool = pooling,
      multithread = N_THREADS
    )
  }
)

# Merge paired-end reads -------------------------------------------------------

merge_pairs_by_strategy <- function(strategy) {
  
  if (strategy == "merge") {
    return(
      mergePairs(
        dada_fwd,
        get_derep_fwd(),
        dada_rev,
        get_derep_rev(),
        minOverlap = mergepairs_min_overlap,
        maxMismatch = mergepairs_max_mismatch,
        verbose = TRUE
      )
    )
  }
  
  if (strategy == "concatenate") {
    return(
      mergePairs(
        dada_fwd,
        get_derep_fwd(),
        dada_rev,
        get_derep_rev(),
        justConcatenate = TRUE,
        verbose = TRUE
      )
    )
  }
  
  stop(
    "Unknown mergepairs_strategy: '",
    strategy,
    "'. Use 'merge' or 'concatenate'."
  )
}

mergers <- compute_and_save(
  file.path(dir_rds, "mergers.rds"),
  function() merge_pairs_by_strategy(mergepairs_strategy)
)

# Build sequence table ---------------------------------------------------------

seqtab <- compute_and_save(
  file.path(dir_rds, "seqtab.rds"),
  function() makeSequenceTable(mergers)
)

message(
  "Sequence table dimensions: ",
  paste(dim(seqtab), collapse = " x ")
)