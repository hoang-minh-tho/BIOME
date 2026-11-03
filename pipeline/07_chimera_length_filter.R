# 7. CHIMERA REMOVAL + ASV LENGTH FILTER

suppressPackageStartupMessages(library(dada2))

seqtab <- readRDS(file.path(dir_rds, "seqtab.rds"))

# Remove chimeras --------------------------------------------------------------

seqtab_nochim <- compute_and_save(
  file.path(dir_rds, "seqtab_nochim.rds"),
  function() {
    if (identical(chimera_method, "none")) {
      message("Chimera removal skipped.")
      return(seqtab)
    }
    # Pooled detection uses mclapply, which cannot fork on Windows.
    chimera_threads <- if (.Platform$OS.type == "windows" &&
                           identical(chimera_method, "pooled")) {
      message("Windows detected: pooled chimera removal uses one worker.")
      FALSE
    } else {
      N_THREADS
    }
    removeBimeraDenovo(
      seqtab,
      method = chimera_method,
      multithread = chimera_threads,
      verbose = TRUE
    )
  }
)

message(
  "Fraction of reads kept after chimera removal: ",
  round(sum(seqtab_nochim) / sum(seqtab), 4)
)

# Optionally filter ASVs by length --------------------------------------------

asv_lengths <- nchar(getSequences(seqtab_nochim))

write.csv(
  table(asv_lengths),
  file.path(dir_tables, "asv_length_distribution.csv"),
  row.names = FALSE
)

seqtab_final <- compute_and_save(
  file.path(dir_rds, "seqtab_final.rds"),
  function() {
    if (!isTRUE(asv_length_filter_enabled)) {
      message("ASV length filtering is disabled; retaining all non-chimeric ASVs.")
      return(seqtab_nochim)
    }

    keep_by_length <- asv_lengths >= asv_min_len &
      asv_lengths <= asv_max_len
    
    message(
      sum(keep_by_length), " of ",
      length(keep_by_length),
      " ASVs kept within the ",
      asv_min_len, "-", asv_max_len, " bp window"
    )
    
    seqtab_nochim[, keep_by_length, drop = FALSE]
    
  }
)
