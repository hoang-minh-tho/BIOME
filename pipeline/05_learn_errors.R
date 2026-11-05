# 5. LEARN ERROR RATES

suppressPackageStartupMessages(library(dada2))

fastq_info <- readRDS(file.path(dir_rds, "fastq_info_filtered.rds"))

# Error model settings ---------------------------------------------------------

err_nbases      <- 2e7   # Number of bases used for error learning
err_max_consist <- 10    # Maximum convergence rounds

# Learn forward error model ----------------------------------------------------

err_fwd <- compute_and_save(
  file.path(dir_rds, "err_fwd.rds"),
  function() {
    learnErrors(
      fastq_info$fnF.filt,
      nbases      = err_nbases,
      randomize   = TRUE,
      MAX_CONSIST = err_max_consist,
      multithread = N_THREADS,
      verbose     = FALSE
    )
  }
)

# Learn reverse error model ----------------------------------------------------

err_rev <- compute_and_save(
  file.path(dir_rds, "err_rev.rds"),
  function() {
    learnErrors(
      fastq_info$fnR.filt,
      nbases      = err_nbases,
      randomize   = TRUE,
      MAX_CONSIST = err_max_consist,
      multithread = N_THREADS,
      verbose     = FALSE
    )
  }
)

# Plot learned error models ----------------------------------------------------

save_plot(
  plotErrors(err_fwd, nominalQ = TRUE) +
    ggplot2::ggtitle("Forward error model"),
  "error_model_forward"
)

save_plot(
  plotErrors(err_rev, nominalQ = TRUE) +
    ggplot2::ggtitle("Reverse error model"),
  "error_model_reverse"
)