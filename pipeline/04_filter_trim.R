# 4. FILTER AND TRIM

suppressPackageStartupMessages({
  library(dada2)
  library(ShortRead)
})

if (primer_method == "cutadapt") {
  fastq_info_path <- file.path(dir_rds, "fastq_info_primer_trimmed.rds")
  fastq_info <- readRDS(fastq_info_path)
  fastq_fwd <- fastq_info$fnF.trim
  fastq_rev <- fastq_info$fnR.trim
  trim_left <- c(0L, 0L)
} else {
  fastq_info_path <- file.path(dir_rds, "fastq_info.rds")
  fastq_info <- readRDS(fastq_info_path)
  fastq_fwd <- fastq_info$fnF
  fastq_rev <- fastq_info$fnR
  trim_left <- c(trim_left_fwd, trim_left_rev)
}

# Determine truncation length --------------------------------------------------

determine_trunc_length <- function(
    fastq_files,
    trunc_qmin,
    min_length_fraction = 0.90,
    safety_trim_bp = 0,
    n_sample_files = 5
) {

  if (!length(fastq_files)) {
    stop("At least one FASTQ file is required to determine truncation length.")
  }

  sample_indices <- unique(round(seq(
    1,
    length(fastq_files),
    length.out = min(n_sample_files, length(fastq_files))
  )))
  sampled_files <- fastq_files[sample_indices]

  quality_list <- lapply(sampled_files, function(f) {
    reads <- ShortRead::readFastq(
      ShortRead::FastqSampler(f, n = 25000)
    )
    as(quality(reads), "matrix")
  })

  max_cycle <- max(vapply(quality_list, ncol, integer(1)))

  quality_list <- lapply(quality_list, function(m) {
    if (ncol(m) < max_cycle) {
      m <- cbind(
        m,
        matrix(NA_real_,
               nrow = nrow(m),
               ncol = max_cycle - ncol(m))
      )
    }
    m
  })

  quality_matrix <- do.call(rbind, quality_list)
  read_fraction <- colMeans(!is.na(quality_matrix))
  supported_cycles <- which(read_fraction >= min_length_fraction)

  if (!length(supported_cycles)) {
    stop(
      "No sequencing cycle is reached by ",
      min_length_fraction * 100,
      "% of reads."
    )
  }

  max_supported_cycle <- max(supported_cycles)
  median_quality <- apply(
    quality_matrix[, seq_len(max_supported_cycle), drop = FALSE],
    2,
    median,
    na.rm = TRUE
  )
  first_drop <- which(median_quality < trunc_qmin)[1]

  if (is.na(first_drop)) {
    trunc_length <- max_supported_cycle
    if (safety_trim_bp > 0)
      trunc_length <- max(1, trunc_length - safety_trim_bp)
  } else {
    trunc_length <- max(1, first_drop - 1)
  }

  trunc_length
}

# Automatically choose truncation lengths -------------------------------------

trunc_values <- if (disable_truncation) {
  c(0L, 0L)
} else if (automatic_truncation) {
  c(
    determine_trunc_length(fastq_fwd, trunc_qmin),
    determine_trunc_length(fastq_rev, trunc_qmin)
  )
} else {
  trunc_len
}
trunc_len_fwd <- trunc_values[1]
trunc_len_rev <- trunc_values[2]

message(
  "Selected truncLen: forward = ",
  trunc_len_fwd,
  ", reverse = ",
  trunc_len_rev
)

remaining_overlap <- (trunc_len_fwd - trim_left[1]) +
  (trunc_len_rev - trim_left[2])
if (isTRUE(asv_length_filter_enabled)) {
  minimum_mergeable_length <- asv_min_len + mergepairs_min_overlap
  if (mergepairs_strategy == "merge" &&
      remaining_overlap < minimum_mergeable_length) {
    warning(
      "The post-trim paired read length (", remaining_overlap,
      " bp) is below the minimum expected amplicon plus overlap (",
      minimum_mergeable_length, " bp); paired-end merging may fail."
    )
  }
}

writeLines(
  c(
    paste("trunc_len_fwd:", trunc_len_fwd),
    paste("trunc_len_rev:", trunc_len_rev),
    paste("trim_left_fwd:", trim_left[1]),
    paste("trim_left_rev:", trim_left[2])
  ),
  file.path(dir_logs, "trunc_len_chosen.txt")
)

# Filter and trim --------------------------------------------------------------

dir_filtered <- file.path(dir_rds, "filtered")
dir.create(dir_filtered, showWarnings = FALSE)

fastq_info <- fastq_info %>%
  mutate(
    fnF.filt = file.path(dir_filtered, paste0(SampleID, "_F_filt.fastq.gz")),
    fnR.filt = file.path(dir_filtered, paste0(SampleID, "_R_filt.fastq.gz"))
  )

filter_summary_path <- file.path(dir_rds, "filter_summary.rds")
filter_summary <- compute_and_save(
  filter_summary_path,
  function() {
    n_samples <- length(fastq_fwd)
    filter_results <- vector("list", n_samples)

    for (i in seq_len(n_samples)) {
      started <- Sys.time()
      message(sprintf(
        "[%s] Filtering sample %d/%d: %s",
        format(started, "%Y-%m-%d %H:%M:%S"),
        i,
        n_samples,
        fastq_info$SampleID[i]
      ))

      filter_results[[i]] <- dada2::filterAndTrim(
        fwd         = fastq_fwd[i],
        filt        = fastq_info$fnF.filt[i],
        rev         = fastq_rev[i],
        filt.rev    = fastq_info$fnR.filt[i],
        truncLen    = c(trunc_len_fwd, trunc_len_rev),
        trimLeft    = trim_left,
        maxEE       = max_ee,
        truncQ      = trunc_q,
        maxN        = max_n,
        minLen      = min_len,
        rm.phix     = remove_phix,
        compress    = TRUE,
        multithread = FALSE,
        verbose     = TRUE
      )

      message(sprintf(
        "[%s] Finished sample %d/%d: %s (%.1f seconds)",
        format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        i,
        n_samples,
        fastq_info$SampleID[i],
        as.numeric(difftime(Sys.time(), started, units = "secs"))
      ))
    }
    do.call(rbind, filter_results)
  }
)

rownames(filter_summary) <- fastq_info$SampleID
filtered_info_path <- file.path(dir_rds, "fastq_info_filtered.rds")
saveRDS(fastq_info, filtered_info_path)

write.csv(
  filter_summary,
  file.path(dir_tables, "filter_summary.csv")
)
