# 3. PRIMER REMOVAL

if (primer_method == "trim_left") {
  message("Primer sequences are not specified; removing the supplied primer lengths.")
} else {
suppressPackageStartupMessages({
  library(tidyverse)
  library(Biostrings)
})

fastq_info <- readRDS(file.path(dir_rds, "fastq_info.rds"))

# Check cutadapt -------------------------------------------------------------

cutadapt_command <- if (exists("CUTADAPT_COMMAND")) CUTADAPT_COMMAND else "cutadapt"
if (!file.exists(cutadapt_command) && Sys.which(cutadapt_command) == "") {
  stop("cutadapt is not installed or not available on PATH.")
}

# Output folder --------------------------------------------------------------

dir_primer_trimmed <- file.path(dir_rds, "primer_trimmed")
dir.create(dir_primer_trimmed, showWarnings = FALSE, recursive = TRUE)

primer_removal_log <- file.path(dir_logs, "cutadapt_primer_removal.log")
rds_out <- file.path(dir_rds, "fastq_info_primer_trimmed.rds")

fastq_info <- fastq_info %>%
  mutate(
    fnF.trim = file.path(
      dir_primer_trimmed,
      paste0(SampleID, "_R1_trimmed.fastq.gz")
    ),
    fnR.trim = file.path(
      dir_primer_trimmed,
      paste0(SampleID, "_R2_trimmed.fastq.gz")
    )
  )

# Reverse-complement primers -------------------------------------------------

primer_rc <- function(seq) {
  as.character(reverseComplement(DNAString(seq)))
}

primer_fwd_rc <- primer_rc(primer_fwd)
primer_rev_rc <- primer_rc(primer_rev)

  message("Removing primer sequences from ", nrow(fastq_info), " samples...")
  
  cutadapt_threads <- if (.Platform$OS.type == "windows") {
    message("Windows detected: Cutadapt parallel processing is disabled; using one worker.")
    "1"
  } else if (isTRUE(N_THREADS)) {
    "0"
  } else {
    as.character(N_THREADS)
  }
  threads_arg <- c("-j", cutadapt_threads)
  
  for (i in seq_len(nrow(fastq_info))) {
    
    args <- c(
      threads_arg,
      "-g", primer_fwd,
      "-G", primer_rev,
      "-a", primer_rev_rc,
      "-A", primer_fwd_rc,
      "-e", as.character(cutadapt_error_rate),
      if (discard_untrimmed) "--discard-untrimmed",
      "-n", "2",
      "-o", fastq_info$fnF.trim[i],
      "-p", fastq_info$fnR.trim[i],
      fastq_info$fnF[i],
      fastq_info$fnR[i]
    )
    
    result <- system2(
      cutadapt_command,
      vapply(args, shQuote, character(1)),
      stdout = TRUE,
      stderr = TRUE
    )
    
    cat(
      c(
        paste("### Sample:", fastq_info$SampleID[i]),
        result,
        ""
      ),
      file = primer_removal_log,
      sep = "\n",
      append = TRUE
    )
    status <- attr(result, "status")
    if (!is.null(status) && status != 0L) {
      stop("Primer removal failed for sample ", fastq_info$SampleID[i],
           ". See ", primer_removal_log)
    }
    if (!all(file.exists(c(fastq_info$fnF.trim[i], fastq_info$fnR.trim[i])))) {
      stop("Primer removal did not create both output files for ", fastq_info$SampleID[i])
    }
  }
  

message("Primer removal completed.")
saveRDS(fastq_info, rds_out)
}
