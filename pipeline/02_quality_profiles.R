# 2. QUALITY PROFILES ----------------------------------------------------------

suppressPackageStartupMessages({
  library(dada2)
  library(tidyverse)
})

input_files <- if (exists("pipeline_input_files", inherits = TRUE)) {
  pipeline_input_files
} else {
  list.files(project_dir, full.names = TRUE, recursive = TRUE)
}

fastq_info <- compute_and_save(
  file.path(dir_rds, "fastq_info.rds"),
  function() {
    
    # Find paired-end FASTQ files -----------------------------------------------
    
    pattern_F <- "(_R1|_1)(_001)?\\.(fastq|fq)(\\.gz)?$"
    pattern_R <- "(_R2|_2)(_001)?\\.(fastq|fq)(\\.gz)?$"
    
    # Python supplies exact selected files; never scan pipeline output.
    fnF <- sort(input_files[grepl(pattern_F, basename(input_files))])
    fnR <- sort(input_files[grepl(pattern_R, basename(input_files))])
    
    raw_sample_names <- sub(pattern_F, "", basename(fnF))
    raw_reverse_names <- sub(pattern_R, "", basename(fnR))
    
    if (!length(fnF) ||
        !identical(raw_sample_names, raw_reverse_names)) {
      stop("FASTQ files not found or forward/reverse pairs do not match.")
    }
    sample_names <- unname(pipeline_sample_id_map[raw_sample_names])
    if (anyNA(sample_names) || anyDuplicated(sample_names)) {
      stop("FASTQ sample names do not have unique metadata SampleID mappings.")
    }
    
    data.frame(
      SampleID = sample_names,
      fnF = fnF,
      fnR = fnR,
      stringsAsFactors = FALSE
    )
  }
)

if (anyDuplicated(metadata$SampleID)) {
  stop("Duplicate SampleIDs found in metadata file: ", metadata_file)
}
missing <- setdiff(fastq_info$SampleID, metadata$SampleID)
if (length(missing)) {
  stop(
    "The following FASTQ samples are missing from the selected metadata file:\n",
    paste(missing, collapse = ", ")
  )
}
fastq_info <- fastq_info %>% inner_join(metadata, by = "SampleID")
message("Successfully matched ", nrow(fastq_info),
        " samples with the selected metadata file.")

save_plot(
  plotQualityProfile(fastq_info$fnF[1]) +
    ggtitle(paste("Forward:", fastq_info$SampleID[1])),
  "quality_profile_forward"
)
save_plot(
  plotQualityProfile(fastq_info$fnR[1]) +
    ggtitle(paste("Reverse:", fastq_info$SampleID[1])),
  "quality_profile_reverse"
)
