# 8. READ TRACKING TABLE

suppressPackageStartupMessages(library(dada2))

fastq_info      <- readRDS(file.path(dir_rds, "fastq_info_filtered.rds"))
filter_summary  <- readRDS(file.path(dir_rds, "filter_summary.rds"))
dada_fwd        <- readRDS(file.path(dir_rds, "dada_fwd.rds"))
mergers         <- readRDS(file.path(dir_rds, "mergers.rds"))
seqtab          <- readRDS(file.path(dir_rds, "seqtab.rds"))
seqtab_nochim   <- readRDS(file.path(dir_rds, "seqtab_nochim.rds"))
seqtab_final    <- readRDS(file.path(dir_rds, "seqtab_final.rds"))

get_n <- function(x) sum(getUniques(x))

read_tracking <- compute_and_save(
  file.path(dir_rds, "read_tracking.rds"),
  function() {
    
    data.frame(
      input            = filter_summary[, "reads.in"],
      filtered         = filter_summary[, "reads.out"],
      denoised         = sapply(dada_fwd, get_n),
      merged           = sapply(mergers, get_n),
      tabled           = rowSums(seqtab),
      nonchim          = rowSums(seqtab_nochim),
      length_filtered  = rowSums(seqtab_final),
      row.names = fastq_info$SampleID
    )
    
  }
)

write.csv(
  read_tracking,
  file.path(dir_tables, "read_tracking.csv")
)

print(read_tracking)