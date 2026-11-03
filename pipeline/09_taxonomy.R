# 9. TAXONOMY ASSIGNMENT

suppressPackageStartupMessages(library(dada2))

if (!file.exists(reference_db_path)) {
  stop(
    "Taxonomy reference not found at: ",
    reference_db_path,
    ". Check the selected reference database path."
  )
}

seqtab_final <- readRDS(file.path(dir_rds, "seqtab_final.rds"))

taxonomy <- compute_and_save(
  file.path(dir_rds, "taxonomy.rds"),
  function() {
    
    assignTaxonomy(
      seqtab_final,
      refFasta    = reference_db_path,
      taxLevels   = taxonomy_ranks,
      tryRC       = TRUE,
      multithread = N_THREADS
    )
    
  }
)

asv_ids <- setNames(paste0("ASV", seq_len(ncol(seqtab_final))), colnames(seqtab_final))
taxonomy_out <- data.frame(
  ASV_ID = unname(asv_ids[rownames(taxonomy)]),
  Sequence = rownames(taxonomy),
  taxonomy, check.names = FALSE
)

write.csv(
  taxonomy_out,
  file.path(dir_tables, "taxonomy.csv"),
  row.names = FALSE
)
