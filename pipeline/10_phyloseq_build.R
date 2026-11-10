# 10. BUILD THE PHYLOSEQ OBJECT
# Combines ASV table, taxonomy table, and metadata

suppressPackageStartupMessages({
  library(phyloseq)
  library(Biostrings)
})

seqtab_final <- readRDS(file.path(dir_rds, "seqtab_final.rds"))
asv_ids <- setNames(paste0("ASV", seq_len(ncol(seqtab_final))), colnames(seqtab_final))
taxonomy     <- readRDS(file.path(dir_rds, "taxonomy.rds"))
fastq_info   <- readRDS(file.path(dir_rds, "fastq_info.rds"))
fastq_info   <- inner_join(fastq_info, metadata, by = "SampleID")
seqtab_decontam_path <- file.path(dir_rds, "seqtab_decontam.rds")
decontam_removed_path <- file.path(dir_rds, "decontam_removed_asvs.rds")

update_decontam_read_tracking <- function(seqtab, dir_rds, dir_tables) {
  path <- file.path(dir_rds, "read_tracking.rds")
  tracking <- readRDS(path)
  counts <- rowSums(seqtab)
  if (anyDuplicated(rownames(tracking)) ||
      !setequal(rownames(tracking), names(counts))) {
    stop("Sample IDs in decontaminated counts do not match the read-tracking table.")
  }
  tracking$decontaminated <- unname(counts[rownames(tracking)])
  tracking$decontam_removed <- tracking$length_filtered - tracking$decontaminated
  saveRDS(tracking, path)
  write.csv(tracking, file.path(dir_tables, "read_tracking.csv"))
  message("Read tracking updated after decontam: ", sum(tracking$decontam_removed),
          " reads removed across ", nrow(tracking), " samples.")
  invisible(tracking)
}

update_taxonomy_read_tracking <- function(counts, keep, mitochondria, chloroplast,
                                          dir_rds, dir_tables) {
  path <- file.path(dir_rds, "read_tracking.rds")
  tracking <- readRDS(path)
  if (anyDuplicated(rownames(tracking)) ||
      !setequal(rownames(tracking), rownames(counts))) {
    stop("Sample IDs in taxonomic counts do not match the read-tracking table.")
  }
  counts <- counts[rownames(tracking), , drop = FALSE]
  # Count each removed read once, even if an ASV has both organelle labels.
  mitochondria <- mitochondria & !keep
  chloroplast <- chloroplast & !keep & !mitochondria
  other <- !keep & !(mitochondria | chloroplast)
  tracking$mitochondria_removed <- rowSums(counts[, mitochondria, drop = FALSE])
  tracking$chloroplast_removed <- rowSums(counts[, chloroplast, drop = FALSE])
  tracking$other_taxa_removed <- rowSums(counts[, other, drop = FALSE])
  tracking$taxonomy_filtered <- rowSums(counts[, keep, drop = FALSE])
  saveRDS(tracking, path)
  write.csv(tracking, file.path(dir_tables, "read_tracking.csv"))
  message("Read tracking updated after taxonomic filtering: ",
          sum(tracking$mitochondria_removed), " mitochondrial reads, ",
          sum(tracking$chloroplast_removed), " chloroplast reads, and ",
          sum(tracking$other_taxa_removed), " other reads removed.")
  invisible(tracking)
}

physeq <- compute_and_save(
  file.path(dir_rds, "physeq.rds"),
  function() {
    is_neg <- is_neg_control[rownames(seqtab_final)]
    if (anyNA(is_neg)) {
      stop("Negative-control classification is missing for samples in seqtab_final.")
    }
    if (run_decontam) {
      if (!any(is_neg) || all(is_neg)) {
        stop("decontam requires at least one Negative Control and one non-negative sample.")
      }
      message(">> Identifying contaminant ASVs from negative controls...")
      contam_df <- decontam::isContaminant(
        seqtab_final,
        neg = is_neg,
        method = "prevalence"
      )
      contam_asvs <- rownames(contam_df)[
        !is.na(contam_df$contaminant) & contam_df$contaminant
      ]
      if (length(contam_asvs) >= ncol(seqtab_final)) {
        stop("decontam marked every ASV as a contaminant; review control classification and data.")
      }
      message("decontam identified ", length(contam_asvs), " contaminant ASVs.")
      saveRDS(
        contam_df[contam_df$contaminant %in% TRUE, , drop = FALSE],
        decontam_removed_path
      )
      seqtab_final <- seqtab_final[, !colnames(seqtab_final) %in% contam_asvs, drop = FALSE]
      taxonomy <- taxonomy[!rownames(taxonomy) %in% contam_asvs, , drop = FALSE]
      saveRDS(seqtab_final, seqtab_decontam_path)
      update_decontam_read_tracking(seqtab_final, dir_rds, dir_tables)
    }
    
# Build phyloseq object -----------------------------------------------
    
    sample_data_df <- fastq_info
    rownames(sample_data_df) <- sample_data_df$SampleID
    
    physeq_raw <- phyloseq(
      otu_table(seqtab_final, taxa_are_rows = FALSE),
      tax_table(taxonomy),
      sample_data(sample_data_df)
    )
    
    ## Preserve original ASV sequences
    asv_sequences <- DNAStringSet(taxa_names(physeq_raw))
    names(asv_sequences) <- unname(asv_ids[taxa_names(physeq_raw)])
    
    taxa_names(physeq_raw) <- names(asv_sequences)
    physeq_raw <- merge_phyloseq(physeq_raw, asv_sequences)
    
# Taxonomic filtering -------------------------------------------------
    
    tax_df <- as.data.frame(tax_table(physeq_raw))
    if (!(target_domain_rank %in% colnames(tax_df))) {
      stop(
        "The selected taxonomy profile expects rank '", target_domain_rank,
        "', but the reference produced: ", paste(colnames(tax_df), collapse = ", ")
      )
    }
    if (length(contaminant_taxa) && !(contaminant_rank %in% colnames(tax_df))) {
      stop(
        "The selected taxonomy profile expects contaminant rank '",
        contaminant_rank, "', but that rank is missing from the taxonomy table."
      )
    }

    keep <- rep(TRUE, nrow(tax_df))
    
    if (!is.null(target_domain)) {
      keep <- keep & tax_df[[target_domain_rank]] %in% target_domain
    }
    
    if (remove_unassigned) {
      keep <- keep & !is.na(tax_df[[target_domain_rank]])
    }
    
    # Organelle labels can occur at different ranks in taxonomy references.
    is_organelle <- rep(FALSE, nrow(tax_df))
    is_mitochondria <- is_chloroplast <- is_organelle
    for (rank in colnames(tax_df)) {
      labels <- tolower(trimws(as.character(tax_df[[rank]])))
      excluded <- labels %in% tolower(organelle_taxa)
      is_organelle <- is_organelle | excluded
      is_mitochondria <- is_mitochondria | (excluded & labels %in% "mitochondria")
      is_chloroplast <- is_chloroplast | (excluded & labels %in% "chloroplast")
    }
    keep <- keep & !is_organelle

    if (length(contaminant_taxa)) {
      message("Explicit genus exclusions: ", paste(contaminant_taxa, collapse = ", "))
      keep <- keep & !(tax_df[[contaminant_rank]] %in% contaminant_taxa)
    }
    
    counts <- as(otu_table(physeq_raw), "matrix")
    if (taxa_are_rows(physeq_raw)) counts <- t(counts)
    counts <- counts[, rownames(tax_df), drop = FALSE]
    update_taxonomy_read_tracking(
      counts, keep, is_mitochondria, is_chloroplast, dir_rds, dir_tables
    )
    physeq <- prune_taxa(keep, physeq_raw)
    
    message(
      sprintf(
        "Taxonomic filtering retained %d of %d taxa.",
        ntaxa(physeq),
        ntaxa(physeq_raw)
      )
    )

    if (exclude_controls || exclude_mock_controls) {
      keep_samples <- if (exclude_controls) {
        !(is_neg_control | is_pos_control)
      } else {
        !is_pos_control
      }
      names(keep_samples) <- names(is_neg_control)
      keep_samples <- keep_samples[sample_names(physeq)]
      if (anyNA(keep_samples)) {
        stop("Control classification is missing for samples in the phyloseq object.")
      }
      physeq <- prune_samples(keep_samples, physeq)
      message("Excluded ", sum(!keep_samples), " mock/community control samples before downstream analyses.")
    }
    physeq
  }
)

seqtab_for_fasta <- if (run_decontam) {
  if (!file.exists(seqtab_decontam_path)) {
    stop("Decontaminated sequence table is missing: ", seqtab_decontam_path)
  }
  readRDS(seqtab_decontam_path)
} else {
  seqtab_final
}
asv_sequences <- DNAStringSet(dada2::getSequences(seqtab_for_fasta))
names(asv_sequences) <- unname(asv_ids[dada2::getSequences(seqtab_for_fasta)])
writeXStringSet(asv_sequences, file.path(dir_tables, "asv_sequences.fasta"))

if (run_decontam) {
  if (!file.exists(decontam_removed_path)) {
    stop("Decontam results are missing: ", decontam_removed_path)
  }
  write.csv(
    readRDS(decontam_removed_path),
    file.path(dir_tables, "decontam_removed_asvs.csv"),
    row.names = TRUE
  )
}
