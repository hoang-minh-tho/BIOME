# 11. PHYLOGENETIC TREE (for UniFrac)

suppressPackageStartupMessages({
  library(phyloseq)
  library(DECIPHER)
  library(phangorn)
})

physeq <- readRDS(file.path(dir_rds, "physeq.rds"))

tree_rearrangement <- "none"

seqs <- refseq(physeq)
seqs <- seqs[taxa_names(physeq)]

message("Building tree from ", length(seqs), " ASVs")

# Multiple sequence alignment --------------------------------------------------

alignment <- compute_and_save(
  file.path(dir_rds, "tree_alignment.rds"),
  function() {
    
    message("[tree 1/5] aligning ASV sequences...")
    AlignSeqs(
      seqs,
      anchor = NA,
      processors = NULL,
      verbose = FALSE
    )
    
  }
)

message("[tree 1/5] alignment done")

# Convert alignment ------------------------------------------------------------

phang_align <- compute_and_save(
  file.path(dir_rds, "tree_phang_align.rds"),
  function() {
    
    message("[tree 2/5] converting alignment...")
    phyDat(as(alignment, "matrix"), type = "DNA")
    
  }
)

dist_ml <- compute_and_save(
  file.path(dir_rds, "tree_dist_ml.rds"),
  function() {
    
    message("[tree 3/5] computing distance matrix...")
    dist.ml(phang_align)
    
  }
)

message("[tree 3/5] distance matrix done")

# Starting tree ----------------------------------------------------------------

tree_nj <- compute_and_save(
  file.path(dir_rds, "tree_nj.rds"),
  function() {
    
    message("[tree 4/5] building neighbour-joining tree...")
    NJ(dist_ml)
    
  }
)

message("[tree 4/5] neighbour-joining tree done")

# Maximum-likelihood tree ------------------------------------------------------

tree <- compute_and_save(
  file.path(dir_rds, "phy_tree.rds"),
  function() {
    
    message("[tree 5/5] fitting GTR tree...")
    
    fit <- pml(tree_nj, data = phang_align)
    
    fit_gtr <- update(fit, k = 4, inv = 0.2)
    
    fit_gtr <- optim.pml(
      fit_gtr,
      model = "GTR",
      optInv = TRUE,
      optGamma = TRUE,
      rearrangement = tree_rearrangement,
      control = pml.control(trace = 0)
    )
    
    fit_gtr$tree
    
  }
)

message("[tree 5/5] tree fitting done")

# Attach tree to phyloseq ------------------------------------------------------

physeq <- compute_and_save(
  file.path(dir_rds, "physeq.rds"),
  function() {
    merge_phyloseq(physeq, phy_tree(tree))
  }
)