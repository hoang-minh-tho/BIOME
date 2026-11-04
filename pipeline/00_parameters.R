# 0. Parameters
# Edit only this file before running the pipeline.

# Configuration ----------------------------------------------------------------

set.seed(1234)
SEED         <- 1234
N_THREADS    <- TRUE

# Input ------------------------------------------------------------------------

# Change this to the path of your project directory containing the FASTQ files, metadata CSV, and taxonomy reference.
project_dir  <- "C:/Users/Grace_Hoang/OneDrive/Desktop/bioinfor/Samples"

# Default for standalone runs; the desktop app passes the selected CSV path.
metadata_file <- Sys.getenv(
  "BIOME_METADATA_FILE",
  unset = Sys.getenv("DATABIO16S_METADATA_FILE",
                     unset = file.path(project_dir, "metadata.csv"))
)

# Metadata column containing the groups to compare
# (e.g., "Group", "Treatment", or "Location")
group_column <- "Group"

# Reference level for differential abundance
reference_group <- "Lagoon_Amuricata"


# Taxonomy reference -----------------------------------------------------------

reference_db_path <- if (is.null(marker_profile$default_reference_relative_path)) {
  ""
} else {
  file.path(project_dir, marker_profile$default_reference_relative_path)
}

# Primers ----------------------------------------------------------------------

primer_method <- "cutadapt"
primer_fwd <- marker_profile$standalone_primers$forward
primer_rev <- marker_profile$standalone_primers$reverse
trim_left_fwd <- marker_profile$default_trim_left[[1]]
trim_left_rev <- marker_profile$default_trim_left[[2]]

# DADA2 ------------------------------------------------------------------------

trunc_qmin <- 30

mergepairs_strategy     <- "merge"
mergepairs_min_overlap  <- 12
mergepairs_max_mismatch <- 0

chimera_method <- "consensus"
max_ee         <- c(2, 2)
max_n          <- 0
trunc_q        <- 2
remove_phix    <- TRUE

# Legacy ASV-length defaults are used only when an older config has no toggle.
asv_min_len <- marker_profile$legacy_asv_length_filter$minimum
asv_max_len <- marker_profile$legacy_asv_length_filter$maximum
asv_length_filter_enabled <- marker_profile$legacy_asv_length_filter$enabled

# Taxonomic filtering ----------------------------------------------------------

remove_unassigned <- TRUE

# Rarefaction ------------------------------------------------------------------
# Choose rarefaction mode:
# 1. "min": Rarefies all samples to the minimum sequencing depth.
# 2. "custom": Rarefies all samples to a user-defined sequencing depth.
# 3. "auto_rarefy": Meets sample-loss and plateau constraints, then maximizes
#    retained samples, breaking ties by greatest sequencing depth. If no depth
#    meets the plateau criterion, uses the same ranking within the loss limit.

rarefaction_mode <- "auto_rarefy"

# If rarefaction_mode = "custom":
custom_rarefaction_depth <- 10000

# If rarefaction_mode = "auto_rarefy":
# Most users should leave these at their default values.

## Maximum sample loss (0.10 = 10%)
max_sample_loss <- 0.10

## Plateau gain threshold 
plateau_gain    <- 0.01

## Minimum samples reaching plateau
plateau_samples <- 0.90

## Plateau evaluation step size (reads)
read_step        <- 1000

# Diversity --------------------------------------------------------------------

alpha_measures <- c("Observed", "Chao1", "Shannon", "Simpson")
beta_distances <- c("bray", "wunifrac", "unifrac")
permanova_permutations <- 999

# Relative abundance -----------------------------------------------------------
# Taxa below these thresholds are grouped as "Other"

min_rank1_abundance <- 0.01
min_rank2_abundance <- 0.02
