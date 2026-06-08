# -----------------------------------------------------------------------------#
# STEP 2: Load 10x Genomics data
# -----------------------------------------------------------------------------#

# Define paths
sample_name <- "Healthy_1"
cellranger_output <- "/nfs/jbailey5/baileyweb/colabs/kwamae/software/scrna-seq-analysis/cellranger-matrix-counts"

# --- THE CRITICAL MATRIX DECISION: RAW VS FILTERED ---
# Cell Ranger (10x software) outputs two entirely distinct count matrices:
#
# 1. RAW (`raw_feature_bc_matrix`):
#    Includes EVERY single droplet that passed through the microfluidic chip
#    (often 100,000+ barcodes). This includes true cells, empty background
#    droplets, ambient "soup" RNA, and doublets.
#
# 2. FILTERED (`filtered_feature_bc_matrix`):
#    Includes only the barcodes that Cell Ranger's internal, conservative
#    algorithm guessed were actual cells (usually ~3,000 - 10,000 cells).
#
# WHY THIS TUTORIAL USES OPTION A (RAW):
# Cell Ranger's built-in filtering often accidentally throws away real, ultra-
# small cells (like resting lymphocytes) or retains large dead cell debris.
# Loading the RAW matrix allows you to run modern, superior downstream tools
# like EmptyDrops (Step 4), SoupX (Step 5), and scDblFinder (Step 6) to perform
# customized, high-precision purification that outperforms Cell Ranger.

# OPTION A: Load RAW matrix (what this tutorial demonstrates)
counts <- Read10X(
  data.dir = file.path(cellranger_output, "raw_feature_bc_matrix")
)

# # OPTION B: Load FILTERED matrix (alternative)
# counts <- Read10X(
#   data.dir = file.path(cellranger_output, "filtered_feature_bc_matrix")
# )


# --- Create Seurat Object ---
# We initialize the object with ZERO pre-filtering constraints.
# Setting `min.cells = 0` and `min.features = 0` is required because we want to
# keep the ambient background noise intact for now. If we filtered genes or
# cells here, tools like SoupX and EmptyDrops would lack the background
# baseline data they need to estimate contamination accurately.

seurat_obj <- CreateSeuratObject(
  counts = counts,
  project = sample_name,
  min.cells = 0, # Explicitly do NOT filter genes yet
  min.features = 0 # Explicitly do NOT filter cells yet
)


# --- Add Experimental Metadata ---
# Annotates the dataset with project-specific structural metadata. This is
# essential for downstream analyses when merging multiple patient batches or
# performing multi-condition differential expression (e.g., Healthy vs Disease).

# Using base R
# seurat_obj$sample_id <- "Healthy_1"
# seurat_obj$sra_id <- "SRR14575500"
# seurat_obj$condition <- "Healthy"
# seurat_obj$patient_id <- "Donor_1"
# seurat_obj$time_point <- NA

# PERFORMANCE NOTE:
# Mutating the `@meta.data` dataframe directly using dplyr is significantly
# faster than assigning individual columns one-by-one using the `$` operator,
# as it executes all annotations in a single, unified memory block update.

seurat_obj@meta.data <- seurat_obj@meta.data %>%
  mutate(
    sample_id  = "Healthy_1",
    sra_id     = "SRR14575500",
    condition  = "Healthy",
    patient_id = "Donor_1",
    time_point = NA_character_ # Explicit typed NA for column consistency
  )


# --- Dimensions Quick Check ---
cat("Loaded:", ncol(seurat_obj), "droplets ×", nrow(seurat_obj), "genes\n")
cat("(Most are empty - EmptyDrops will filter in Step 4)\n")
