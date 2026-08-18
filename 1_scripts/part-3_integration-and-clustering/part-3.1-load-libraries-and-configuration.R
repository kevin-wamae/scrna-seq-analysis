# ****************************************************************************#
# STEP 1: Load required libraries
# ****************************************************************************#


# Declare pipeline variables
# ----------------------------------------------------------------------------#
# The date-stamped job run ID (from HPC scheduler) for tracking
RUN_ID <- "2026_06_09_brown_job_3058993" 


# Terminal output colouring (load first so errors/warnings are coloured)
library(colorout)


# Core single-cell analysis
library(Seurat)
library(SeuratObject)


# Integration methods
library(harmony)
library(batchelor)
library(SeuratWrappers)


# Visualization and data manipulation
library(ggplot2)
library(ggrepel)
library(dplyr)
library(patchwork)
library(RColorBrewer)
library(viridis)
library(reshape2)


# Quality metrics and utilities
library(FNN)              # K-nearest neighbor calculations for mixing metrics
library(cluster)          # Silhouette scores for clustering quality


# Parallel processing
library(future)
library(furrr)
options(future.globals.maxSize = 20 * 1024^3)  # Increase to 20GB for large datasets

# # Set working directory (adjust to your path)
# setwd("~/GSE174609_scRNA/integration_analysis")

# --- DYNAMIC OUTPUT PATH REGISTRATION ---
# ****************************************************************************#
# Dynamically construct paths inside a date-stamped job run ID.
PLOTS_OUT_DIR <- file.path(
  "3_output", RUN_ID, "integration_and_clustering", "plots"
)
PLOTS_COMPARISON_DIR <- file.path(PLOTS_OUT_DIR, "integration_comparison")
PLOTS_CLUSTERING_DIR <- file.path(PLOTS_OUT_DIR, "clustering")
DATA_OUT_DIR <- file.path(
  "3_output", RUN_ID, "integration_and_clustering", "integrated_data"
)
METADATA_OUT_DIR <- file.path(
  "3_output", RUN_ID, "integration_and_clustering", "metadata"
)


# --- AUTOMATED DIRECTORY PROVISIONING ---
# ****************************************************************************#
# Construct the folder trees on disk:
#  - `recursive = TRUE`: Instructs R to build missing parent paths on the fly
#  - `showWarnings = FALSE`: Prevents the script from halting or throwing noise
#    if the directories were already initialized by an earlier script step

dir.create(PLOTS_OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOTS_COMPARISON_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOTS_CLUSTERING_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(METADATA_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("✓ Output directories securely synchronized to disk:\n")
cat("  • Visual Diagnostic Canvas:", PLOTS_OUT_DIR, "\n")
cat("  •   ├─ Integration Comparison:", PLOTS_COMPARISON_DIR, "\n")
cat("  •   └─ Clustering:", PLOTS_CLUSTERING_DIR, "\n")
cat("  • Integrated Data Archive:", DATA_OUT_DIR, "\n")
cat("  • Metadata Ledger:", METADATA_OUT_DIR, "\n\n")

# Set random seed for reproducibility
set.seed(42)

# Configure plotting defaults
theme_set(theme_classic(base_size = 12))

cat("Environment setup complete\n")
cat("Seurat version:", as.character(packageVersion("Seurat")), "\n")
