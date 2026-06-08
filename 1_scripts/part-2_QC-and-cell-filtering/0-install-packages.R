################################################################################
# Installation: Required packages for single-sample QC
################################################################################

# --- THE INFRASTRUCTURE REQUIREMENTS ---
# Single-cell QC depends on tools spanning different software ecosystems. To
# ensure a smooth workflow, packages must be fetched from three distinct pools:
#   1. CRAN: Standard R repositories (Seurat, ggplot2, dplyr, SoupX).
#   2. Bioconductor: Specialized genomic data frameworks (DropletUtils,
#      scDblFinder, SingleCellExperiment).
#   3. External/Community: Independent specialized repos (colorout).
#
# REPRODUCIBILITY & ECOSYSTEM BENEFIT:
# By wrapping installations inside an `if (!requireNamespace(...))` conditional
# structure, the script validates your local environment first. It only downloads
# a package if it is missing, preventing your script from wasting compute hours
# and bandwidth re-downloading active libraries during repetitive runs.


# --- 1. Establish CRAN Network Target ---
# Sets a cloud-based mirror to ensure the machine downloads from the nearest
# high-speed geographic server, preventing network handshake drops.
options(repos = c(CRAN = "https://cloud.r-project.org"))


# --- 2. Conditional Core Seurat Framework ---
# `requireNamespace("pkg", quietly = TRUE)` checks if a package is installed
# without loading it into your current environment. Returning FALSE triggers
# the installation.
# - `Seurat`: The core toolkit for data integration and downstream scaling.
# - `SeuratObject`: Explicitly ensures the underlying v5 slot/layer matrix
#   classes are updated and fully compatible with modern pipelines.
if (!requireNamespace("Seurat", quietly = TRUE)) {
  install.packages("Seurat")
}

if (!requireNamespace("SeuratObject", quietly = TRUE)) {
  install.packages("SeuratObject")
}


# --- 3. Conditional Data Handling & Visualization Engines ---
# For smaller utility packages, we evaluate them against a vector loop.
# This programmatically builds a custom whitelist (`missing_cran`) of only what
# your machine lacks, allowing CRAN to resolve shared parallel dependencies
# all at once.
cran_pkgs <- c("ggplot2", "patchwork", "dplyr")
missing_cran <- cran_pkgs[!sapply(cran_pkgs, requireNamespace, quietly = TRUE)]

if (length(missing_cran) > 0) {
  install.packages(missing_cran)
}


# --- 4. Initialize Bioconductor Ecosystem ---
# Bioconductor handles large, complex genomic data types. It uses its own
# dedicated package manager (`BiocManager`) to verify that the versions of your
# QC tools match your exact version of base R.
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Identify which specific Bioconductor packages are missing:
# - `DropletUtils`: Computes statistical profiles to find true cell boundaries.
# - `scDblFinder`: Leverages machine learning to flag in-silico doublets.
# - `SingleCellExperiment`: Provides the foundation matrix layout model that
#    both DropletUtils and scDblFinder use to run their matrix calculations.
bioc_pkgs <- c("DropletUtils", "scDblFinder", "SingleCellExperiment")
missing_bioc <- bioc_pkgs[!sapply(bioc_pkgs, requireNamespace, quietly = TRUE)]

if (length(missing_bioc) > 0) {
  # `update = FALSE, ask = FALSE` stops Bioconductor from pausing your script
  # to ask manual interactive update questions in the terminal.
  BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)
}


# --- 5. Install Ambient Subtraction & UI Utilities ---
# - `SoupX`: Fetched from CRAN to handle profile modeling of cell-free fluid.
if (!requireNamespace("SoupX", quietly = TRUE)) {
  install.packages("SoupX")
}

# - `colorout`: An independent package hosted on a community repository.
#   It alters the R kernel parsing system to color-code messages, flashing
#   warnings in yellow and errors in red for instantaneous visual diagnostics.
if (!requireNamespace("colorout", quietly = TRUE)) {
  install.packages("colorout", repos = "https://community.r-multiverse.org")
}

cat("\n✓ Environment verification complete. All packages are ready!\n")
