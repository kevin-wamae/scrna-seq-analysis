# ****************************************************************************#
# STEP 1: Load libraries and configure environment
# ****************************************************************************#


# Single-cell RNA-seq data arrives as a massive, noisy matrix of counts. Before
# jumping into biological clustering, the data must pass through a linear
# sequence of "clean-up" steps to strip away technical and physical noise.
# Each library loaded below handles highly specific stages of this purification:
#   1. Empty Droplet Filtering -> DropletUtils (Identifies soup vs cells)
#   2. Ambient Noise Scrubbing -> SoupX (Subtracts background mRNA drift)
#   3. Artificial Multi-cell Eradication -> scDblFinder (Excisise doublets)
#   4. Downstream Downstream Integration -> Seurat (Core structural toolkit)


# --- REPRODUCIBILITY MECHANICS ---
# ****************************************************************************#
# `set.seed(100)` pins the pseudorandom number generator. Complex algorithms
# (like UMAP layout dimensions or scDblFinder in-silico synthetic blending)
# utilize stochastic processes. Setting a static seed guarantees that any
# analysis performed anywhere will generate identical data clusters.
set.seed(100)


# Terminal output colouring (load first so errors/warnings are coloured)
library(colorout)


# --- Core Single-Cell Infrastructure ---
# ****************************************************************************#
# - `Seurat`: The industry-standard framework for managing, normalising, and
#   clustering single-cell datasets.
# - `SeuratObject`: Houses the underlying object architecture data models.
library(Seurat)
library(SeuratObject)


# --- QC-Specific Packages ---
# ****************************************************************************#
# - `DropletUtils`: Provides the statistical framework (`EmptyDrops`) needed to
#   salvage small, genuine cell types from cell-free droplets.
# - `SoupX`: Profiles background mRNA profiles and applies mathematical
#   subtraction to cleanse ambient contamination from true cell data layers.
# - `scDblFinder`: Generates in-silico hybrid cell profiles to train a model
#   capable of finding and removing hidden cellular multiplexes (doublets).
# - `SingleCellExperiment`: A standard Bioconductor object model required by
#   bioconductor packages (SoupX, scDblFinder) to read and compute fast matrices.
library(DropletUtils)
library(scDblFinder)
library(SoupX)
library(SingleCellExperiment)


# --- Data Visualization & Manipulation ---
# ****************************************************************************#
# - `ggplot2`: The foundational graphics syntax library engine used to build
#   and customize publication-ready quality distribution plots.
# - `patchwork`: Combines separate plot structures together into single panels.
# - `dplyr`: Provides optimized C++ syntax for lightning-fast metadata column
#   manipulation, cell annotations, and metadata row mutation.
# - `scales`: Provides x/y axis scaling/formatting for ggplot2.
library(ggplot2)
library(patchwork)
library(dplyr)
library(scales)


# --- REPRODUCIBILITY MECHANICS ---
# ****************************************************************************#
# `set.seed(100)` pins the pseudorandom number generator. Complex algorithms
# (like UMAP layout dimensions or scDblFinder in-silico synthetic blending)
# utilize stochastic processes. Setting a static seed guarantees that any
# analysis performed anywhere will generate identical data clusters.
set.seed(100)


# --- SEURAT 5 COMPUTE FRAMEWORKS ---
# ****************************************************************************#
# SEURAT 5 VS OLD ARCHITECTURES:
# Seurat 5 introduced a complete redesign of data slot storage known as the
# "Layer Architecture".
#   - Legacy objects store data in static, monolithic matrices (`counts`,
#     `data`, `scale.data`) which easily trigger out-of-memory RAM crashes.
#   - Seurat 5 objects split matrices dynamically into granular multi-layered
#     slots. This allows your script to retain uncorrected raw counts in one
#     independent layer while applying background mathematical corrections in
#     another parallel layer within the exact same object footprint.

if (packageVersion("Seurat") >= "5.0.0") {
    cat("✓ Seurat 5 detected - layer-based architecture available\n")
} else {
    warning("Please upgrade to Seurat 5 for this tutorial")
}
