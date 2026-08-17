# ****************************************************************************#
# STEP 1: Load libraries and configure environment
# ****************************************************************************#

# TODO: Speed up processes
# TODO: Identify processes that can be improved with bpcells package

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


# --- SET WORKER/CPU COUNT ---
# ****************************************************************************#
# Determine worker/CPU count for downstream parallel processing:
#   - If running under a job scheduler (SLURM), respect the allocation and
#     leave 1 core free for OS/monitoring overhead.
#   - Otherwise assume this is a shared interactive machine (e.g. a lab
#     workstation) and cap usage at 4 cores, regardless of how many the
#     machine physically has, to avoid starving other users.

# Count number of CPUs available in SLURM scheduler
SLURM_ALLOC <- Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA)


# Set number of workers/CPUs depending on environment
if (!is.na(SLURM_ALLOC)) {
    N_WORKERS <- max(1, as.integer(SLURM_ALLOC) - 1)
    cat("Scheduler detected (SLURM) - using", N_WORKERS, "worker(s)\n")
} else {
    N_WORKERS <- min(4, max(1, parallel::detectCores() - 1))
    cat("No scheduler detected - shared environment cap:", N_WORKERS, "worker(s)\n")
}


# --- PIPELINE LOGGING UTILITY ---
# ****************************************************************************#
# Several steps in this pipeline (EmptyDrops, SoupX, scDblFinder, Read10X on
# large matrices) can run for several minutes with no console output. Without
# feedback, this looks indistinguishable from a hang. LOG_STEP() prints a
# start message, runs the supplied code, times it, and prints a "Done in Xs"
# message when it finishes — flushing stdout immediately so messages appear
# on screen right away, even when running non-interactively (Rscript, slurm
# logs), where R would otherwise buffer output until the call finishes.
#
# USAGE:
#   COUNTS <- LOG_STEP("Loading 10x matrix...", {
#     Read10X(data.dir = file.path(CELLRANGER_INPUT, "raw_feature_bc_matrix"))
#   })
LOG_STEP <- function(msg, expr) {
    message(msg)
    flush(stdout())

    ELAPSED <- system.time({
        RESULT <- eval.parent(substitute(expr))
    })

    message(sprintf("✓ Done in %.1f seconds", ELAPSED["elapsed"]))
    flush(stdout())

    RESULT
}
