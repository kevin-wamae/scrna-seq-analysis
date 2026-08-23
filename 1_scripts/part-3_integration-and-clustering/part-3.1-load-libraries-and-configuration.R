# ****************************************************************************#
# STEP 3.1: Load required libraries and configure environment
# ****************************************************************************#


# --- DECLARE PIPELINE VARIABLES ---
# ****************************************************************************#
# The date-stamped job run ID (from HPC scheduler) for tracking. Every output
# path constructed later in this script nests under this ID, so all files
# from a given cluster job land in one traceable folder tree.
RUN_ID <- "2026_06_09_brown_job_3058993"


# --- REPRODUCIBILITY MECHANICS ---
# ****************************************************************************#
# `set.seed(100)` pins the pseudorandom number generator before anything in
# this pipeline can draw on it. Harmony's stochastic optimization, FastMNN's
# nearest-neighbour searches, and UMAP layout all rely on random draws — set
# this early, once, so every downstream step (including parallel workers
# spawned later via `future`) inherits a reproducible baseline.
set.seed(100)


# --- TERMINAL OUTPUT COLOURING ---
# ****************************************************************************#
# Loaded first so that any errors/warnings thrown by the libraries below are
# already colour-coded when they print.
library(colorout)


# --- CORE SINGLE-CELL INFRASTRUCTURE ---
# ****************************************************************************#
# - `Seurat`: The core toolkit for managing, normalising, and clustering
#   single-cell datasets.
# - `SeuratObject`: Houses the underlying object architecture data models.
library(Seurat)
library(SeuratObject)


# --- INTEGRATION METHODS ---
# ****************************************************************************#
# - `harmony`: Fast, probabilistic batch correction via iterative clustering.
# - `batchelor`: Provides FastMNN, used internally by Seurat's
#   `IntegrateLayers()` as an alternative integration method.
# - `SeuratWrappers`: Community-maintained wrapper functions that let Seurat
#   call third-party integration methods (Harmony, FastMNN, etc.) through a
#   consistent `IntegrateLayers()` interface.
library(harmony)
library(batchelor)
library(SeuratWrappers)


# --- DATA VISUALIZATION & MANIPULATION ---
# ****************************************************************************#
# - `ggplot2` / `patchwork`: Core plotting engine and multi-panel assembly.
# - `ggrepel`: Non-overlapping text labels on crowded plots (e.g. marker
#   gene labels on variable feature plots).
# - `RColorBrewer` / `viridis`: Colour palettes for categorical and
#   perceptually-uniform continuous scales, respectively.
# - `dplyr`: Fast, readable metadata manipulation and filtering.
# - `reshape2`: Reshaping data (wide/long) for custom cluster quality plots.
library(ggplot2)
library(ggrepel)
library(dplyr)
library(patchwork)
library(RColorBrewer)
library(viridis)
library(reshape2)


# --- QUALITY METRICS & UTILITIES ---
# ****************************************************************************#
# - `FNN`: K-nearest neighbour calculations, used for batch-mixing metrics
#   that assess how well integration removed technical separation.
# - `cluster`: Silhouette scores, used to assess clustering quality/cohesion.
library(FNN)
library(cluster)


# --- PARALLEL PROCESSING ---
# ****************************************************************************#
# - `future`: Parallel/asynchronous processing backend.
# - `furrr` / `purrr`: purrr-flavored parallel map() functions built on
#   `future`, used for parallelizing per-sample operations (e.g. loading
#   QC-filtered `.rds` files in Step 2).
library(future)
library(furrr)
library(purrr)

# `future.globals.maxSize` caps how much data `future` can serialize and
# export to each parallel worker in one go. 20GB is a placeholder — right-size
# this to your actual merged object size once samples are loaded (Step 2):
#   OBJ_SIZE_GB <- as.numeric(object.size(SEURAT_OBJ)) / 1024^3
#   options(future.globals.maxSize = ceiling(OBJ_SIZE_GB * 3) * 1024^3)
# and keep it comfortably under whatever --mem your srun/sbatch job requests.
options(future.globals.maxSize = 20 * 1024^3)


# --- PIPELINE LOGGING UTILITY ---
# ****************************************************************************#
# Several steps in this pipeline (integration, clustering, UMAP on large
# matrices) can run for several minutes with no console output. Without
# feedback, this looks indistinguishable from a hang. LOG_STEP() prints a
# start message, runs the supplied code, times it, and prints a "Done in Xs"
# message when it finishes — flushing stdout immediately so messages appear
# on screen right away, even when running non-interactively (Rscript, slurm
# logs), where R would otherwise buffer output until the call finishes.
#
# USAGE:
#   SEURAT_INTEGRATED <- LOG_STEP("Running Harmony integration...", {
#     IntegrateLayers(SEURAT_OBJ, method = HarmonyIntegration)
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

# Checkpoint archive: intermediate objects (merged_naive, merged_seurat,
# and the four per-method integrated objects) saved so a crashed/killed
# HPC job can resume mid-pipeline without re-running earlier integrations.
# Kept separate from DATA_OUT_DIR so `integrated_data/` holds only the
# final deliverable, and this folder can be deleted once a run is validated.
DATA_CHECKPOINT_DIR <- file.path(
  "3_output", RUN_ID, "integration_and_clustering", "checkpoints"
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
dir.create(DATA_CHECKPOINT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("✓ Output directories securely synchronized to disk:\n")
cat("  • Visual Diagnostic Canvas:", PLOTS_OUT_DIR, "\n")
cat("  •   ├─ Integration Comparison:", PLOTS_COMPARISON_DIR, "\n")
cat("  •   └─ Clustering:", PLOTS_CLUSTERING_DIR, "\n")
cat("  • Integrated Data Archive:", DATA_OUT_DIR, "\n")
cat("  • Metadata Ledger:", METADATA_OUT_DIR, "\n")
cat("  • Checkpoint Vault:", DATA_CHECKPOINT_DIR, "\n\n")


# --- PLOTTING DEFAULTS ---
# ****************************************************************************#
# Applies a consistent, publication-friendly theme across every plot
# generated downstream, so individual plotting blocks don't need to repeat it.
theme_set(theme_classic(base_size = 12))


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


# --- ENVIRONMENT VERIFICATION ---
# ****************************************************************************#
cat("Environment setup complete\n")
cat("Seurat version:", as.character(packageVersion("Seurat")), "\n")