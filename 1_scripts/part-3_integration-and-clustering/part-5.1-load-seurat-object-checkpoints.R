# ****************************************************************************#
# STEP 5.1: Load Seurat object checkpoints
# ****************************************************************************#

# --- WHY THIS STEP EXISTS ---
# ****************************************************************************#
#   Steps 3.3A, 4.2A-D produced five checkpointed objects on disk: the naive,
#   uncorrected baseline plus four independently-run integration methods
#   (CCA, RPCA, Harmony, FastMNN). Both the visual comparison (Step 5.2A)
#   and the quantitative mixing metrics (Step 5.2B) need all five of these
#   objects loaded into memory, so the loading logic lives here once instead
#   of being duplicated — and potentially drifting out of sync — across both
#   downstream scripts.
#
#   Pulling this into its own script also makes failures easier to debug: if
#   a checkpoint is missing or corrupted, that surfaces immediately in a
#   focused, single-purpose step rather than partway through a 200+ line
#   plotting or metrics script where the cause is less obvious.


# --- 1. Load All Five Checkpointed Objects ---
# ****************************************************************************#
#   Loads the naive baseline (Step 3.3A) plus all four integration methods
#   (Steps 4.2A-D) back from the checkpoint vault built in Step 3.1
#   (`DATA_CHECKPOINT_DIR`).
#
#   FORMAT HANDLING: `CHECKPOINT_FORMAT` (set once in Step 3.1: 1 = qs2,
#   2 = rds) governs every checkpoint read/write in this pipeline. This
#   loader respects that same switch rather than hardcoding `qs2::qs_read()`,
#   so a run configured for `rds` doesn't silently fail here.

# Create a named vector of checkpoint names to load
CHECKPOINT_NAMES <- c(
    naive   = "01_merged_naive",
    cca     = "02_integrated_cca",
    rpca    = "03_integrated_rpca",
    harmony = "04_integrated_harmony",
    fastmnn = "05_integrated_fastmnn"
)

# Get the checkpoint extension based on the checkpoint format
CHECKPOINT_EXT <- if (CHECKPOINT_FORMAT == 1) {
    "qs2"
} else if (CHECKPOINT_FORMAT == 2) {
    "rds"
} else {
    stop("CHECKPOINT_FORMAT must be 1 (qs2) or 2 (rds)")
}

# Helper function to load a checkpoint
load_checkpoint <- function(name) {
    path <- file.path(DATA_CHECKPOINT_DIR, paste0(name, ".", CHECKPOINT_EXT))
    if (CHECKPOINT_FORMAT == 1) {
        qs2::qs_read(path)
    } else {
        readRDS(path)
    }
}

# Load all checkpoints
CHECKPOINTS <- LOG_STEP("Loading all 5 checkpointed objects for comparison...", {
    setNames(
        lapply(CHECKPOINT_NAMES, load_checkpoint),
        names(CHECKPOINT_NAMES)
    )
})


# --- 2. Extract Objects into Named Variables ---
# ****************************************************************************#
#   Unpacks the `CHECKPOINTS` list into the individual object names
#   (`merged_naive`, `integrated_cca`, etc.) that Steps 5.2A and 5.2B
#   reference directly, so downstream scripts don't need to know or care
#   that the objects were loaded via a named list in the first place.
merged_naive        <- CHECKPOINTS$naive
integrated_cca      <- CHECKPOINTS$cca
integrated_rpca     <- CHECKPOINTS$rpca
integrated_harmony  <- CHECKPOINTS$harmony
integrated_fastmnn  <- CHECKPOINTS$fastmnn

# Print summary of loaded checkpoints
cat("✓ All 5 checkpoints loaded\n\n")


# ****************************************************************************#
# SUMMARY & PIPELINE MILESTONE TRANSITION
# ****************************************************************************#
# WHERE WE STARTED:
#   Five checkpointed objects sitting on disk from Steps 3.3A and 4.2A-D,
#   with no single, reusable way to bring them all back into memory for
#   comparison.
#
# WHAT WE HAVE ACCOMPLISHED:
#   We loaded the naive baseline and all four integration methods from
#   their checkpoints (respecting whichever `CHECKPOINT_FORMAT` this run is
#   configured for) and unpacked them into the five named objects
#   (`merged_naive`, `integrated_cca`, `integrated_rpca`,
#   `integrated_harmony`, `integrated_fastmnn`) that the comparison steps
#   expect to find in memory.
#
# WHERE WE ARE HEADING (STEP 5.2A — VISUAL INTEGRATION COMPARISON, THEN
# STEP 5.2B — QUANTITATIVE INTEGRATION METRICS):
#   With all five objects loaded once here, Step 5.2A builds side-by-side
#   UMAP grids (colored by sample, then by condition) to visually judge
#   mixing and biological preservation, and Step 5.2B follows with the
#   numeric mixing-score companion to that comparison. Both steps assume
#   this script has already been run in the same session.
# ****************************************************************************#