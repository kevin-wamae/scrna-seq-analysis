# ****************************************************************************#
# STEP 6.3: Load clustering checkpoints
# ****************************************************************************#


# --- PREREQUISITES (scripts that must run before this one in this session) ---
# ****************************************************************************#
#   • part-3.0-install-packages.R — packages for library() calls.
#   • part-3.1-load-libraries-and-configuration.R — RUN_ID, CHECKPOINT_FORMAT,
#     DATA_CHECKPOINT_DIR, and LOG_STEP used to load the checkpoints below.
#   The dependency manager is always re-sourced here (base-R only, cheap) so
#   the latest DEP_TABLE is loaded on every run; ensure_dependencies() does
#   nothing if those objects are already in the environment, otherwise it
#   offers to source them (interactive) or stops with a clear message
#   (non-interactive/batch). When everything is already satisfied, an
#   interactive session is asked whether to re-source fresh or proceed as-is.
source("1_scripts/part-3_integration-and-clustering/part-3.0-dependencies.R")
ensure_dependencies(step = "part-6.3-load-clustering-checkpoints.R")

# NOTE: The checkpoints themselves must already exist on disk. They are
#       written by Step 6.2B (after Steps 6.1 and 6.2A have run). This script loads
#       them — it does not re-run the expensive clustering or silhouette
#       computation.


# --- WHY THIS STEP EXISTS ---
# ****************************************************************************#
#   Reaching the clustered final object is the slowest segment of Part 6: on
#   a fresh session it requires re-sourcing Steps 5.3B-6.2B, which means
#   reloading the integration checkpoints, recomputing mixing metrics, and
#   re-running both the five-resolution clustering sweep (6.1) and the
#   silhouette loop (6.2B). Step 6.2B already persisted everything downstream
#   needs to disk:
#
#     - `06_clustered_final`     : the final `integrated_final` (all
#       `clusters_res_*` columns AND `seurat_clusters` at the chosen
#       resolution).
#     - `07_clustering_metrics`  : `reduction_final`, `resolutions`,
#       `optimal_resolution`, `resolution_comparison`.
#
#   This step loads both back into memory in seconds and restores the exact
#   session state Step 6.2B left behind — so Steps 6.4, 6.5, and the final
#   save never have to recompute anything. This mirrors how Step 5.1 loads
#   the integration checkpoints so the comparison steps never re-run the
#   producers (Steps 3.3A, 4.2A-D).
#
#   CONSUMER NOTE: downstream steps (6.4, and the later final-visualization
#   and save steps) list THIS step as their prerequisite rather than 6.1/6.2B,
#   which is what keeps their in-session dependency chain short. Step 6.2A
#   (visualization) deliberately still requires 6.1 — it is upstream of the
#   checkpoints, so it has nothing to load from them.


# --- 1. Define the Checkpoints to Load ---
# ****************************************************************************#
#   The same named-vector approach as Step 5.1's `CHECKPOINT_NAMES`, and the
#   same `CHECKPOINT_FORMAT` switch (1 = qs2, 2 = rds) set once in Step 3.1,
#   so the loader honours whatever serialization the run is configured for.
CHECKPOINT_EXT <- if (CHECKPOINT_FORMAT == 1) {
    "qs2"
} else if (CHECKPOINT_FORMAT == 2) {
    "rds"
} else {
    stop("CHECKPOINT_FORMAT must be 1 (qs2) or 2 (rds)")
}

# Helper to load a single checkpoint, failing loudly with a build hint if the
# file is missing on disk (rather than letting qs2/readRDS throw a cryptic
# "cannot open file" deep inside this script).
load_checkpoint <- function(slug) {
    path <- file.path(DATA_CHECKPOINT_DIR, paste0(slug, ".", CHECKPOINT_EXT))
    if (!file.exists(path)) {
        stop(
            "Clustering checkpoint not found: ", path,
            "\nRun Steps 6.1 -> 6.2B first to build it (Step 6.2B writes this file)."
        )
    }
    if (CHECKPOINT_FORMAT == 1) qs2::qs_read(path) else readRDS(path)
}


# --- 2. Load Both Checkpoints & Restore Session State ---
# ****************************************************************************#
checked_clusters <- LOG_STEP("Loading clustering checkpoints...", {
    list(
        object  = load_checkpoint("06_clustered_final"),
        metrics = load_checkpoint("07_clustering_metrics")
    )
})

integrated_final     <- checked_clusters$object
reduction_final      <- checked_clusters$metrics$reduction_final
resolutions          <- checked_clusters$metrics$resolutions
optimal_resolution   <- checked_clusters$metrics$optimal_resolution
resolution_comparison<- checked_clusters$metrics$resolution_comparison

# `seurat_clusters` was saved at `optimal_resolution` in Step 6.2B; re-assert
# it as the active identity in case the object carried a different Idents
# state on disk.
Idents(integrated_final) <- "seurat_clusters"

cat("✓ Clustering checkpoints loaded\n")
cat(sprintf(
    "  • %d clusters at resolution %.1f (optimal), from %s UMAP\n",
    nlevels(Idents(integrated_final)), optimal_resolution, reduction_final
))


# ****************************************************************************#
# SUMMARY & PIPELINE MILESTONE TRANSITION
# ****************************************************************************#
# WHERE WE STARTED:
#   Step 6.2B had written the final clustered object (`06_clustered_final`)
#   and its supporting metrics (`07_clustering_metrics`) to disk, but a fresh
#   session had no fast, single way to bring them back — any downstream step
#   would otherwise have re-sourced the entire expensive 5.3B-6.2B chain.
#
# WHAT WE HAVE ACCOMPLISHED:
#   We loaded both checkpoints (honouring `CHECKPOINT_FORMAT`) and restored
#   `integrated_final`, `reduction_final`, `resolutions`,
#   `optimal_resolution`, and `resolution_comparison` into the session,
#   re-asserting `seurat_clusters` as the active identity. Downstream steps
#   now reach this state in seconds instead of re-running the heavy chain.
#
# WHERE WE ARE HEADING (NEXT: CLUSTER QUALITY ASSESSMENT, THEN FINAL
# VISUALIZATION):
#   Step 6.4 assesses cluster quality/stability on this restored object
#   (cluster sizes, sample/condition composition), then guide §7.8's final
#   integrated visualization and save conclude Part 3. As long as the
#   checkpoints exist, none of those steps need to recompute clustering.
# ****************************************************************************#