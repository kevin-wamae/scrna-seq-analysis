# ****************************************************************************#
# STEP 6.2B: Evaluate optimal clustering resolution
# ****************************************************************************#


# --- PREREQUISITES (scripts that must run before this one in this session) ---
# ****************************************************************************#
#   • part-3.0-install-packages.R — packages for library() calls (`cluster`,
#     loaded in Step 3.1, supplies `silhouette()`).
#   • part-3.1-load-libraries-and-configuration.R — RUN_ID, LOG_STEP, and
#     output-directory variables.
#   • part-5.3B-select-best-integration-method.R — supplies `reduction_final`,
#     used below to look up the matching integrated (pre-UMAP) embedding.
#   • part-6.1-clustering-multiple-resolutions.R — supplies `integrated_final`
#     with its `clusters_res_<res>` columns, and `resolutions`, the vector of
#     resolutions scored below.
#   The dependency manager is always re-sourced here (base-R only, cheap) so
#   the latest DEP_TABLE is loaded on every run; ensure_dependencies() does
#   nothing if those objects are already in the environment, otherwise it
#   offers to source them (interactive) or stops with a clear message
#   (non-interactive/batch). When everything is already satisfied, an
#   interactive session is asked whether to re-source fresh or proceed as-is.
source("1_scripts/part-3_integration-and-clustering/part-3.0-dependencies.R")
ensure_dependencies(step = "part-6.2B-evaluate-optimal-clustering-resolution.R")

# NOTE: Requires Steps 5.3B and 6.1 to have already run in this session.


# --- THE SILHOUETTE SCORE METRIC (FOR BEGINNERS) ---
# ****************************************************************************#
#   Step 6.2A gave us a visual first impression of the five resolutions —
#   useful, but "this one looks about right" isn't something you can report
#   or compare objectively. The silhouette score is the quantitative
#   counterpart, used here to assess cluster quality at each resolution.
#
#   WHAT IT MEASURES: for each cell, we calculate its average distance to
#   other cells in its OWN cluster versus its average distance to cells in
#   the nearest DIFFERENT cluster, then average that comparison across all
#   cells to get one score per resolution. Scores range from -1 to +1.
#
#   INTERPRETATION:
#     - Strong separation (>0.5)     : clusters are well-separated, cells
#       clearly belong to their clusters — good ✓
#     - Moderate separation (0.3-0.5): acceptable clusters, some overlap
#       between groups
#     - Weak separation (<0.3)       : clusters poorly defined, may
#       indicate over-clustering
#
#   THE TRADE-OFF:
#     - Too few clusters (low resolution)  → biology oversimplified,
#       different cell types lumped together
#     - Too many clusters (high resolution) → over-splitting of natural
#       cell populations
#
#   SELECTION STRATEGY: choose the resolution with the highest silhouette
#   score WHILE maintaining a biologically interpretable cluster number —
#   the two criteria together, not the score alone. As a rough benchmark
#   for PBMC data like ours, 8-15 clusters is typical once major cell types
#   and their subtypes are both accounted for; a "winning" resolution that
#   falls far outside that range is worth treating as a flag to double-check
#   against Step 6.2A's UMAP grid, not a reason to override the score
#   automatically.
#
#   OUR STRATEGY, STEP BY STEP:
#     1. Test multiple resolutions — the sweep already done in Step 6.1.
#     2. Calculate a silhouette score for each, using a 5,000-cell
#        subsample (below).
#     3. Look for the resolution with a high silhouette score AND an
#        interpretable cluster number — not just whichever score happens
#        to be numerically highest.
#     4. Visually validate that the chosen resolution makes biological
#        sense — this is what Step 6.2A's comparison grid was for, and is
#        worth revisiting once a candidate resolution is picked here.
#
#   NUMERIC RECOMMENDATION, NOT BIOLOGICAL VERDICT:
#   The silhouette score supplies a reproducible comparison of separation in
#   the integrated space. It cannot determine whether a cluster represents a
#   real cell type, a cell state, or an over-split continuum. The
#   broad PBMC ranges (for example, roughly 8-15 clusters for major types plus
#   subtypes) are orientation only; they should not be used as a universal
#   pass/fail filter or forced onto a different cohort.
#
#   Therefore the resolution printed below is best read as the top numeric
#   candidate. A defensible final choice also revisits the Step 6.2A UMAP grid,
#   cluster sizes and sample composition from Step 6.4, and—once annotation
#   begins—marker-gene coherence. A silhouette value below the rough
#   0.3 reference is a prompt to investigate overlap or over-splitting, not an
#   automatic rejection; values depend on the embedding, distance scale, and
#   biological complexity of the dataset.
#
#   WHY THE INTEGRATED EMBEDDING, NOT THE UMAP: silhouette scores are
#   calculated on `integrated_reduction` (e.g. "harmony", "integrated.cca")
#   rather than the 2D UMAP used for plotting. UMAP is a visualisation
#   layout optimised to look good in 2 dimensions, which can distort true
#   distances; the 30+-dimension integrated embedding is the actual space
#   clustering was performed in, so it is what the score should measure.


# --- 1. Map the UMAP Reduction Back to Its Integrated Embedding ---
# ****************************************************************************#
#   `reduction_final` (set in Step 5.3B) names the UMAP layout used for
#   plotting; the silhouette score instead needs the higher-dimensional
#   embedding that UMAP itself was built from. This lookup mirrors the
#   `reduction.name` given to each method's `RunUMAP()` call back in
#   Section 4 (e.g. "umap.cca" was built from "integrated.cca").
integrated_reduction <- switch(reduction_final,
    "umap.cca"     = "integrated.cca",
    "umap.rpca"    = "integrated.rpca",
    "umap.harmony" = "harmony",
    "umap.mnn"     = "integrated.mnn",
    "umap"         = "pca"   # Fallback for the naive merge
)

cat("\nEvaluating clustering resolutions using", integrated_reduction, "space\n")


# --- 2. Subsample for Computational Feasibility ---
# ****************************************************************************#
#   Silhouette scoring needs a full pairwise distance matrix between cells,
#   which scales quadratically with cell count — at our ~72,000 cells that
#   would require over 25 GB of memory just to hold that matrix, making it
#   computationally prohibitive. Subsampling to 5,000 cells (our strategy's
#   step 2, above) keeps the calculation tractable — this is standard
#   practice, not a shortcut that compromises the result: silhouette scores
#   from a well-chosen subsample reliably reflect the separation quality of
#   the full clustering.
#
#   NOTE ON THE LOCAL SEED: `set.seed(100)` was already called once, in Step
#   3.1, to seed the pipeline as a whole. The `set.seed(42)` here is a
#   deliberate, separate seed scoped only to this subsample draw, so the
#   exact 5,000 cells chosen are reproducible on their own regardless of
#   how much other random state the rest of the session has consumed by
#   this point.
n_cells <- ncol(integrated_final)
max_cells_for_silhouette <- 5000

if (n_cells > max_cells_for_silhouette) {
    cat("Dataset has", n_cells, "cells - subsampling",
        max_cells_for_silhouette, "cells for silhouette calculation\n")
    set.seed(42)
    subsample_idx <- sample(1:n_cells, max_cells_for_silhouette)
} else {
    subsample_idx <- 1:n_cells
}


# --- 3. Calculate Silhouette Score for Each Resolution ---
# ****************************************************************************#
#   One score per resolution, computed on the SAME subsampled cells and the
#   SAME integrated embedding each time, so the five scores are directly
#   comparable to one another. Capped at the first 30 dimensions of the
#   embedding to match the dimensionality used during clustering itself
#   (`dims = 1:30` throughout Section 4/5).
silhouette_scores <- LOG_STEP("Calculating silhouette scores across resolutions...", {
    sapply(resolutions, function(res) {
        cluster_col <- paste0("clusters_res_", res)

        # Cluster assignments for the subsampled cells
        clusters <- as.numeric(integrated_final@meta.data[[cluster_col]][subsample_idx])

        # Integrated coordinates for the same cells (first 30 dims)
        coords <- Embeddings(integrated_final, reduction = integrated_reduction)[subsample_idx, ]
        coords <- coords[, 1:min(30, ncol(coords))]

        # Silhouette is undefined with a single cluster — guard against that
        if (length(unique(clusters)) > 1) {
            dist_matrix <- dist(coords)
            sil <- cluster::silhouette(clusters, dist_matrix)
            mean(sil[, 3])
        } else {
            NA   # Only one cluster at this resolution
        }
    })
})


# --- 4. Build & Save the Resolution Comparison Table ---
# ****************************************************************************#
#   One row per resolution: its cluster count (from Step 6.1) alongside its
#   silhouette score (calculated above), so both numbers driving the
#   eventual choice sit side by side rather than in separate objects.
resolution_comparison <- data.frame(
    resolution = resolutions,
    n_clusters = sapply(resolutions, function(res) {
        cluster_col <- paste0("clusters_res_", res)
        length(unique(integrated_final@meta.data[[cluster_col]]))
    }),
    silhouette_score = silhouette_scores
)

cat("\nClustering resolution comparison:\n")
print(resolution_comparison)

write.csv(resolution_comparison,
    file.path(METADATA_OUT_DIR, "resolution_comparison_metrics.csv"),
    row.names = FALSE
)
cat("→ Saved:", file.path(METADATA_OUT_DIR, "resolution_comparison_metrics.csv"), "\n")


# --- 5. Recommend the Optimal Resolution ---
# ****************************************************************************#
#   Strategy step 3 above calls for the resolution with a high silhouette
#   score AND an interpretable cluster number — not simply whichever score
#   is numerically highest. `which.max()` is used as a starting point
#   because it is deterministic and reproducible, but the printed table
#   above is there so this pick can be sanity-checked by eye: for PBMC data
#   like ours, 8-15 clusters is the typical interpretable range (major cell
#   types plus subtypes), so a top-scoring resolution landing well outside
#   that range is worth a second look before trusting it outright.
#
#   Strategy step 4 — visually validating that the chosen resolution makes
#   biological sense — is not re-done here; it is exactly what Step 6.2A's
#   comparison grid already gives you. Revisit that figure for whichever
#   resolution is recommended below before treating it as final.
#
#   `Idents()` and `seurat_clusters` are set to this resolution's column so
#   that any code downstream that expects Seurat's usual "current active
#   clustering" (rather than an explicit `clusters_res_<res>` column) picks
#   up the recommended resolution by default.
optimal_idx <- which.max(resolution_comparison$silhouette_score)
optimal_resolution <- resolution_comparison$resolution[optimal_idx]

cat("\nRecommended resolution:", optimal_resolution,
    "(", resolution_comparison$n_clusters[optimal_idx], "clusters )\n")
cat("→ Sanity-check this against 06_multi_resolution_clustering.png",
    "(Step 6.2A) before treating it as final.\n\n")

integrated_final$seurat_clusters <- integrated_final@meta.data[[paste0("clusters_res_", optimal_resolution)]]
Idents(integrated_final) <- "seurat_clusters"


# --- 6. Checkpoint the Final Clustered Object & Selection ---
# ****************************************************************************#
#   Reaching this point is the most expensive segment of the pipeline after
#   integration: Steps 5.3B-6.1 require reloading the integration
#   checkpoints, recomputing mixing metrics, and re-clustering at five
#   resolutions, and the silhouette loop above is the slowest single
#   computation in Part 6. Everything downstream needs — the clustered
#   object, the chosen resolution, its supporting objects — is now in this
#   session, so we persist two checkpoints honoring `CHECKPOINT_FORMAT`
#   (set once in Step 3.1: 1 = qs2, 2 = rds):
#
#     1. `06_clustered_final` — `integrated_final` exactly as it stands
#        right now: every `clusters_res_*` column from Step 6.1 AND
#        `seurat_clusters` set to the chosen `optimal_resolution`.
#     2. `07_clustering_metrics` — the small supporting objects needed to
#        reconstruct the session: `reduction_final`, `resolutions`,
#        `optimal_resolution`, and `resolution_comparison`. Kept as a
#        separate, light file so loading back from disk never requires
#        re-reading a multi-GB object just to recover a few scalars.
#
#   Step 6.3 exists to load both back; with them on disk this script only
#   ever needs to run once.
clustering_metrics <- list(
    reduction_final       = reduction_final,
    resolutions           = resolutions,
    optimal_resolution    = optimal_resolution,
    resolution_comparison = resolution_comparison
)

if (CHECKPOINT_FORMAT == 1) {
    CLUSTERED_CP  <- file.path(DATA_CHECKPOINT_DIR, "06_clustered_final.qs2")
    METRICS_CP    <- file.path(DATA_CHECKPOINT_DIR, "07_clustering_metrics.qs2")
    LOG_STEP("Saving clustered object + metrics checkpoints (qs2)...", {
        qs2::qs_save(integrated_final, CLUSTERED_CP, nthreads = N_WORKERS)
        qs2::qs_save(clustering_metrics, METRICS_CP, nthreads = N_WORKERS)
    })
} else if (CHECKPOINT_FORMAT == 2) {
    CLUSTERED_CP  <- file.path(DATA_CHECKPOINT_DIR, "06_clustered_final.rds")
    METRICS_CP    <- file.path(DATA_CHECKPOINT_DIR, "07_clustering_metrics.rds")
    LOG_STEP("Saving clustered object + metrics checkpoints (rds, single-threaded)...", {
        saveRDS(integrated_final, CLUSTERED_CP)
        saveRDS(clustering_metrics, METRICS_CP)
    })
} else {
    stop("CHECKPOINT_FORMAT must be 1 (qs2) or 2 (rds)")
}

cat("✓ Checkpoint written:", CLUSTERED_CP, "\n")
cat("✓ Checkpoint written:", METRICS_CP, "\n")


# ****************************************************************************#
# SUMMARY & PIPELINE MILESTONE TRANSITION
# ****************************************************************************#
# WHERE WE STARTED:
#   Step 6.2A gave us a visual, side-by-side comparison of five resolutions —
#   persuasive as a first look, but not a number that could be reported or
#   used to make a decision on its own.
#
# WHAT WE HAVE ACCOMPLISHED:
#   We calculated a silhouette score per resolution on a reproducible
#   5,000-cell subsample in the integrated embedding, saved the full
#   resolution/cluster-count/silhouette-score table to
#   resolution_comparison_metrics.csv, and set `integrated_final`'s active
#   clustering (`seurat_clusters` / `Idents()`) to `optimal_resolution` —
#   the top-scoring candidate, flagged here for a visual sanity check
#   against Step 6.2A's comparison grid rather than accepted blind. We then
#   checkpointed the final clustered object (`06_clustered_final`) together
#   with its supporting metrics (`07_clustering_metrics`) to
#   `DATA_CHECKPOINT_DIR`, so this expensive step only ever runs once.
#
# WHERE WE ARE HEADING (NEXT: STEP 6.3 LOADER, THEN CLUSTER QUALITY):
#   Step 6.3 loads the two checkpoints back into memory in seconds,
#   letting every downstream step (cluster quality assessment, the final
#   integrated visualization, and the eventual save) pick
#   up exactly this state without re-running Steps 5.3B-6.2B. From there
#   Step 6.4 checks cluster sizes and each cluster's sample composition at
#   `optimal_resolution`, flagging anything that still looks sample-
#   dominated before this clustering is treated as final.
# ****************************************************************************#