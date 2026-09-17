# ****************************************************************************#
# STEP 6.1: Perform clustering at multiple resolutions
# ****************************************************************************#


# --- PREREQUISITES (scripts that must run before this one in this session) ---
# ****************************************************************************#
#   • part-3.0-install-packages.R — packages for library() calls.
#   • part-3.1-load-libraries-and-configuration.R — RUN_ID, LOG_STEP, and the
#     output-directory variables.
#   • part-5.3B-select-best-integration-method.R — supplies `integrated_final`
#     (the winning method's Seurat object) and `reduction_final` (its UMAP
#     reduction name), committed to as the single integration result every
#     downstream step builds on.
#   If `integrated_final` and `reduction_final` are already in the
#   environment, this block does nothing; otherwise it offers to source the
#   prerequisites (interactive) or stops with a clear message
#   (non-interactive/batch).
if (!exists("ensure_dependencies", inherits = TRUE)) {
    source("1_scripts/part-3_integration-and-clustering/part-3.0-dependencies.R")
}
ensure_dependencies(step = "part-6.1-clustering-multiple-resolutions.R")

# NOTE: Requires Step 5.3B to have already run in this session — that step is
#       what selects `integrated_final` / `reduction_final` out of
#       `methods_list`.

# ****************************************************************************#
# --- WHY THIS STEP EXISTS ---
# ****************************************************************************#
#   Every integration method was clustered at a single, arbitrary resolution
#   (0.6 — Seurat's standard default) purely so each method could carry a
#   comparable `seurat_clusters` column for the mixing and preservation
#   metrics. But that resolution is a decision, not a fact: Louvain modularity
#   clustering's `resolution` parameter directly controls how coarsely or
#   finely the shared-nearest-neighbour (SNN) graph is partitioned
#   (higher → more, smaller communities; lower → fewer, coarser ones). No
#   single value is "correct" for every dataset, so before committing to one
#   partitioning we sweep a range and record what the data actually looks like
#   at each granularity.
#
#   THE RESOLUTION / GRANULARITY TRADE-OFF:
#     - Low resolutions (~0.4): broad cell-type families (T cells, B cells,
#       monocytes) — a useful overview, but genuinely distinct types can be
#       lumped together.
#     - Mid resolutions (~0.6–0.8): the typical working range — separated
#       cell types plus coarse subtypes, stable run-to-run.
#     - High resolutions (~1.0–1.2): fine-grained cell states and rare
#       populations — useful, but be alert to over-splitting a single true
#       cell type into artificial fragments.
#   The point of the sweep is not to pick "the" right resolution here — it is
#   to generate, in one pass, the cluster assignments future steps need (UMAP
#   grids, silhouette-based optimal-resolution selection, quality checks)
#   without re-running clustering repeatedly.
#
#   NOTE ON `seurat_clusters`: `FindClusters()` overwrites the object's active
#   `seurat_clusters` column on every run, so after this step it reflects only
#   the LAST resolution (1.2). Downstream steps must therefore read the
#   explicit `clusters_res_*` columns written below, never `seurat_clusters`.


# --- 1. Define the Resolution Sweep ---
# ****************************************************************************#
#   Five evenly-spaced values spanning the coarse-to-fine working range
#   (guide §7.4). Resolution sits in log-like space, so 0.4 → 1.2 crosses the
#   typical low/mid/high boundaries of the granularity trade-off above.
resolutions <- c(0.4, 0.6, 0.8, 1.0, 1.2)


# --- 2. Cluster at Each Resolution ---
# ****************************************************************************#
#   Re-partitions the SAME SNN graph built during integration (Step 3.3A /
#   4.2A-D) — `FindClusters()` consumes an existing graph, it never recomputes
#   nearest neighbours, so sweeping five resolutions is cheap.
#
#   COLUMN NAMING: `FindClusters()` stores each run under Seurat's internal
#   `RNA_snn_res.<res>` label. Rather than renaming all metadata columns with a
#   trailing `colnames()`/`gsub()` pass (fragile: it assumes the
#   `RNA_snn_res.` prefix and rewrites every column name), we store each
#   partition immediately under a stable, reader-friendly `clusters_res_<res>`
#   alias via `AddMetaData()`. Downstream code can then reference the columns
#   directly instead of reconstructing names from string prefixes. The count
#   is read back from the object's idents (which `FindClusters()` sets to the
#   latest partition), not from a string-reconstructed column name.
for (res in resolutions) {
    integrated_final <- FindClusters(
        integrated_final,
        resolution = res,
        verbose = FALSE
    )

    cluster_idents <- Idents(integrated_final)
    integrated_final <- AddMetaData(
        integrated_final,
        metadata = as.character(cluster_idents),
        col.name = paste0("clusters_res_", res)
    )

    cat(sprintf("Resolution %.1f: %d clusters\n", res, nlevels(cluster_idents)))
}

cat("\nCluster assignments stored in integrated_final$clusters_res_* \n")
cat("(one column per resolution; seurat_clusters now reflects res = 1.2 only)\n")


# ****************************************************************************#
# SUMMARY & PIPELINE MILESTONE TRANSITION
# ****************************************************************************#
# WHERE WE STARTED:
#   Step 5.3B committed the pipeline to one integrated object
#   (`integrated_final`), but the object carried only a single clustering from
#   a fixed, arbitrary resolution (0.6) that downstream steps would have
#   silently inherited.
#
# WHAT WE HAVE ACCOMPLISHED:
#   We re-partitioned `integrated_final`'s SNN graph at five resolutions
#   (0.4, 0.6, 0.8, 1.0, 1.2), printing the community count each resolution
#   recovers and storing every partition under an explicit, stable
#   `clusters_res_<res>` metadata column. The multi-resolution landscape now
#   lives in one object, ready to be examined without re-running clustering.
#
# WHERE WE ARE HEADING (VISUALIZE → SELECT → SAVE):
#   The guide continues with three steps built directly on these columns:
#   (1) visualize each resolution's UMAP colored by its `clusters_res_*`
#   assignment to judge coarse-vs-fine structure side by side (guide §7.5,
#   "STEP 16"); (2) compute a silhouette score per resolution — on a
#   subsample, in the integrated embedding — to objectively select the
#   optimal resolution that balances well-separated, interpretable clusters
#   (guide §7.6, "STEP 17"); and (3) assess cluster quality/stability, then
#   save the final integrated-and-clustered object (guide §7.7–7.8, §8).
#   From there the annotated, clustered data feeds cell-type annotation and
#   within-cell-type differential expression between Healthy and
#   Periodontitis_Post_Treatment.
# ****************************************************************************#