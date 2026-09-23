# ****************************************************************************#
# STEP 6.2: Visualize clustering results across resolutions
# ****************************************************************************#


# --- PREREQUISITES (scripts that must run before this one in this session) ---
# ****************************************************************************#
#   • part-3.0-install-packages.R — packages for library() calls.
#   • part-3.1-load-libraries-and-configuration.R — RUN_ID, LOG_STEP, and
#     output-directory variables.
#   • part-6.1-clustering-multiple-resolutions.R — supplies `integrated_final`
#     with its `clusters_res_<res>` columns, and `resolutions`, the vector of
#     resolutions swept there.
#   The dependency manager is always re-sourced here (base-R only, cheap) so
#   the latest DEP_TABLE is loaded on every run; ensure_dependencies() does
#   nothing if those objects are already in the environment, otherwise it
#   offers to source them (interactive) or stops with a clear message
#   (non-interactive/batch). When everything is already satisfied, an
#   interactive session is asked whether to re-source fresh or proceed as-is.
source("1_scripts/part-3_integration-and-clustering/part-3.0-dependencies.R")
ensure_dependencies(step = "part-6.2-view-multiple-resolutions-clustering.R")

# NOTE: Requires Step 6.1 to have already run in this session — that step is
#       what writes the `clusters_res_*` columns plotted below.


# --- WHY THIS STEP EXISTS ---
# ****************************************************************************#
#   Step 6.1 swept five resolutions and printed a cluster COUNT for each
#   one, but a count alone can't tell you whether that partitioning makes
#   biological sense — "18 clusters" could mean 18 genuinely distinct cell
#   populations, or it could mean one real cell type got needlessly
#   fragmented into near-duplicate pieces while two others stayed lumped
#   together. The only way to actually judge that is to look at where the
#   cluster boundaries fall on the UMAP.
#
#   This step renders all five resolutions as UMAP panels side by side, so
#   the coarse-to-fine progression from Step 6.1 becomes something you can
#   visually inspect in one figure rather than five separate ones — and so
#   the eventual choice of resolution (guide §7.6) is based on seeing the
#   clusters, not just counting them.


# --- 1. Build One UMAP Panel Per Resolution ---
# ****************************************************************************#
#   `reduction_final` (set in Step 5.3B) is used for every panel — all five
#   resolutions were computed on the SAME integrated embedding, so they are
#   only ever different partitionings of that one space, not different
#   spaces. `label = TRUE` numbers each cluster directly on the plot, which
#   is more readable than a legend once cluster counts get into the double
#   digits, especially at the higher resolutions.
plot_list <- LOG_STEP("Building UMAP panels for each resolution...", {
  lapply(resolutions, function(res) {
    cluster_col <- paste0("clusters_res_", res)
    n_clusters <- length(unique(integrated_final@meta.data[[cluster_col]]))

    DimPlot(integrated_final,
            reduction = reduction_final,
            group.by = cluster_col,
            label = TRUE,
            label.size = 4,
            pt.size = 0.3) +
      ggtitle(paste0("Resolution ", res, " (", n_clusters, " clusters)")) +
      NoLegend()
  })
})


# --- 2. Assemble & Save the Comparison Grid ---
# ****************************************************************************#
#   `wrap_plots(ncol = 2)` lays the five panels out 2-per-row (with the
#   sixth grid cell left empty), so coarse (0.4) and fine (1.2) resolutions
#   can be compared directly by eye rather than flipping between files.
combined_resolutions <- wrap_plots(plot_list, ncol = 2)

ggsave(
    file.path(PLOTS_CLUSTERING_DIR, "06_multi_resolution_clustering.png"),
    plot = combined_resolutions, width = 14, height = 14, dpi = 300
)

cat("→ Saved:", file.path(PLOTS_CLUSTERING_DIR, "06_multi_resolution_clustering.png"), "\n")
cat("   Look for: clusters that stay stable across resolutions (real\n")
cat("   structure) vs clusters that only appear at high resolution and\n")
cat("   split off a larger neighbour (possible over-splitting).\n\n")


# ****************************************************************************#
# SUMMARY & PIPELINE MILESTONE TRANSITION
# ****************************************************************************#
# WHERE WE STARTED:
#   Step 6.1 left us with five parallel cluster assignments stored as
#   `clusters_res_<res>` columns — a coarse-to-fine landscape of options,
#   but only visible as five cluster counts, not as anything you could
#   evaluate by eye.
#
# WHAT WE HAVE ACCOMPLISHED:
#   We rendered each resolution's clustering as its own labeled UMAP panel
#   and combined all five into a single comparison grid
#   (06_multi_resolution_clustering.png). This turns the abstract
#   coarse/medium/fine distinction from Step 6.1 into something concrete:
#   you can now see which cluster boundaries persist across resolutions and
#   which ones only appear once the resolution gets pushed high.
#
# WHERE WE ARE HEADING (NEXT: DETERMINING OPTIMAL RESOLUTION):
#   A visual grid is useful for a first impression, but "resolution 0.8
#   looks about right" isn't a number you can report or defend on its own.
#   The next step computes a silhouette score for each resolution — on a
#   subsample, in the integrated embedding — to objectively quantify how
#   well-separated each partitioning's clusters actually are (guide §7.6),
#   giving this visual comparison a numeric counterpart before a final
#   resolution is chosen.
# ****************************************************************************#