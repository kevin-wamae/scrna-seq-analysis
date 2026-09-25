# ****************************************************************************#
# STEP 6.5: Create comprehensive visualization of final results
# ****************************************************************************#


# --- PREREQUISITES (scripts that must run before this one in this session) ---
# ****************************************************************************#
#   • part-3.0-install-packages.R — packages for library() calls.
#   • part-3.1-load-libraries-and-configuration.R — RUN_ID, PLOTS_CLUSTERING_DIR,
#     and the shared `sample_colors` / `condition_colors` palettes used below.
#   • part-5.3B-select-best-integration-method.R — supplies `best_method`, the
#     integration method that produced the embedding, named in the panel title.
#   • part-6.3-load-clustering-checkpoints.R — supplies `integrated_final` with
#     `seurat_clusters` committed to `optimal_resolution`, plus `reduction_final`.
#   The dependency manager is always re-sourced here (base-R only, cheap) so
#   the latest DEP_TABLE is loaded on every run; ensure_dependencies() does
#   nothing if those objects are already in the environment, otherwise it
#   offers to source them (interactive) or stops with a clear message
#   (non-interactive/batch). When everything is already satisfied, an
#   interactive session is asked whether to re-source fresh or proceed as-is.
source("1_scripts/part-3_integration-and-clustering/part-3.0-dependencies.R")
ensure_dependencies(step = "part-6.5-final-results-comprehensive-visualization.R")

# NOTE: Requires Step 6.3 to have loaded the clustering checkpoints in this
#       session (that load is what supplies `integrated_final` with
#       `seurat_clusters` set), and Step 5.3B for `best_method`. Both are wired
#       into DEP_TABLE, so this script no longer needs Step 6.4 to have run
#       first just to pick up `best_method`.


# --- WHY THIS STEP EXISTS ---
# ****************************************************************************#
#   Step 6.4 vetted the chosen partition numerically (cluster sizes, sample and
#   condition composition). This step is the visual bookend to that work: a
#   single figure presenting the final integrated-and-clustered object the way
#   a reader first meets it (guide §7.8, STEP 19). Four panels, each answering
#   a different question:
#
#     1. FINAL CLUSTERS (top left): what the chosen `optimal_resolution`
#        partition actually looks like on the integrated UMAP. Labels are
#        Seurat cluster IDs only — no cell-type claims are made here;
#        annotation is a downstream task.
#     2. BIOLOGICAL CONDITIONS (top right): the two conditions coloured on the
#        same embedding. Some condition-specific structure SHOULD survive
#        integration — Healthy vs Periodontitis_Post_Treatment is the
#        biological axis this study is about. What we do not want is a single
#        condition owning a whole cluster outright.
#     3. SAMPLE MIXING (bottom left): cells coloured by `sample_id`. This is
#        the integration-quality readout — a well-corrected embedding shows the
#        8 samples interleaved within each cluster, not segregated into
#        sample-shaped islands.
#     4. CLUSTERS BY CONDITION (separate figure): the same clustered UMAP
#        faceted into one panel per condition, so a cluster appearing in only
#        one condition is immediately visible.
#
#   INTERPRETATION (guide §9.1):
#     - Good integration : samples intermingle; condition signal present but
#       not driving the global structure.
#     - Red flags        : sample-shaped islands (under-integration), or a
#       perfectly mixed condition layer everywhere (over-integration).
#   These are eyeball checks that complement, not replace, the numbers in
#   `cluster_sample_composition.csv` / `cluster_condition_composition.csv`
#   from Step 6.4.
#
#   CHECKPOINT DECISION: none. This is a plotting-only step — its entire output
#   is two PNGs, and it creates no heavy intermediate that a later step would
#   otherwise recompute. Persisting a checkpoint here would add a file with no
#   consumer.


# --- 1. Main UMAP: Clusters on the Final Embedding ---
# ****************************************************************************#
#   `group.by = "seurat_clusters"` uses the partition Step 6.2B committed to
#   disk (and Step 6.3 restored), not whatever Idents state the object happens
#   to carry, so the labels are unambiguous. The title records `best_method`
#   so the figure is self-describing when read on its own.
p_final_clusters <- DimPlot(integrated_final,
                            reduction = reduction_final,
                            group.by = "seurat_clusters",
                            label = TRUE,
                            label.size = 5,
                            pt.size = 0.05) +
  ggtitle(paste0("Final Clustering (", best_method, " Integration)")) +
  theme(plot.title = element_text(face = "bold", size = 14))


# --- 2. Clusters Split by Condition ---
# ****************************************************************************#
#   `split.by = "condition"` facets the SAME cluster layout into one panel per
#   condition (ncol = 2), keeping cluster colours consistent across panels. A
#   cluster whose cells appear in one facet but not the other is
#   condition-restricted — expected for some biology, and worth explaining if
#   it dominates.
p_clusters_by_condition_split <- DimPlot(integrated_final,
                                         reduction = reduction_final,
                                         group.by = "seurat_clusters",
                                         split.by = "condition",
                                         label = TRUE,
                                         label.size = 4,
                                         pt.size = 0.05,
                                         ncol = 2) +
  ggtitle("Clusters by Treatment Condition") +
  theme(strip.text = element_text(face = "bold"))


# --- 3. Sample Mixing (Integration Quality) ---
# ****************************************************************************#
#   Coloured with the shared `sample_colors` from Step 3.1, so each sample
#   keeps the same colour it had in Steps 3.3B and 5.2A. Read this panel as the
#   visual counterpart to the >70% single-sample flag in Step 6.4.
p_sample_mixing <- DimPlot(integrated_final,
                           reduction = reduction_final,
                           group.by = "sample_id",
                           pt.size = 0.05,
                           cols = sample_colors) +
  ggtitle("Sample Mixing (Integration Quality)") +
  theme(legend.text = element_text(size = 7))


# --- 4. Biological Conditions ---
# ****************************************************************************#
#   The same embedding coloured by `condition`, using the shared
#   `condition_colors` (identical palette to Steps 3.3B and 5.2A), so the
#   biological signal can be judged at a glance without the cluster grid.
p_conditions <- DimPlot(integrated_final,
                        reduction = reduction_final,
                        group.by = "condition",
                        pt.size = 0.05,
                        cols = condition_colors) +
  ggtitle("Treatment Conditions")


# --- 5. Assemble & Save the Final Figures ---
# ****************************************************************************#
#   `patchwork` assembles the four panels: clusters and conditions on top,
#   sample mixing bottom-left with a spacer bottom-right. Two saves:
#   the combined overview (09) and the larger condition-split view (10).
#   The run-order numbering continues from Step 6.4's plots 07/08, so the
#   clustering figures read in pipeline order.
p_final_integrated_clustered <- (p_final_clusters | p_conditions) /
                                (p_sample_mixing | plot_spacer())

LOG_STEP("Saving final integrated visualization (09)...", {
    ggsave(
        file.path(PLOTS_CLUSTERING_DIR, "09_final_integrated_clustered.png"),
        plot = p_final_integrated_clustered, width = 16, height = 12, dpi = 300
    )
})
cat("→ Saved:", file.path(PLOTS_CLUSTERING_DIR, "09_final_integrated_clustered.png"), "\n")

LOG_STEP("Saving clusters-by-condition split (10)...", {
    ggsave(
        file.path(PLOTS_CLUSTERING_DIR, "10_clusters_by_condition_split.png"),
        plot = p_clusters_by_condition_split, width = 18, height = 6, dpi = 300
    )
})
cat("→ Saved:", file.path(PLOTS_CLUSTERING_DIR, "10_clusters_by_condition_split.png"), "\n")


# ****************************************************************************#
# SUMMARY & PIPELINE MILESTONE TRANSITION
# ****************************************************************************#
# WHERE WE STARTED:
#   Step 6.4 had sanity-checked the chosen partition in numbers — cluster
#   sizes and sample/condition composition — but the final integrated object
#   still had no single figure a reader could look at to see the result.
#
# WHAT WE HAVE ACCOMPLISHED:
#   We built the four-panel closing visualization of the final
#   `integrated_final` object at `optimal_resolution`: clusters on the
#   integrated UMAP, conditions on the same embedding, sample mixing as the
#   integration-quality readout, and a per-condition faceted cluster view.
#   The combined overview is saved as 09_final_integrated_clustered.png and
#   the condition split as 10_clusters_by_condition_split.png, both coloured
#   with the shared palettes so they match every earlier figure. No checkpoint
#   was written — this step only produces PNGs.
#
# WHERE WE ARE HEADING (NEXT: SAVE THE FINAL OBJECT):
#   Step 6.6 produces the actual Part 3 deliverable: it joins the assay layers
#   for a portable object, saves it to `integrated_data/` (honouring
#   CHECKPOINT_FORMAT), and writes the per-cell metadata and the one-row
#   integration summary to `metadata/` (guide §8, STEP 20), so downstream
#   cell-type annotation can start from a clean file.
# ****************************************************************************#
