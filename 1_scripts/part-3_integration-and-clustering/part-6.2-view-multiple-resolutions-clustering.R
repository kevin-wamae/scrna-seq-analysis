#-----------------------------------------------
# STEP 6.2: Visualize clustering results across resolutions
#-----------------------------------------------

# Create UMAP plots for each resolution
plot_list <- lapply(resolutions, function(res) {
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

# Combine plots
combined_resolutions <- wrap_plots(plot_list, ncol = 2)


       # Save comparison by condition image
ggsave(
    file.path(PLOTS_CLUSTERING_DIR, "06_multi_resolution_clustering.png"),
    plot = combined_resolutions, width = 14, height = 14, dpi = 300
)
