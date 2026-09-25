#-----------------------------------------------
# STEP 6.5: Create comprehensive visualization of final results
#-----------------------------------------------
# TODO: `best_method` seems to be loaded in 5.3B, so running 6.4 before 6.5
#       doesn't pick this up and maybe more, needs dependency-fixing to make this optimal

# Main UMAP with clusters
p_final_clusters <- DimPlot(integrated_final,
                            reduction = reduction_final,
                            group.by = "seurat_clusters",
                            label = TRUE,
                            label.size = 5,
                            pt.size = 0.05) +
  ggtitle(paste0("Final Clustering (", best_method, " Integration)")) +
  theme(plot.title = element_text(face = "bold", size = 14))

# Split by condition to show preservation
p_by_condition <- DimPlot(integrated_final,
                          reduction = reduction_final,
                          group.by = "seurat_clusters",
                          split.by = "condition",
                          label = TRUE,
                          label.size = 4,
                          pt.size = 0.05,
                          ncol = 2) +
  ggtitle("Clusters by Treatment Condition") +
  theme(strip.text = element_text(face = "bold"))

# Colored by sample (integration quality check)
p_by_sample <- DimPlot(integrated_final,
                       reduction = reduction_final,
                       group.by = "sample_id",
                       pt.size = 0.05,
                       cols = sample_colors) +
  ggtitle("Sample Mixing (Integration Quality)") +
  theme(legend.text = element_text(size = 7))

# Colored by condition
p_by_condition_single <- DimPlot(integrated_final,
                                 reduction = reduction_final,
                                 group.by = "condition",
                                 pt.size = 0.05,
                                 cols = condition_colors) +
  ggtitle("Treatment Conditions")

# Combine final visualizations
combined_final <- (p_final_clusters | p_by_condition_single) /
                  (p_by_sample | plot_spacer())


ggsave(
    file.path(PLOTS_CLUSTERING_DIR, "09_final_integrated_clustered.png"),
    plot = combined_final, width = 16, height = 12, dpi = 300
)

ggsave(
    file.path(PLOTS_CLUSTERING_DIR, "10_clusters_by_condition_split.png"),
    plot = p_by_condition, width = 18, height = 6, dpi = 300
)
