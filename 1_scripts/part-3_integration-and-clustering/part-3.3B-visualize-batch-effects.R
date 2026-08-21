#-----------------------------------------------
# STEP 4: Visualize batch effects in naive merge
#-----------------------------------------------

# Define color palettes for visualization
# 8 samples need 8 distinct, visible colors
sample_colors <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",  # Healthy 1-4: Red, Blue, Green, Purple
  "#FF7F00", "#A65628", "#F781BF", "#999999"# Post_Patient 1-4: Orange, Brown, Pink, Gray
)
names(sample_colors) <- sample_metadata$sample_id

# 2 conditions need 2 distinct colors
condition_colors <- c(
  "Healthy"= "#2E86AB",       # Blue
  "Post_Treatment"= "#F18F01"# Orange
)

# UMAP colored by sample - shows batch effects
p1_naive <- DimPlot(merged_naive, reduction = "umap", group.by = "sample_id",
                     pt.size = 0.05, cols = sample_colors) +
  ggtitle("Naive Merge: Colored by Sample") +
  theme(legend.position = "right", legend.text = element_text(size = 8))

# UMAP colored by condition
p2_naive <- DimPlot(merged_naive, reduction = "umap", group.by = "condition",
                     pt.size = 0.05, cols = condition_colors) +
  ggtitle("Naive Merge: Colored by Condition")

# UMAP colored by clusters
p3_naive <- DimPlot(merged_naive, reduction = "umap", group.by = "seurat_clusters",
                     pt.size = 0.05, label = TRUE, label.size = 5) +
  ggtitle("Naive Merge: Clusters") +
  NoLegend()

# Split by sample to see separation
p4_naive <- DimPlot(merged_naive, reduction = "umap", group.by = "condition",
                     split.by = "sample_id", pt.size = 0.05, ncol = 4,
                     cols = condition_colors) +
  ggtitle("Naive Merge: Split by Sample") +
  theme(strip.text = element_text(size = 9, face = "bold"))

# Combine plots
combined_naive <- (p1_naive | p2_naive) / (p3_naive | p4_naive)
ggsave(
  file.path(PLOTS_OUT_DIR, "01_naive_merge_batch_effects.png"),
  plot = combined_naive, width = 16, height = 12, dpi = 300
)
