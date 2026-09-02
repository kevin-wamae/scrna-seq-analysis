#-----------------------------------------------
# STEP 12: Calculate integration quality metrics
#-----------------------------------------------

# Function to calculate mixing metric (local inverse Simpson's Index)
# This measures sample diversity in each cell's k-nearest neighborhood
calculate_mixing_metric <- function(seurat_obj, group_by = "sample_id", reduction = "umap", k = 50) {
  # Get embedding
  embedding <- Embeddings(seurat_obj, reduction = reduction)

  # Calculate k-nearest neighbors using FNN package
  nn_result <- FNN::get.knn(embedding[, 1:2], k = k)

  # For each cell, calculate diversity of samples in its neighborhood
  group_labels <- seurat_obj@meta.data[[group_by]]
  mixing_scores <- sapply(1:nrow(nn_result$nn.index), function(i) {
    neighbors <- nn_result$nn.index[i, ]
    neighbor_groups <- group_labels[neighbors]
    props <- table(neighbor_groups) / length(neighbor_groups)
    # Inverse Simpson's Index (higher = more mixed)
    1 / sum(props^2)
  })

  return(mean(mixing_scores))
}

# Calculate for all methods
methods_list <- list(
  "Naive"= list(obj = merged_naive, reduction = "umap"),
  "CCA"= list(obj = integrated_cca, reduction = "umap.cca"),
  "RPCA"= list(obj = integrated_rpca, reduction = "umap.rpca"),
  "Harmony"= list(obj = integrated_harmony, reduction = "umap.harmony"),
  "FastMNN"= list(obj = integrated_fastmnn, reduction = "umap.mnn")
)

# Calculate mixing metrics
mixing_results <- data.frame(
  method = names(methods_list),
  mixing_score = sapply(methods_list, function(x) {
    calculate_mixing_metric(x$obj, group_by = "sample_id", reduction = x$reduction)
  }),
  n_clusters = sapply(methods_list, function(x) {
    length(unique(x$obj$seurat_clusters))
  })
)

# Display results
cat("\nIntegration Quality Metrics:\n")
cat("(Higher mixing score = better sample mixing)\n\n")
print(mixing_results)

# Save metrics
write.csv(mixing_results, "3_output/2026_06_09_brown_job_3058993/integration_and_clustering/metadata/integration_comparison_metrics.csv", row.names = FALSE)

# Visualize mixing scores
p_mixing <- ggplot(mixing_results, aes(x = reorder(method, -mixing_score), y = mixing_score, fill = method)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = round(mixing_score, 2)), vjust = -0.5) +
  labs(title = "Integration Method Comparison: Sample Mixing",
       subtitle = "Higher score = better integration (samples mix within cell types)",
       x = "Integration Method",
       y = "Mixing Score (Local Inverse Simpson's Index)") +
  theme(legend.position = "none") +
  scale_fill_brewer(palette = "Set2")

ggsave("3_output/2026_06_09_brown_job_3058993/integration_and_clustering/plots/integration_comparison/04_mixing_scores.png", p_mixing,
       width = 10, height = 6, dpi = 300)
