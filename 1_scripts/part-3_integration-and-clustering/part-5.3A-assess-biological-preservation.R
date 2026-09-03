#-----------------------------------------------
# STEP 13: Assess biological signal preservation
#-----------------------------------------------

# Function to test if conditions are still distinguishable after integration
# Uses UMAP distances between conditions
assess_condition_separation <- function(seurat_obj, reduction) {
  umap_coords <- seurat_obj[[reduction]]@cell.embeddings[, 1:2]
  conditions <- seurat_obj$condition

  # Calculate mean UMAP position for each condition
  condition_centers <- aggregate(umap_coords, by = list(conditions), FUN = mean)

  # Calculate pairwise distances between condition centers
  dist_matrix <- dist(condition_centers[, -1])
  mean_dist <- mean(dist_matrix)

  return(mean_dist)
}

# Calculate for all methods
separation_results <- data.frame(
  method = names(methods_list),
  condition_separation = sapply(methods_list, function(x) {
    assess_condition_separation(x$obj, x$reduction)
  })
)

mixing_results$condition_separation <- separation_results$condition_separation

cat("\nBiological Signal Preservation:\n")
cat("(Higher separation = conditions remain distinguishable)\n\n")
print(mixing_results[, c("method", "mixing_score", "condition_separation", "n_clusters")])

# Ideal integration: High mixing + Moderate separation
# Over-integration: High mixing + Low separation (conditions lost)
# Under-integration: Low mixing + High separation (batches not corrected)

# Visualize the trade-off
p_tradeoff <- ggplot(mixing_results, aes(x = mixing_score, y = condition_separation, label = method)) +
  geom_point(size = 4, aes(color = method)) +
  geom_text_repel(size = 4, box.padding = 0.5) +
  labs(title = "Integration Quality: Mixing vs Biological Preservation",
       subtitle = "Ideal: Upper right (high mixing + preserved biology)",
       x = "Sample Mixing Score (higher = better integration)",
       y = "Condition Separation (higher = preserved biology)") +
  theme(legend.position = "none") +
  scale_color_brewer(palette = "Set2")

ggsave("3_output/2026_06_09_brown_job_3058993/integration_and_clustering/plots/integration_comparison/05_mixing_vs_preservation.png", p_tradeoff,
       width = 10, height = 7, dpi = 300)