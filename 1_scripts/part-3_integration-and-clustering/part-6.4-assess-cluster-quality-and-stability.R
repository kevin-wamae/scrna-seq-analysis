#-----------------------------------------------
# STEP 18: Assess cluster quality and stability
#-----------------------------------------------

# Check cluster sizes
cluster_sizes <- table(integrated_final$seurat_clusters)
cat("\nCluster sizes:\n")
print(cluster_sizes)

# Identify very small clusters (<1% of cells)
min_cluster_size <- 0.01 * ncol(integrated_final)
small_clusters <- which(cluster_sizes < min_cluster_size)

if (length(small_clusters) > 0) {
  cat("\n⚠  Warning: Small clusters detected (<1% of cells):\n")
  print(cluster_sizes[small_clusters])
  cat("   These may represent rare cell types or over-clustering\n")
}

# Check sample distribution across clusters
sample_cluster_table <- table(integrated_final$seurat_clusters, integrated_final$sample_id)
sample_cluster_pct <- prop.table(sample_cluster_table, margin = 1) * 100

# Identify sample-specific clusters (>70% from one sample)
sample_specific <- apply(sample_cluster_pct, 1, max) > 70
if (any(sample_specific)) {
  cat("\n⚠  Warning: Sample-dominated clusters detected (>70% from single sample):\n")
  print(names(which(sample_specific)))
  cat("   These may indicate incomplete batch correction\n")
}

# Visualize sample distribution in clusters
sample_dist_data <- as.data.frame.matrix(sample_cluster_pct)
sample_dist_data$cluster <- rownames(sample_dist_data)
sample_dist_long <- reshape2::melt(sample_dist_data, id.vars = "cluster",
                                    variable.name = "sample", value.name = "percentage")

p_sample_dist <- ggplot(sample_dist_long, aes(x = cluster, y = percentage, fill = sample)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Sample Distribution Across Clusters",
       subtitle = "Check for sample-dominated clusters (poor integration)",
       x = "Cluster", y = "Percentage of Cells") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

ggsave(
    file.path(PLOTS_CLUSTERING_DIR, "07_sample_distribution_clusters.png"),
    plot = p_sample_dist, width = 12, height = 6, dpi = 300
)

# Check condition distribution
condition_dist <- table(integrated_final$seurat_clusters, integrated_final$condition)
condition_dist_pct <- prop.table(condition_dist, margin = 1) * 100

cat("\nCondition distribution across clusters:\n")
print(round(condition_dist_pct, 1))

# Visualize condition distribution
condition_dist_data <- as.data.frame.matrix(condition_dist_pct)
condition_dist_data$cluster <- rownames(condition_dist_data)
condition_dist_long <- reshape2::melt(condition_dist_data, id.vars = "cluster",
                                      variable.name = "condition", value.name = "percentage")

p_condition_dist <- ggplot(condition_dist_long, aes(x = cluster, y = percentage, fill = condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Condition Distribution Across Clusters",
       subtitle = "Check if biological conditions are preserved",
       x = "Cluster", y = "Percentage of Cells") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("Healthy"= "#2E86AB",
                                "Post_Treatment"= "#F18F01"))

ggsave(
    file.path(PLOTS_CLUSTERING_DIR, "08_condition_distribution_clusters.png"),
    plot = p_condition_dist, width = 12, height = 6, dpi = 300
)
