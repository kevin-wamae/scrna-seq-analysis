#-----------------------------------------------
# STEP 5: Quantify batch effects with mixing metrics
#-----------------------------------------------

# Calculate per-cluster sample composition
cluster_composition <- table(merged_naive$seurat_clusters, merged_naive$sample_id)
cluster_composition_pct <- prop.table(cluster_composition, margin = 1) * 100

# Find clusters dominated by single samples (>50% from one sample = batch-driven)
dominant_sample_clusters <- apply(cluster_composition_pct, 1, max) > 50
n_batch_clusters <- sum(dominant_sample_clusters)

cat("\nBatch effect assessment:\n")
cat("  Clusters dominated by single sample (>50%):", n_batch_clusters,
    "/", nrow(cluster_composition), "\n")
# Clusters dominated by single sample (>50%): 1 / 19

if (n_batch_clusters > 0) {
  cat("  ⚠  Naive merge shows significant batch effects\n")
  cat("  → Integration is necessary\n")
}
# → Integration is necessary