#-----------------------------------------------
# STEP 6.3: Evaluate optimal clustering resolution
#-----------------------------------------------

# We'll use silhouette scores to assess cluster quality at each resolution
# Silhouette score measures how similar cells are to their own cluster vs other clusters
# Values range from -1 to 1 (higher is better)


# Map UMAP reduction name to corresponding integrated reduction name
# We use integrated space (not UMAP) for silhouette calculation because:
# - Integrated space preserves more information (30+ dimensions vs UMAP's 2)
# - Silhouette scores are more meaningful in the space where clustering was performed
integrated_reduction <- switch(reduction_final,
  "umap.cca"= "integrated.cca",
  "umap.rpca"= "integrated.rpca",
  "umap.harmony"= "harmony",
  "umap.mnn"= "integrated.mnn",
  "umap"= "pca"# Fallback for naive merge
)

cat("\nEvaluating clustering resolutions using", integrated_reduction, "space\n")

# For large datasets (>10,000 cells), subsample for silhouette calculation (use 5000 cells in our case)
# Computing distance matrices for all cells is computationally prohibitive
n_cells <- ncol(integrated_final)
max_cells_for_silhouette <- 5000

if (n_cells > max_cells_for_silhouette) {
  cat("Dataset has", n_cells, "cells - subsampling", max_cells_for_silhouette, "cells for silhouette calculation\n")
  set.seed(42)
  subsample_idx <- sample(1:n_cells, max_cells_for_silhouette)
} else{
  subsample_idx <- 1:n_cells
}

silhouette_scores <- sapply(resolutions, function(res) {
  cluster_col <- paste0("clusters_res_", res)

  # Get clusters for subsampled cells
  clusters <- as.numeric(integrated_final@meta.data[[cluster_col]][subsample_idx])

  # Get integrated coordinates for subsampled cells (use first 30 dimensions)
  coords <- Embeddings(integrated_final, reduction = integrated_reduction)[subsample_idx, ]
  coords <- coords[, 1:min(30, ncol(coords))]

  # Calculate silhouette (only if we have at least 2 clusters)
  if (length(unique(clusters)) > 1) {
    dist_matrix <- dist(coords)
    sil <- silhouette(clusters, dist_matrix)
    mean(sil[, 3])
  } else{
    NA# Return NA if only 1 cluster
  }
})

# Create resolution comparison table
resolution_comparison <- data.frame(
  resolution = resolutions,
  n_clusters = sapply(resolutions, function(res) {
    cluster_col <- paste0("clusters_res_", res)
    length(unique(integrated_final@meta.data[[cluster_col]]))
  }),
  silhouette_score = silhouette_scores
)

cat("\nClustering resolution comparison:\n")
print(resolution_comparison)

##   resolution n_clusters silhouette_score
## 1        0.4         17        0.2040787
## 2        0.6         23        0.2013915
## 3        0.8         24        0.1951003
## 4        1.0         26        0.1962814
## 5        1.2         30        0.1718389