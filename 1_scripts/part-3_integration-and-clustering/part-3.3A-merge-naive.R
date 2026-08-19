#-----------------------------------------------
# STEP 3: Naive merge without integration
#-----------------------------------------------

# Merge all samples
merged_naive <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = names(seurat_list),
  project = "GSE174609_Naive_Merge"
)

cat("Merged dataset:", ncol(merged_naive), "cells ×", nrow(merged_naive), "genes\n")
# Merged dataset: 72649 cells × 18861 genes

# Standard Seurat workflow
merged_naive <- NormalizeData(merged_naive, verbose = FALSE)
merged_naive <- FindVariableFeatures(merged_naive, nfeatures = 2000, verbose = FALSE)
merged_naive <- ScaleData(merged_naive, verbose = FALSE)
merged_naive <- RunPCA(merged_naive, npcs = 50, verbose = FALSE)
merged_naive <- RunUMAP(merged_naive, dims = 1:30, reduction = "pca", verbose = FALSE)
merged_naive <- FindNeighbors(merged_naive, dims = 1:30, verbose = FALSE)
merged_naive <- FindClusters(merged_naive, resolution = 0.6, verbose = FALSE)

cat("Clustering complete:", length(unique(merged_naive$seurat_clusters)), "clusters identified\n")
# Clustering complete: 19 clusters identified