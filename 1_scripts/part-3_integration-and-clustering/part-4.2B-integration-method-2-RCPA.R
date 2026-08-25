#-----------------------------------------------
# STEP 8: RPCA Integration
#-----------------------------------------------

# Integrate using RPCA
integrated_rpca <- IntegrateLayers(
  object = merged_seurat,
  method = RPCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.rpca",
  dims = 1:30,
  verbose = FALSE
)

# Downstream workflow
integrated_rpca <- FindNeighbors(integrated_rpca, reduction = "integrated.rpca", dims = 1:30, verbose = FALSE)
integrated_rpca <- FindClusters(integrated_rpca, resolution = 0.6, verbose = FALSE)
integrated_rpca <- RunUMAP(integrated_rpca, reduction = "integrated.rpca", dims = 1:30,
                           reduction.name = "umap.rpca", verbose = FALSE)

cat("RPCA integration complete:", length(unique(integrated_rpca$seurat_clusters)), "clusters\n")
# RPCA integration complete: 22 clusters