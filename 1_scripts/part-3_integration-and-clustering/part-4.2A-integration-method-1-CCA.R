#-----------------------------------------------
# STEP 7: CCA Integration
#-----------------------------------------------

# Integrate using CCA
integrated_cca <- IntegrateLayers(
  object = merged_seurat,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca",
  dims = 1:30,
  verbose = FALSE
)

# Standard downstream workflow
integrated_cca <- FindNeighbors(integrated_cca, reduction = "integrated.cca", dims = 1:30, verbose = FALSE)
integrated_cca <- FindClusters(integrated_cca, resolution = 0.6, verbose = FALSE)
integrated_cca <- RunUMAP(integrated_cca, reduction = "integrated.cca", dims = 1:30,
                          reduction.name = "umap.cca", verbose = FALSE)

cat("CCA integration complete:", length(unique(integrated_cca$seurat_clusters)), "clusters\n")
# CCA integration complete: 23 clusters