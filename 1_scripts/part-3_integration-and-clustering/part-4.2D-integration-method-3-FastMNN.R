#-----------------------------------------------
# STEP 10: FastMNN Integration
#-----------------------------------------------

# Integrate using FastMNN (works on split layers)
integrated_fastmnn <- IntegrateLayers(
  object = merged_seurat,
  method = FastMNNIntegration,
  new.reduction = "integrated.mnn",
  verbose = FALSE
)

cat("FastMNN integration complete\n")

# Downstream workflow
integrated_fastmnn <- FindNeighbors(integrated_fastmnn, reduction = "integrated.mnn", dims = 1:30, verbose = FALSE)
integrated_fastmnn <- FindClusters(integrated_fastmnn, resolution = 0.6, verbose = FALSE)
integrated_fastmnn <- RunUMAP(integrated_fastmnn, reduction = "integrated.mnn", dims = 1:30,
                              reduction.name = "umap.mnn", verbose = FALSE)

cat("Clustering complete:", length(unique(integrated_fastmnn$seurat_clusters)), "clusters\n")
# FastMNN integration complete: 20 clusters