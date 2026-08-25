#-----------------------------------------------
# STEP 9: Harmony Integration
#-----------------------------------------------

# Harmony works directly on PCA embedding
# First, join layers and run standard workflow
merged_harmony <- JoinLayers(merged_seurat)
merged_harmony <- RunPCA(merged_harmony, npcs = 50, verbose = FALSE)

# Run Harmony (corrects batch effects on PCA embedding)
integrated_harmony <- RunHarmony(merged_harmony, "sample_id")

# Downstream workflow
integrated_harmony <- FindNeighbors(integrated_harmony, reduction = "harmony", dims = 1:30, verbose = FALSE)
integrated_harmony <- FindClusters(integrated_harmony, resolution = 0.6, verbose = FALSE)
integrated_harmony <- RunUMAP(integrated_harmony, reduction = "harmony", dims = 1:30,
                                reduction.name = "umap.harmony", verbose = FALSE)

cat("Harmony integration complete:", length(unique(integrated_harmony$seurat_clusters)), "clusters\n")
# Harmony integration complete: 21 clusters