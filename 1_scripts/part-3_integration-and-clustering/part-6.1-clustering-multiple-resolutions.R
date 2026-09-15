#-----------------------------------------------
# STEP 15: Perform clustering at multiple resolutions
#-----------------------------------------------

# Test multiple resolutions
resolutions <- c(0.4, 0.6, 0.8, 1.0, 1.2)

for (res in resolutions) {
  integrated_final <- FindClusters(
    integrated_final,
    resolution = res,
    verbose = FALSE
  )

  n_clusters <- length(unique(integrated_final@meta.data[[paste0("RNA_snn_res.", res)]]))
  cat(sprintf("Resolution %.1f: %d clusters\n", res, n_clusters))
}

# Rename cluster columns for clarity
colnames(integrated_final@meta.data) <- gsub("RNA_snn_res\\.", "clusters_res_",
                                              colnames(integrated_final@meta.data))