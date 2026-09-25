#-----------------------------------------------
# STEP 20: Save integrated and clustered data
#-----------------------------------------------

cat("\n=== Saving Final Results ===\n")

# Join layers for downstream analysis
integrated_final <- JoinLayers(integrated_final)

# Save the final integrated object
qs2::qs_save(integrated_final, file.path(DATA_OUT_DIR, "integrated_clustered_seurat.qs2"), nthreads = N_WORKERS)

# Save metadata
write.csv(integrated_final@meta.data,
          file.path(METADATA_OUT_DIR, "integrated_cell_metadata.csv"),
          row.names = TRUE)

# Create comprehensive summary
integration_summary <- data.frame(
  dataset = "GSE174609",
  n_samples = length(unique(integrated_final$sample_id)),
  n_cells_total = ncol(integrated_final),
  n_genes = nrow(integrated_final),
  integration_method = best_method,
  optimal_resolution = optimal_resolution,
  n_clusters = length(unique(integrated_final$seurat_clusters)),
  mixing_score = mixing_results$mixing_score[mixing_results$method == best_method],
  condition_separation = mixing_results$condition_separation[mixing_results$method == best_method]
)

write.csv(integration_summary, file.path(METADATA_OUT_DIR, "integration_summary.csv"), row.names = FALSE)
