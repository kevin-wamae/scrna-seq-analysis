#-----------------------------------------------
# STEP 14: Select optimal integration method
#-----------------------------------------------

# Rank methods by mixing score (primary criterion)
mixing_results <- mixing_results[order(-mixing_results$mixing_score), ]

cat("\nMethod ranking (by sample mixing):\n")
for (i in 1:nrow(mixing_results)) {
  cat(sprintf("  %d. %s (mixing: %.2f, separation: %.2f)\n",
              i,
              mixing_results$method[i],
              mixing_results$mixing_score[i],
              mixing_results$condition_separation[i]))
}

# Select top method
best_method <- mixing_results$method[1]
cat("\nRecommended method:", best_method, "\n")

# For subsequent analysis, we'll use the best method
# In most cases, this will be Harmony or CCA
if (best_method == "Harmony") {
  integrated_final <- integrated_harmony
  reduction_final <- "umap.harmony"
} else if (best_method == "CCA") {
  integrated_final <- integrated_cca
  reduction_final <- "umap.cca"
} else if (best_method == "RPCA") {
  integrated_final <- integrated_rpca
  reduction_final <- "umap.rpca"
} else if (best_method == "FastMNN") {
  integrated_final <- integrated_fastmnn
  reduction_final <- "umap.mnn"
} else{
  # Default to CCA if method not found
  integrated_final <- integrated_cca
  reduction_final <- "umap.cca"
  cat("Warning: Method not recognized, defaulting to CCA\n")
}

cat("\nUsing", best_method, "for downstream clustering analysis\n")
# Using CCA for downstream clustering analysis