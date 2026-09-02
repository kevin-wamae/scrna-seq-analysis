#-----------------------------------------------
# STEP 11: Generate comparison UMAPs
#-----------------------------------------------

cat("\n=== Generating Comparison Plots ===\n")

# Create a function to make standardized UMAP plots
make_comparison_plot <- function(seurat_obj, reduction_name, title, group_by = "sample_id", colors = NULL) {
  p <- DimPlot(seurat_obj, reduction = reduction_name, group.by = group_by, pt.size = 0.05) +
    ggtitle(title) +
    theme(plot.title = element_text(face = "bold", size = 12),
          legend.text = element_text(size = 7),
          legend.key.size = unit(0.3, "cm"))

  # Apply custom colors if provided
  if (!is.null(colors)) {
    p <- p + scale_color_manual(values = colors)
  }

  return(p)
}

# Read the sample metadata
sample_metadata <- read.delim("2_input/sample-metadata/sample_names.tsv",
    stringsAsFactors = FALSE) %>%
    select(sample_id, condition, patient_id) %>%
    filter(condition != "Periodontitis_Pre_Treatment")

# Sanity check that the metadata looks right
cat("Sample metadata loaded:\n")
cat("  Samples:   ", nrow(sample_metadata), "\n")
cat("  Conditions:", paste(unique(sample_metadata$condition), collapse = ", "), "\n\n")


# Define color palettes for visualization
# 8 samples need 8 distinct, visible colors
sample_colors <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",  # Healthy 1-4: Red, Blue, Green, Purple
  "#FF7F00", "#A65628", "#F781BF", "#999999"   # Post_Patient 1-4: Orange, Brown, Pink, Gray
)
names(sample_colors) <- sample_metadata$sample_id

# 2 conditions need 2 distinct colors
condition_colors <- c(
  "Healthy"= "#2E86AB",       # Blue
  "Post_Treatment"= "#F18F01"# Orange
)

# Load processed objects from checkpoint files
merged_naive <- qs2::qs_read(file.path(
    DATA_CHECKPOINT_DIR,
    "01_merged_naive.qs2"
))
integrated_cca <- qs2::qs_read(file.path(
    DATA_CHECKPOINT_DIR,
    "02_integrated_cca.qs2"
))
integrated_rpca <- qs2::qs_read(file.path(
    DATA_CHECKPOINT_DIR,
    "03_integrated_rpca.qs2"
))
integrated_harmony <- qs2::qs_read(file.path(
    DATA_CHECKPOINT_DIR,
    "04_integrated_harmony.qs2"
))
integrated_fastmnn <- qs2::qs_read(file.path(
    DATA_CHECKPOINT_DIR,
    "05_integrated_fastmnn.qs2"
))

# UMAPs colored by sample (assesses mixing)
p_sample_naive <- make_comparison_plot(
    merged_naive,
    "umap",
    "Naive Merge",
    colors = sample_colors
)
p_sample_cca <- make_comparison_plot(
    integrated_cca,
    "umap.cca",
    "CCA",
    colors = sample_colors
)
p_sample_rpca <- make_comparison_plot(
    integrated_rpca,
    "umap.rpca",
    "RPCA",
    colors = sample_colors
)
p_sample_harmony <- make_comparison_plot(
    integrated_harmony,
    "umap.harmony",
    "Harmony",
    colors = sample_colors
)
p_sample_mnn <- make_comparison_plot(
    integrated_fastmnn,
    "umap.mnn",
    "FastMNN",
    colors = sample_colors
)

# Combine sample-colored UMAPs
combined_samples <- (p_sample_naive | p_sample_cca | p_sample_rpca) /
                    (p_sample_harmony | p_sample_mnn | plot_spacer())

ggsave(
  file.path(PLOTS_OUT_DIR, "integration_comparison/02_integration_by_sample.png"),
  plot = combined_samples, width = 16, height = 12, dpi = 300
)

# UMAPs colored by condition (assesses biological preservation)
p_cond_naive <- make_comparison_plot(merged_naive, "umap", "Naive Merge", "condition",
                                     colors = condition_colors)
p_cond_cca <- make_comparison_plot(integrated_cca, "umap.cca", "CCA", "condition",
                                   colors = condition_colors)
p_cond_rpca <- make_comparison_plot(integrated_rpca, "umap.rpca", "RPCA", "condition",
                                    colors = condition_colors)
p_cond_harmony <- make_comparison_plot(integrated_harmony, "umap.harmony", "Harmony", "condition",
                                       colors = condition_colors)
p_cond_mnn <- make_comparison_plot(integrated_fastmnn, "umap.mnn", "FastMNN", "condition",
                                   colors = condition_colors)

# Combine condition-colored UMAPs
combined_conditions <- (p_cond_naive | p_cond_cca | p_cond_rpca) /
                       (p_cond_harmony | p_cond_mnn | plot_spacer())

ggsave(
  file.path(PLOTS_OUT_DIR, "integration_comparison/03_integration_by_condition.png"),
  plot = combined_conditions, width = 16, height = 12, dpi = 300
)
