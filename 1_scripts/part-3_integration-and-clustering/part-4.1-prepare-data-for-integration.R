#-----------------------------------------------
# STEP 4.1: Prepare data for integration
#-----------------------------------------------

# Merge all samples with split layers (required for Seurat 5 integration)
merged_seurat <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = names(seurat_list),
  project = "GSE174609_Integration"
)

# Add all metadata
# Extract sample_id from cell names (added via add.cell.ids)
# Cell names format: "SampleID_BARCODE-1"
cell_names <- colnames(merged_seurat)
merged_seurat$sample_id <- gsub("_[ACGT].*$", "", cell_names)  # Remove barcode, keep sample_id

# Add other metadata by matching sample_id
merged_seurat$condition <- sample_metadata$condition[match(merged_seurat$sample_id, sample_metadata$sample_id)]
merged_seurat$patient_id <- sample_metadata$patient_id[match(merged_seurat$sample_id, sample_metadata$sample_id)]

# Check if layers are already split
current_layers <- Layers(merged_seurat[["RNA"]])
if (length(current_layers) > 1) {
  cat("Layers already split by merge operation (", length(current_layers), " layers)\n", sep = "")
} else{
  # Split layers by sample if not already split
  cat("Splitting layers by sample\n")
  merged_seurat[["RNA"]] <- split(merged_seurat[["RNA"]], f = merged_seurat$sample_id)
}

# Normalize and find variable features on split layers
merged_seurat <- NormalizeData(merged_seurat, verbose = FALSE)
merged_seurat <- FindVariableFeatures(merged_seurat, nfeatures = 2000, verbose = FALSE)

# Scale data and run PCA on split layers
merged_seurat <- ScaleData(merged_seurat, verbose = FALSE)
merged_seurat <- RunPCA(merged_seurat, npcs = 50, verbose = FALSE)

cat("Data prepared for integration\n")
cat("Layers:", length(Layers(merged_seurat, search = "data")), "samples\n")
# Layers: 8 samples