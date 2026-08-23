#-----------------------------------------------
# STEP 4.1: Prepare data for integration
#-----------------------------------------------
# WHY THIS STEP EXISTS:
#   Section 3 established that our naive merge has batch-driven clustering
#   (1/19 clusters dominated by a single sample). Before we can run any
#   integration method (CCA, RPCA, Harmony, FastMNN), Seurat 5's integration
#   workflow needs the data in a specific shape: one merged object, but with
#   each sample's counts kept in SEPARATE layers rather than pooled into one.
#   Integration methods learn a correction by comparing samples layer-by-
#   layer — if the layers were pre-merged into one pooled layer, there'd be
#   nothing for the algorithm to compare.

# --- 1. Merge all samples into a single object ---
# ****************************************************************************#
#   `add.cell.ids` prefixes every cell barcode with its sample name (e.g.
#   "Healthy_1_AAACCCAAGAAACCAT-1"), which prevents barcode collisions across
#   samples and is also what lets us recover sample_id below by parsing the
#   cell name. This merge does NOT pool counts into one matrix — Seurat 5
#   keeps each sample's counts as its own layer inside the RNA assay (more on
#   this in step 3).
merged_seurat <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = names(seurat_list),
  project = "GSE174609_Integration"
)

# --- 2. Recover sample-level metadata from cell names ---
# ****************************************************************************#
#   `merge()` does not carry over our external metadata table automatically,
#   so we rebuild it here. Cell names now look like "SampleID_BARCODE-1"
#   (e.g. "Healthy_1_AAACCCAAGAAACCAT-1") thanks to `add.cell.ids` above.
#   The regex strips everything from the first nucleotide-code barcode
#   character onward, leaving just the sample_id prefix.
#   CAVEAT: this regex assumes sample IDs never themselves start with a run
#   of A/C/G/T characters directly followed by more barcode-looking text —
#   true for names like "Healthy_1", but worth re-checking if sample IDs
#   ever change format.
cell_names <- colnames(merged_seurat)
merged_seurat$sample_id <- gsub("_[ACGT].*$", "", cell_names)  # strip barcode, keep sample_id

# With sample_id recovered, join the rest of the per-sample metadata
# (condition, patient_id) from the external sample sheet via a lookup match.
merged_seurat$condition   <- sample_metadata$condition[match(merged_seurat$sample_id, sample_metadata$sample_id)]
merged_seurat$patient_id  <- sample_metadata$patient_id[match(merged_seurat$sample_id, sample_metadata$sample_id)]

# --- 3. Ensure RNA assay layers are split by sample ---
# ****************************************************************************#
#   Seurat 5 integration methods operate on a MULTI-LAYER assay — one layer
#   of counts per sample — rather than one pooled layer. Depending on Seurat
#   version/settings, `merge()` may already leave the assay split; this check
#   avoids redundantly re-splitting (and the associated recompute cost) if
#   it's already in the right shape.
current_layers <- Layers(merged_seurat[["RNA"]])

if (length(current_layers) > 1) {
  cat("Layers already split by merge operation (", length(current_layers), " layers)\n", sep = "")
} else {
  cat("Splitting layers by sample\n")
  merged_seurat[["RNA"]] <- split(merged_seurat[["RNA"]], f = merged_seurat$sample_id)
}

# --- 4. Normalize, find variable features, scale, and run PCA — per layer ---
# ****************************************************************************#
#   Run on the SPLIT layers (not pooled data): each of these steps is
#   computed independently within each sample's layer, so normalization and
#   HVG selection aren't influenced by between-sample technical differences
#   before integration has a chance to correct for them. This mirrors the
#   naive-merge workflow from Part 3, but the layer-splitting above is what
#   makes it integration-ready rather than just another naive merge.
merged_seurat <- LOG_STEP("Preparing data for integration (normalize, HVGs, scale, PCA)...", {
  merged_seurat %>%
    NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(nfeatures = 2000, verbose = FALSE) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(npcs = 50, verbose = FALSE)
})

cat("Data prepared for integration\n")
cat("Layers:", length(Layers(merged_seurat, search = "data")), "samples\n")
# Layers: 8 samples