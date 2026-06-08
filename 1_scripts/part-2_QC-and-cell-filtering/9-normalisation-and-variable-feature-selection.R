#-----------------------------------------------
# STEP 9: Normalization and variable features
#-----------------------------------------------

# --- 1. Log-Normalization Execution ---
# Why we normalize: Raw sequencing depth varies across single cells due to
# technical capturing efficiency, not just biology. If Cell A has 10,000 UMIs
# and Cell B has 2,000 UMIs, a gene will appear 5 times more expressed in Cell A
# purely because it was sequenced deeper.
#
# Under the hood ("LogNormalize"):
#   - It scales each cell's gene counts to a universal constant (10,000 UMIs).
#   - It then applies a natural log transformation log(count + 1). This squashes
#     the extreme values of high-flying genes so they do not dominate downstream
#     mathematical computations like principal component analyses (PCA).
seurat_obj <- NormalizeData(seurat_obj,
    normalization.method = "LogNormalize",
    scale.factor = 10000, verbose = FALSE
)

# --- 2. Highly Variable Feature Selection ---
# Why we isolate genes: Out of ~20,000 genes, the vast majority are either flat
# baseline background or household genes expressed equally across all cell types
# (like actin). Keeping all of them adds massive computational dead-weight.
#
# Under the hood ("vst"):
#   - The Variance Stabilizing Transformation (VST) calculates a mean-variance
#     relationship model across your dataset.
#   - It ranks features by their standardized variance. This isolates the top
#     2,000 genes whose expression swings wildly between different cells. These
#     genes are the true markers that define different cell populations.
seurat_obj <- FindVariableFeatures(seurat_obj,
    selection.method = "vst",
    nfeatures = 2000, verbose = FALSE
)

# --- 3. Extracting Key Lineage Drivers ---
# Pulls the top 10 most variable genes out of the object to look for dominant
# biological pathways. In immune datasets, these are typically immunoglobulin
# chains (e.g., IGLC2, IGKC, IGHM) driving cell-type discrimination.
top10 <- head(VariableFeatures(seurat_obj), 10)
cat("Top 10 variable genes:", paste(top10, collapse = ", "), "\n")

# --- 4. High-Variance Feature Visualization ---
# Generates your '07_variable_features.png' scatter plot.
# How to interpret the visual output and understand the statistical jargon:
#
#   - Average Expression (x-axis): This measures the baseline abundance of a
#     gene's transcripts across your entire dataset. Because single-cell gene
#     counts are highly skewed (a few genes have millions of reads while others
#     have zero), this axis uses an exponential log scale to stretch out the
#     data so you can see low, medium, and high-abundance genes clearly.
#
#   - Standardized Variance (y-axis): This is the core "biological signal
#     sensor." In standard statistics, highly expressed genes naturally look
#     like they vary more just because their raw numbers are larger. Seurat
#     uses a mathematical model (vst) to strip away this technical library-size
#     bias. What remains on this axis is purely "standardized variance"—a score
#     showing how aggressively a gene's volume button is turned up or down
#     between different cells, completely independent of its absolute size.
#
#   - Black Dots (The Non-Variable Background): The dense black shelf at the
#     bottom represents genes that do not change much between cell states. They
#     are either flat computational noise or crucial housekeeping genes (like
#     cellular structural components) that are turned on equally in every cell.
#
#   - Red Dots (The 2,000 Highly Variable Features): These genes actively swing
#     in expression from cell to cell. They represent the active biological
#     variance driving your dataset.
#
#   - Labeled Outliers (The Lineage Markers): The high-flying text points
#     at the very top have maximum standardized variance. They are heavily
#     turned on in one specific cell type and completely turned off in others,
#     making them the primary drivers for downstream clustering.
p10 <- VariableFeaturePlot(seurat_obj)
p10 <- LabelPoints(plot = p10, points = top10, repel = TRUE)
ggsave("plots/07_variable_features.png", p10, width = 10, height = 7, dpi = 300)
