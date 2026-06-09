# ****************************************************************************#
# STEP 4: Empty droplet detection
# ****************************************************************************#


# --- THE EMPTY DROPLET PROBLEM ---
# ****************************************************************************#
# Not every droplet captured by a microfluidic system contains a real cell.
# A raw single-cell run typically generates three distinct types of droplets:
#   1. Cell-containing droplets: High UMI counts, diverse gene expression.
#   2. Empty droplets: Low UMI counts, containing only ambient background RNA.
#   3. Damaged cells: Intermediate UMI counts (which we filter later in Step 7).
#
# WHY WE USE EMPTYDROPS INSTEAD OF A HARD CUTOFF:
# Simple UMI thresholds (e.g., discarding anything <1000 UMIs) bias against
# small cells. Highly active cells (macrophages) yield high counts, but tiny,
# quiescent cells (resting T-cells, stem cells) naturally yield very few
# transcripts. Simple thresholds permanently kill this rare biology.
#
# STATISTICAL DETECTION MECHANICS (PROFILE VS. VOLUME):
# EmptyDrops avoids threshold bias by testing the expression profile (gene mix)
# instead of total count numbers. It:
#   1. Estimates the ambient RNA profile ("the soup") from low-count droplets.
#   2. Tests each droplet: "Is this expression profile different from the soup?"
#   3. Assigns a False Discovery Rate (FDR) p-value to each tested barcode.
#
# CONCEPTUAL PROFILE EXAMPLE:
# Imagine two different droplets that both contain exactly 500 UMIs:
#   - Droplet A: Expresses 50 genes in a highly specific pattern (e.g., marker
#     genes matching a rare cell type). EmptyDrops identifies this as a CELL.
#   - Droplet B: Expresses genes randomly matching standard background ambient
#     housekeeping noise. EmptyDrops identifies this as an EMPTY DROPLET.

# --- WORKFLOW CONDITIONAL PREREQUISITE ---
# NOTE: If you chose to load a Cell Ranger FILTERED matrix instead of a RAW
# matrix in Step 2, you MUST skip this step (Step 4) and the Ambient RNA
# correction (Step 5). Cell Ranger's internal script has already completed a
# basic empty droplet evacuation. Skip directly to Step 6 (Doublet Detection).


# Convert Seurat object to SingleCellExperiment for EmptyDrops processing
sce <- as.SingleCellExperiment(seurat_obj)

# Track initial benchmarks for downstream pipeline tracking
initial_droplet_count <- ncol(sce)
initial_gene_count <- nrow(sce)


# --- Run EmptyDrops ---
# - `lower = 100`: The empirical "sweet spot" for 10x chemistry. Barcodes with
#   <100 UMIs are assumed to be cell-free empty droplets. They are pooled to
#   model the exact baseline composition of the ambient background soup.
# - `test.ambient = TRUE`: Forces the algorithm to evaluate droplets below the
#   lower threshold anyway. If a droplet with 80 UMIs has a specific gene mix
#   that differs drastically from the soup, it can be rescued as a valid cell.
set.seed(100)
empty_results <- emptyDrops(
  m = counts(sce),
  lower = 100,
  niters = 10000,
  test.ambient = TRUE
)


# --- STATISTICAL CLASSIFICATION ---
# ****************************************************************************#
# 1. SET AN ERROR CEILING: Evaluate the False Discovery Rate (FDR). An FDR
#    under 0.01 means there is less than a 1% chance the droplet is just soup.
#    If it passes this strict math test, it gets flagged as TRUE (a cell).
is_cell <- empty_results$FDR < 0.01

# 2. CLEAN UP SKIPPED DROPLETS: Barcodes with ultra-low counts (e.g., <100 UMIs)
#    are skipped by the algorithm to save computing power, returning an NA.
#    We explicitly turn these blanks into FALSE so they don't crash our code.
is_cell[is.na(is_cell)] <- FALSE

# 3. GRAB THE WHITELIST: Look at the matrix column names and pull out the
#    exact text barcodes of the droplets that successfully scored a TRUE.
validated_barcodes <- colnames(sce)[is_cell]


# --- UNDERSTANDING RUN TIME METRICS & RANGES ---
# ****************************************************************************#
# EXPECTED RESULTS (Based on raw 10x Genomics runs):
# - Total Droplets Tested: ~1,000,000 to 1,500,000 barcodes.
# - Cells Called: ~3,000 to 12,000 true cells (depending on target loading).
# - Yield Percentage: Typically 0.5% to 2.0% of all loaded raw droplets.
#
# INTERPRETING YOUR OUTPUT:
# - Droplets tested: 1,389,510 -> Normal partition capacity.
# - Cells called: 9,668 (0.7%) -> Great yield for a standard single channel.
# - Empty droplets: 1,379,842 (99.3%) -> Normal. Most captures contain soup.

cat("Droplets tested:", ncol(sce), "\n")
cat(
  "Cells called:", sum(is_cell),
  "(", round(sum(is_cell) / ncol(sce) * 100, 1), "%)\n"
)
cat(
  "Empty droplets removed:", sum(!is_cell),
  "(", round(sum(!is_cell) / ncol(sce) * 100, 1), "%)\n"
)


# --- VISUAL VALIDATION & DECODING LOG10(UMI + 1) ---
# ****************************************************************************#
# UNDERSTANDING LOG SCALE AXES:
# Single-cell UMI ranges are too vast for standard linear graphing. A plot with
# a cell at 50,000 UMIs squashes a cell with 100 UMIs flat against the y-axis.
# Log10 conversion spreads out the scale by counting the number of "zeros":
#   - Log value of 0 -> 1 UMI      (10^0)
#   - Log value of 1 -> 10 UMIs    (10^1)
#   - Log value of 2 -> 100 UMIs   (10^2) -> The EmptyDrops `lower` baseline
#   - Log value of 3 -> 1,000 UMIs (10^3) -> Small cells start appearing here
#   - Log value of 4 -> 10,000 UMIs(10^4) -> High-transcript, healthy cells
#
# HOW TO INTERPRET YOUR DATASET'S HISTOGRAM:
# - The Left Side (Typically Log 0 to ~2.0+): Expect massive, towering bars
#   classified as empty droplets (FALSE). These represent background matrix
#   noise and negligible ambient transcripts captured from the cell suspension.
# - The Right Side (Typically Log ~3.0 to 4.5+): Look for a distinct, broad
#   hill classified as true cells (TRUE). This shape represents your genuine,
#   biologically diverse population of low- to high-transcript cells.
# - The Mid-Plot Valley: A clean separation or "valley" between these two mass
#   distributions confirms a high-quality sample preparation, indicating that
#   the algorithm successfully distinguished background soup from healthy cells.

empty_df <- data.frame(
  total_umi = colSums(counts(sce)),
  is_cell = is_cell
)

p1 <- ggplot(empty_df, aes(x = log10(total_umi + 1), fill = is_cell)) +
  geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
  scale_fill_manual(
    values = c("TRUE" = "#2E86AB", "FALSE" = "#A23B72"),
    labels = c("Empty Droplet", "Cell"),
    name = "Classification"
  ) +
  labs(
    title = "EmptyDrops: Cell vs Empty Droplet Detection",
    subtitle = paste(sum(is_cell), "cells called from", ncol(sce), "droplets"),
    x = "log10(UMI + 1)",
    y = "Number of Droplets"
  ) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", size = 14))

ggsave(
  file.path(
    "3_output", RUN_ID, "qc-and-filtering",
    "plots", META_SAMPLE_NAME, "01_empty_droplets.png"
  ),
  plot = p1, width = 10, height = 6, dpi = 300
)


# --- PHYSICAL FILTERING ---
# ****************************************************************************#
# CRITICAL PIPELINE PRESERVATION STEP:
# We extract and save the `raw_counts_all_droplets` matrix BEFORE filtering.
# This unfiltered background matrix is strictly required by SoupX in Step 5,
# which needs to inspect the total "soup" background matrix to figure out
# exactly which ambient genes to subtract from your true cells.
raw_counts_all_droplets <- LayerData(seurat_obj, layer = "counts")

# Subset the Seurat object to securely retain only the validated cell barcodes
seurat_obj <- subset(seurat_obj, cells = validated_barcodes)

cat("After EmptyDrops:", ncol(seurat_obj), "cells retained\n")


# ****************************************************************************#
# SUMMARY & PIPELINE MILESTONE TRANSITION
# ****************************************************************************#
# WHERE WE STARTED:
# We loaded a raw matrix containing every single barcode generated by the 10x
# microfluidic channels. At that stage, we were looking at an unvetted data
# block containing over a million barcodes, where true cellular signals were
# entirely overwhelmed by an enormous volume of background technical noise.
#
# WHAT WE HAVE ACCOMPLISHED:
# In this step, we successfully executed our first major data purification
# milestone by running the statistical empty droplet filter (EmptyDrops). Out
# of the massive pool of raw barcodes, the algorithm modeled the expression
# pattern of the ambient background soup and ran a customized statistical test
# on every droplet. By setting a false discovery rate (FDR < 0.01) ceiling,
# we safely distinguished true cell profiles from fluid-only droplets. This
# allowed us to rescue rare, low-transcript quiescent cells (like resting
# lymphocytes) while confidently purging over 99% of empty background spaces.
#
# WHERE WE ARE HEADING (STEP 5):
# Right now, our dataset has been radically refined down from over a million
# empty spaces to a highly targeted, clean cellular matrix. However, while we
# have successfully removed completely empty droplets, the cells we saved are
# still technically contaminated by extracellular noise.
#
# Because thousands of fragile cells ruptured during initial tissue dissociation,
# their contents created a free-floating background molecular soup that was
# drawn into every single droplet—including our healthy ones. In Step 5, we will
# pass this filtered matrix to SoupX. We will use the un-filtered background
# count matrix we purposely cached right before filtering (`raw_counts_all_droplets`)
# to profile this background soup, calculate a sample-specific contamination
# fraction, and mathematically wash ambient transcripts directly out of our cells.
# ****************************************************************************#
