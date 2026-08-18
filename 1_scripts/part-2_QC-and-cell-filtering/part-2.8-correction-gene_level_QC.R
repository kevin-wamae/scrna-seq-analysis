# ****************************************************************************#
# STEP 8: Gene-level QC with detection threshold
# ****************************************************************************#

# Having successfully purified our cell population, we now pivot to cleaning our
# gene features. This phase shifts our focus from cell health to feature utility
#
# INTERPRETATION FRAMEWORK (GENERIC BALANCING BENCHMARKS):
# ==============================================================================
# 1. THE PROBLEM OF SPARSITY & "DROPOUTS"
#    - Single-cell RNA-seq captures only a small fraction of the true mRNA pool
#      present inside a cell. Consequently, many genes register a zero count
#      purely due to sampling luck rather than absent biology.
#    - However, if a gene registers as zero across almost your entire tissue
#      sample, it provides no mathematical variance. Keeping it introduces
#      pervasive background noise that misleads clustering algorithms
#
# 2. SELECTING A STRATEGY BASED ON YOUR DISCOVERY GOAL:
#    - Lenient Approach (0.1% of cells): Ideal for heterogeneous tissues or
#      tumor microenvironments. It ensures that rare transcripts unique to
#      highly specialized or low-abundance cell types are not accidentally wiped
#      from the dataset.
#    - Standard Approach (1.0% of cells): Ideal for uniform populations (like
#      pure cell lines or heavily characterized profiles). It aggressively
#      strips out background transcripts, accelerating downstream compute times.
#    - Conservative Approach (Fixed count, e.g., >=3 cells): A baseline filter
#      used to remove extreme singleton artifacts or mapping errors
#
# 3. HEMOGLOBIN EXCLUSION EXPLANATION (is_hb)
#    - Hemoglobin genes (^HB[AB]) originate from red blood cells. Mature RBCs
#      lack nuclei and should not be captured as intact single cells. High Hb
#      signals represent structural contamination or ambient cellular lysis
#      soup. We explicitly drop them to clean our downstream biological signal
# ==============================================================================

# TODO: include check for filtering from step 7 before continuing

# --- 1. Define the Detection Floor ---
#   Assign your minimum percentage filter (Adjust dynamically per project).
#   This establishes the baseline frequency required for a gene to be kept.
MIN_PCT_CELLS <- 0.1

# Calculate gene detection across remaining cell columns
COUNTS_MATRIX <- LayerData(SEURAT_OBJ, layer = "counts")
GENE_DETECTION <- rowSums(COUNTS_MATRIX > 0)

# Build diagnostic feature tracking spreadsheet
GENE_QC <- data.frame(
  gene               = rownames(SEURAT_OBJ),
  n_cells_detected   = GENE_DETECTION,
  pct_cells_detected = (GENE_DETECTION / ncol(SEURAT_OBJ)) * 100,
  is_mt              = grepl("^MT-", rownames(SEURAT_OBJ)),
  is_ribo            = grepl("^RP[SL]", rownames(SEURAT_OBJ)),
  is_hb              = grepl("^HB[AB]", rownames(SEURAT_OBJ))
)

# Print summary stats
cat("Total genes in matrix:", nrow(GENE_QC), "\n")
cat(
  "MT genes:", sum(GENE_QC$is_mt),
  "| Ribo genes:", sum(GENE_QC$is_ribo),
  "| Hb genes:", sum(GENE_QC$is_hb), "\n\n"
)

# Visualise of global gene detection distributions
P9 <- ggplot(
  # Filter out genes with zero detection before plotting
  # log10(0) = -Inf which introduces infinite values and drops rows
  GENE_QC %>% filter(pct_cells_detected > 0),
  aes(x = pct_cells_detected)
) +
  geom_histogram(bins = 50, fill = "#118AB2", alpha = 0.7) +
  # Use MIN_PCT_CELLS variable instead of hardcoded value so the line
  # always reflects the threshold actually set in the script
  geom_vline(xintercept = MIN_PCT_CELLS, linetype = "dashed", color = "red") +
  scale_x_log10(
    breaks = c(0.01, 0.1, 1, 10, 100),
    labels = scales::label_number(suffix = "%")
  ) +
  labs(
    title = "Gene Detection Across Cells",
    subtitle = paste0("Red line: ", MIN_PCT_CELLS, "% detection threshold"),
    x = "% of Cells Expressing Gene (log scale)",
    y = "Number of Genes"
  ) +
  theme_classic()

# --- Save Visual Diagnostics to Disk ---
ggsave(
  file.path(PLOTS_OUT_DIR, "06_gene_detection.png"),
  plot = P9, width = 8, height = 5, dpi = 300
)

# --- Console Status Update ---
cat("→ EXAMINE file:", file.path(PLOTS_OUT_DIR, "06_gene_detection.png"), "\n")
cat("   Choose threshold based on where detection drops off\n\n")


# ==============================================================================
# --- GENERIC VISUAL PLOT INTERPRETATION GUIDE ---
# What you are seeing in this histogram profile (06_gene_detection.png):
#   - Axis Scales: The y-axis tracks the total "Number of Genes". The x-axis
#     uses a log10 scale showing the "% of Cells Expressing a Gene", ranging
#     from the extreme rare fraction (left, e.g., 0.01%) up to ubiquitous
#     housekeeping genes expressed by 100% of your cell pool (right, 1e+02).
#   - The Noise Spikes (Far Left Columns): The towering bars on the extreme left
#     represent thousands of rare genes or technical fragments captured in only
#     one or two individual droplets across the entire dataset.
#   - The Biological Core (Right-Hand Waves): The broader wave of bars on
#     the right represents standard structural and tissue-specific genes
#     that characterize stable, functional cellular identities.
#   - The Dashed Red Cutoff Line: This line marks your chosen filter barrier.
#     Everything to the LEFT of this line consists of low-frequency transcripts
#     that will be purged from memory; everything to the RIGHT will be
#     preserved
# ==============================================================================

# ==============================================================================
# --- FEATURE FILTER EXECUTION & AUDITING PHASE ---
# ==============================================================================

# --- 2. Dynamic Metric Conversion ---
#   Single-cell dimensions fluctuate. This print statement calculates and checks
#   exactly how many real cells your percentage threshold translates to on the
#   fly (e.g., a 0.1% threshold across 8k cells means a gene must be in >= 8
#   cells).
cat("Setting gene filter threshold:\n")
cat("  Minimum detection: ≥", MIN_PCT_CELLS, "% of cells\n")
cat(
  "  (Equivalent to ≥", ceiling(ncol(SEURAT_OBJ) * MIN_PCT_CELLS / 100),
  "cells based on active dataset dimensions)\n\n"
)

# --- 3. Multi-Gate Boolean Assessment ---
# To pass into your kept pool, a gene must satisfy two distinct rules at once:
#   * Rule A: It must cross or equal your defined detection floor.
#   * Rule B: It must NOT be a Hemoglobin gene (!is_hb). RBCs lack nuclei, so
#     high Hb signals represent extracellular contamination or lysed background
#     soup. We explicitly drop them here to sanitize downstream biology.
GENES_TO_KEEP <- (GENE_QC$pct_cells_detected >= MIN_PCT_CELLS) & !GENE_QC$is_hb

# --- 4. Quality Control Reporting ---
#   Audits and logs your data loss metrics before making any permanent changes.
cat("  Genes passing filter:", sum(GENES_TO_KEEP), "/", nrow(GENE_QC), "\n")
cat("  Genes removed:\n")
cat("  - Low detection:", sum(GENE_QC$pct_cells_detected < MIN_PCT_CELLS), "\n")
cat("  - Hemoglobin:", sum(GENE_QC$is_hb), "\n\n")

# --- 5. Row-Wise Matrix Subsetting (Execution Phase) ---
#   Unlike cell filtering which deletes matrix columns, this operation slices
#   away matrix rows. Disqualified low-variance genes are permanently deleted
#   from memory, reducing RAM usage and preparing a clean, unpolluted gene
#   matrix for accurate, high-performance variable gene selection in the next
#   step
SEURAT_OBJ <- SEURAT_OBJ[GENE_QC$gene[GENES_TO_KEEP], ]

# --- 6. Final Feature Validation Check ---
# Confirms the final clean matrix dimensions available for downstream analysis.
cat("After gene filtering:", nrow(SEURAT_OBJ), "genes remaining in memory\n")


# ****************************************************************************#
# SUMMARY & PIPELINE MILESTONE TRANSITION
# ****************************************************************************#
# WHERE WE STARTED:
#   Before entering this step, we completed our cell-level quality control
#   (STEP 7) and successfully established a clean, high-viability population
#   of single cell barcodes. However, while our cell columns were pristine,
#   our gene rows were still completely unvetted—containing thousands of rare
#   technical fragments, zero-variance features, and extracellular contaminants
#   like loose Hemoglobin.
#
# WHAT WE HAVE ACCOMPLISHED:
#   In this step, we shifted our focus from cell health to feature utility by
#   executing a robust Gene-Level Quality Control filter. We audited the global
#   expression prevalence of every gene across our verified cell pool. Using a
#   log-scaled distribution histogram (06_gene_detection.png), we identified a
#   massive long-tail of low-frequency transcripts captured due to sampling luck
#   rather than functional biology. By applying a lenient 0.1% detection floor
#   and a hard biological block against red-blood-cell-derived Haemoglobin
#   transcripts (!is_hb), we successfully sliced away thousands of uninformat-
#   ive, zero-variance background rows without risking the loss of rare cell-
#   type markers.
#
# WHERE WE ARE HEADING (STEP 9):
#   Our dataset is now fully purified along both dimensions: we have optimized
#   cell columns and a sanitized, high-signal gene feature matrix. However, the
#   raw counts inside this clean matrix are still heavily biased by variable
#   technical sequencing depths across individual cells.
#
#   In Step 9, we will transition out of the quality control phase and enter
#   downstream processing by running Log-Normalization and Variable Feature
#   Selection. We will scale all cellular transcript volumes to a standard
#   factor of 10,000 UMIs to make them statistically comparable, apply a
#   log-transformation to stabilize variance, and deploy a variance-
#   stabilizing transformation (vst) to select the top 2,000 highly variable
#   genes that drive true biological clustering.
# ****************************************************************************#
