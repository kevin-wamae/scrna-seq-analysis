#-----------------------------------------------
# STEP 8: Gene-level QC with detection threshold
#-----------------------------------------------
# Having successfully purified our cell population, we now pivot to cleaning our
# gene features. This phase shifts our focus from cell health to feature utility.
#
# INTERPRETATION FRAMEWORK (GENERIC BALANCING BENCHMARKS):
# ==============================================================================
# 1. THE PROBLEM OF SPARSITY & "DROPOUTS"
#    - Single-cell RNA-seq captures only a small fraction of the true mRNA pool
#      present inside a cell. Consequently, many genes register a zero count
#      purely due to sampling luck rather than absent biology.
#    - However, if a gene registers as zero across almost your entire tissue
#      sample, it provides no mathematical variance. Keeping it introduces
#      pervasive background noise that misleads clustering algorithms.
#
# 2. SELECTING A STRATEGY BASED ON YOUR DISCOVERY GOAL:
#    - Lenient Approach (0.1% of cells): Ideal for heterogeneous tissues or
#      tumor microenvironments. It ensures that rare transcripts unique to
#      highly specialized or low-abundance cell types are not accidentally wiped
#      from the dataset.
#    - Standard Approach (1.0% of cells): Ideal for uniform populations (like
#      pure cell lines or heavily characterized profiles). It aggressively strips
#      out background transcripts, accelerating downstream compute times.
#    - Conservative Approach (Fixed count, e.g., >=3 cells): A baseline filter
#      used to remove extreme singleton artifacts or mapping errors.
#
# 3. HEMOGLOBIN EXCLUSION EXPLANATION (is_hb)
#    - Hemoglobin genes (^HB[AB]) originate from red blood cells. Mature RBCs
#      lack nuclei and should not be captured as intact single cells. High Hb
#      signals represent structural contamination or ambient cellular lysis soup.
#      We explicitly drop them to clean our downstream biological signal.
# ==============================================================================

# Calculate gene detection across remaining cell columns
counts_matrix <- LayerData(seurat_obj, layer = "counts")
gene_detection <- rowSums(counts_matrix > 0)

# Build diagnostic feature tracking spreadsheet
gene_qc <- data.frame(
  gene               = rownames(seurat_obj),
  n_cells_detected   = gene_detection,
  pct_cells_detected = (gene_detection / ncol(seurat_obj)) * 100,
  is_mt              = grepl("^MT-", rownames(seurat_obj)),
  is_ribo            = grepl("^RP[SL]", rownames(seurat_obj)),
  is_hb              = grepl("^HB[AB]", rownames(seurat_obj))
)

cat("Total genes in matrix:", nrow(gene_qc), "\n")
cat(
  "MT genes:", sum(gene_qc$is_mt),
  "| Ribo genes:", sum(gene_qc$is_ribo),
  "| Hb genes:", sum(gene_qc$is_hb), "\n\n"
)

# Visualize global gene detection distributions
p9 <- ggplot(gene_qc, aes(x = pct_cells_detected)) +
  geom_histogram(bins = 50, fill = "#118AB2", alpha = 0.7) +
  geom_vline(xintercept = 0.1, linetype = "dashed", color = "red") +
  scale_x_log10() +
  labs(
    title = "Gene Detection Across Cells",
    subtitle = "Red line: Defined % of cells threshold",
    x = "% of Cells Expressing Gene (log scale)",
    y = "Number of Genes"
  ) +
  theme_classic()

ggsave("plots/06_gene_detection.png", p9, width = 8, height = 5, dpi = 300)

cat("→ EXAMINE plots/06_gene_detection.png\n")
cat("   Choose threshold based on where detection drops off\n\n")


# ==============================================================================
# --- generic visual plot interpretation guide ---
# What you are seeing in this histogram profile (06_gene_detection.png):
#   - Axis Scales: The y-axis tracks the total "Number of Genes". The x-axis uses
#     a log10 scale showing the "% of Cells Expressing a Gene", ranging from the
#     extreme rare fraction (left, e.g., 0.01%) up to ubiquitous housekeeping
#     genes expressed by 100% of your cell pool (right, 1e+02).
#   - The Noise Spikes (Far Left Columns): The towering bars on the extreme left
#     represent thousands of rare genes or technical fragments captured in only
#     one or two individual droplets across the entire dataset.
#   - The Biological Core (Right-Hand Waves): The broader wave of bars on the
#     right represents standard structural and tissue-specific genes that
#     characterize stable, functional cellular identities.
#   - The Dashed Red Cutoff Line: This line marks your chosen filter barrier.
#     Everything to the LEFT of this line consists of low-frequency transcripts
#     that will be purged from memory; everything to the RIGHT will be preserved.
# ==============================================================================

# ==============================================================================
# --- FEATURE FILTER EXECUTION & AUDITING PHASE ---
# ==============================================================================

# --- 1. Define the Detection Floor ---
# Assign your minimum percentage filter (Adjust dynamically per project).
# This establishes the baseline frequency required for a gene to be kept.
min_pct_cells <- 0.1

# --- 2. Dynamic Metric Conversion ---
# Single-cell dimensions fluctuate. This print statement calculates and checks
# exactly how many real cells your percentage threshold translates to on the fly
# (e.g., a 0.1% threshold across 8,000 cells means a gene must be in >= 8 cells).
cat("Setting gene filter threshold:\n")
cat("  Minimum detection: ≥", min_pct_cells, "% of cells\n")
cat(
  "  (Equivalent to ≥", ceiling(ncol(seurat_obj) * min_pct_cells / 100),
  "cells based on active dataset dimensions)\n\n"
)

# --- 3. Multi-Gate Boolean Assessment ---
# To pass into your kept pool, a gene must satisfy two distinct rules at once:
#   * Rule A: It must cross or equal your defined detection floor.
#   * Rule B: It must NOT be a Hemoglobin gene (!is_hb). RBCs lack nuclei, so
#     high Hb signals represent extracellular contamination or lysed background
#     soup. We explicitly drop them here to sanitize downstream biology.
genes_to_keep <- (gene_qc$pct_cells_detected >= min_pct_cells) & !gene_qc$is_hb

# --- 4. Quality Control Reporting ---
# Audits and logs your data loss metrics before making any permanent changes.
cat("Genes passing filter:", sum(genes_to_keep), "/", nrow(gene_qc), "\n")
cat("Genes removed:\n")
cat("  Low detection:", sum(gene_qc$pct_cells_detected < min_pct_cells), "\n")
cat("  Hemoglobin:", sum(gene_qc$is_hb), "\n\n")

# --- 5. Row-Wise Matrix Subsetting (Execution Phase) ---
# Unlike cell filtering which deletes matrix columns, this operation slices
# away matrix rows. Disqualified low-variance genes are permanently deleted
# from memory, reducing RAM usage and preparing a clean, unpolluted gene matrix
# for accurate, high-performance variable gene selection in the next step.
seurat_obj <- seurat_obj[gene_qc$gene[genes_to_keep], ]

# --- 6. Final Feature Validation Check ---
# Confirms the final clean matrix dimensions available for downstream analysis.
cat("After gene filtering:", nrow(seurat_obj), "genes remaining in memory\n")
