# ==============================================================================
# --- STEP 10: SAVE PROCESSED DATA AND CREATE SUMMARY ---
# ==============================================================================

# --- 1. Save Clean Seurat Object (.rds Binary) ---
# Preserves the entire processed object (raw counts, normalized layers,
# and metadata) in a compressed R binary format. This serves as your
# master file for downstream PCA, UMAP, and clustering analysis.
saveRDS(seurat_obj, file = "filtered_data/Healthy_1_qc_filtered.rds")

# --- 2. Compile Comprehensive QC Summary (dplyr Tibble) ---
# Consolidates sample-wide processing checkpoints into a single macro ledger.
# Tracks step-by-step cell and gene survival counts through the entire QC
# pipeline, capturing final median metrics and active filtering parameters.
#
# EXPECTED LEDGER BENCHMARKS (CROSS-REFERENCE YOUR QC_summary.csv COLUMNS):
#   - [final_cells]: 3,000-5,000 (Typical yielding scale for a PBMC lineage)
#   - [median_umi]: 2,000-10,000 transcripts per cell
#   - [median_genes]: 1,000-3,000 unique features detected per cell
#   - [median_mt_pct]: <3-5% (Confirms structural membrane viability)
#   - [doublet_rate_pct]: <2% remaining predicted physical multiplets
#   - [contamination_pct]: <5% remaining background ambient soup fraction
#
# EXPECTED ATTRITION FILTER RATES FOR CELL DATA AUDITING:
#   - [cell_removal_secondary_pct]: Traced through the following standard zones:
#       * 10.0% - 25.0%: Normal, healthy clean-up step (the golden sweet spot)
#       * <5.0%: Too lenient (danger of retaining low-complexity noise)
#       * 30.0% - 40.0%: Too strict (over-filtering; relax mt_threshold)
#       * >50.0%: Severe failure (catastrophic tissue cell lysis/debris)
#
# EXPECTED ATTRITION FILTER RATES FOR FEATURE (GENE) DATA AUDITING:
#   - [gene_removal_initial_pct]: Typically ranges from 30.0% to 50.0%. This
#     is normal; thousands of genomic features are completely unexpressed
#     or represent stochastic dropouts in quiet immune lineages.
#   - [gene_removal_secondary_pct]: Typically <1.0%. This represents the neat
#     extraction of specialized targeted contaminants like Hemoglobin.
qc_summary <- tibble(
    sample = sample_name,

    # --- Cell Survival Tracking Across Processing Gates ---
    initial_droplets = initial_droplet_count,
    after_emptydrops = length(validated_barcodes),
    after_soupx = length(validated_barcodes),
    after_doublets = cells_before_doublet_removal - n_doublets,
    after_cell_qc = ncol(seurat_obj),
    final_cells = ncol(seurat_obj),

    # --- Step-by-Step Cell Attrition Percentages ---
    # Initial: Massive evacuation of un-encapsulated empty background soup
    cell_removal_initial_pct = round(((initial_droplets - after_emptydrops) /
        initial_droplets) * 100, 2),
    # Secondary: Target pruning of physical doublets and dead ghost cells
    cell_removal_secondary_pct = round(((after_emptydrops - final_cells) /
        after_emptydrops) * 100, 2),

    # --- Feature (Gene) Resolution Dynamics ---
    initial_genes = initial_gene_count,
    # Calculate intermediate count before explicit biological gene purging
    after_frequency_filter = sum(gene_qc$pct_cells_detected >=
        min_pct_cells),
    final_genes = nrow(seurat_obj),

    # --- Step-by-Step Feature Attrition Percentages ---
    # Initial: Percentage of genome lost to low-frequency/dropout filters
    gene_removal_initial_pct = round(((initial_genes -
        after_frequency_filter) /
        initial_genes) * 100, 2),
    # Secondary: Percentage of filtered genes lost to explicit Hb exclusion
    gene_removal_secondary_pct = round(((after_frequency_filter -
        final_genes) /
        after_frequency_filter) * 100, 2),

    # --- Final Median Genomic Footprints ---
    median_umi = median(seurat_obj$nCount_RNA),
    median_genes = median(seurat_obj$nFeature_RNA),
    median_mt_pct = median(seurat_obj$percent.mt),

    # --- Pipeline Math Parameters ---
    contamination_pct = ifelse(is.null(contamination_fraction), 0,
        round(contamination_fraction * 100, 2)
    ),
    doublet_rate_pct = round(doublet_rate, 2),

    # --- Reproducibility Constraint Logging ---
    ncount_min = ncount_min,
    ncount_max = ncount_max,
    nfeature_min = nfeature_min,
    nfeature_max = nfeature_max,
    mt_threshold = mt_thresh
)

# --- 3. Export Macro Quality Audit Trail ---
# Saves the consolidated summary metrics table as a flat CSV file. This
# provides a permanent, human-readable record of sample quality and data
# attrition required for publication transparency and methods sections.
write.csv(qc_summary, "qc_metrics/QC_summary.csv", row.names = FALSE)

# --- 4. Export Cell-Level Metadata Matrix ---
# Exports a flat, barcode-by-barcode spreadsheet containing metrics for
# every individual surviving cell (UMIs, gene counts, MT%, Ribo%). This
# allows for external auditing or custom plot generation in other tools.
write.csv(seurat_obj@meta.data, "filtered_data/cell_metadata.csv")

# --- 5. Log Pipeline Finalization Status ---
# Prints a definitive terminal report summarizing the output file matrix,
# final cell/gene dimensions, and saved visualization diagnostic profiles.
cat("\n=== QC Pipeline Complete ===\n")
cat("Final cells:", ncol(seurat_obj), "\n")
cat("Final genes:", nrow(seurat_obj), "\n")
cat("Plots saved:", length(list.files("plots")), "files\n\n")

cat("Files created:\n")
cat("  • filtered_data/Healthy_1_qc_filtered.rds - Clean Seurat object\n")
cat("  • qc_metrics/QC_summary.csv - Comprehensive QC summary\n")
cat("  • filtered_data/cell_metadata.csv - Cell-level metadata\n")
