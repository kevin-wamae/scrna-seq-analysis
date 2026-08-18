# ****************************************************************************#
# STEP 10: Save processed data and create summary
# ****************************************************************************#

# ****************************************************************************#
# --- LAYMAN CONTEXT & PURPOSE ---
# ****************************************************************************#
#   This final step is our pipeline's "shipping dock." Over the last nine steps,
#   we have taken a raw sequencing file cluttered with empty droplets, backgro-
#   und molecular soup, cell doublets, and dying cellular debris, and transfor-
#   med it into a pristine biological dataset.
#
#   Here, we wrap this clean data up into a compressed R master binary file
#   (.rds) so you never have to repeat the time-consuming clean-up steps again.
#   We also compile an automated "macro ledger" (QC_summary.csv) that acts as
#   a medical chart for your sample, tracking exactly how many cells and genes
#   survived each technical gate. This ensures absolute transparency and
#   mathematical reproducibility for any downstream analysis or manuscript
#   publication.


# --- 1. Save Clean Seurat Object (.rds Binary) ---
#   Preserves the entire processed object (raw counts, normalized layers,
#   and metadata) in a compressed R binary format. This serves as your
#   master file for downstream PCA, UMAP, and clustering analysis.
RDS_FILENAME <- paste0(META_SAMPLE_NAME, "_qc_filtered.rds")

# Print the initialization message
message("\nSaving cleaned Seurat object... ", appendLF = FALSE)
# Save the cleaned Seurat object
saveRDS(SEURAT_OBJ, file = file.path(DATA_OUT_DIR, RDS_FILENAME))
# Print the success confirmation message
message("✓ Done!\n")


# --- 2. Compile Comprehensive QC Summary (dplyr Tibble) ---
#   Consolidates sample-wide processing checkpoints into a single macro ledger.
#   Tracks step-by-step cell and gene survival counts through the entire QC
#   pipeline, capturing final median metrics and active filtering parameters.
#
# EXPECTED LEDGER BENCHMARKS (CROSS-REFERENCE YOUR QC_summary.csv COLUMNS):
#   - [final_cells]: 3,000-12,000 (Typical yielding scale for 10x channels)
#   - [median_umi]: 2,000-10,000 transcripts per cell
#   - [median_genes]: 1,000-3,000 unique features detected per cell
#   - [median_mt_pct]: <3-5% (Confirms structural membrane viability)
#   - [doublet_rate_pct]: Expected rate based on cell load scale
#   - [contamination_pct]: Remaining background ambient soup fraction
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
QC_SUMMARY <- tibble(
    sample = META_SAMPLE_NAME,

    # --- Cell Survival Tracking Across Processing Gates ---
    initial_droplets = INITIAL_DROPLET_COUNT,
    after_emptydrops = length(VALIDATED_BARCODES),
    after_soupx = length(VALIDATED_BARCODES),
    after_doublets = CELLS_BEFORE_DOUBLET_REMOVAL - N_DOUBLETS,
    after_cell_qc = ncol(SEURAT_OBJ),
    final_cells = ncol(SEURAT_OBJ),

    # --- Step-by-Step Cell Attrition Percentages ---
    # Initial: Massive evacuation of un-encapsulated empty background soup
    cell_removal_initial_pct = round(((initial_droplets - after_emptydrops) /
        initial_droplets) * 100, 2),
    # Secondary: Target pruning of physical doublets and dead ghost cells
    cell_removal_secondary_pct = round(((after_emptydrops - final_cells) /
        after_emptydrops) * 100, 2),

    # --- Feature (Gene) Resolution Dynamics ---
    initial_genes = INITIAL_GENE_COUNT,
    # Calculate intermediate count before explicit biological gene purging
    after_frequency_filter = sum(GENE_QC$pct_cells_detected >=
        MIN_PCT_CELLS),
    final_genes = nrow(SEURAT_OBJ),

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
    median_umi = median(SEURAT_OBJ$nCount_RNA),
    median_genes = median(SEURAT_OBJ$nFeature_RNA),
    median_mt_pct = median(SEURAT_OBJ$percent.mt),

    # --- Pipeline Math Parameters ---
    contamination_pct = ifelse(is.null(CONTAMINATION_FRACTION), 0,
        round(CONTAMINATION_FRACTION * 100, 2)
    ),
    doublet_rate_pct = round(DOUBLET_RATE, 2),

    # --- Reproducibility Constraint Logging ---
    ncount_min = NCOUNT_MIN,
    ncount_max = NCOUNT_MAX,
    nfeature_min = NFEATURE_MIN,
    nfeature_max = NFEATURE_MAX,
    mt_threshold = MT_THRESH
)


# --- 3. Export Macro Quality Audit Trail ---
#   Saves the consolidated summary metrics table as a flat CSV file. This
#   provides a permanent, human-readable record of sample quality and data
#   attrition required for publication transparency and methods sections.

SUMMARY_CSV_PATH <- file.path(
    METRICS_OUT_DIR, paste0(META_SAMPLE_NAME, "_QC_summary.csv")
)

write.csv(QC_SUMMARY, SUMMARY_CSV_PATH, row.names = FALSE)


# --- 4. Export Cell-Level Metadata Matrix ---
#   Exports a flat, barcode-by-barcode spreadsheet containing metrics for
#   every individual surviving cell (UMIs, gene counts, MT%, Ribo%). This
#   allows for external auditing or custom plot generation in other tools.
METADATA_CSV_PATH <- file.path(DATA_OUT_DIR, "cell_metadata.csv")
write.csv(SEURAT_OBJ@meta.data, METADATA_CSV_PATH)


# --- 5. Log Pipeline Finalization Status ---
# Prints a definitive terminal report summarizing the output file matrix,
# final cell/gene dimensions, and saved visualization diagnostic profiles.
cat("\n=== QC Pipeline Complete ===\n")
cat("Final cells:", ncol(SEURAT_OBJ), "\n")
cat("Final genes:", nrow(SEURAT_OBJ), "\n")
cat("Plots saved:", length(list.files(PLOTS_OUT_DIR)), "files\n\n")

cat("Files created:\n")
cat("  •", file.path(DATA_OUT_DIR, RDS_FILENAME), "- Clean Seurat object\n")
cat("  •", SUMMARY_CSV_PATH, "- Comprehensive QC summary table\n")
cat("  •", METADATA_CSV_PATH, "- Cell-level metadata matrix\n")


# ****************************************************************************#
# SUMMARY & PIPELINE MILESTONE TRANSITION
# ****************************************************************************#
# WHERE WE STARTED:
#   We began this analytical journey with a raw, unvetted sequencing matrix
#   generated by the 10x microfluidic chip—a massive block containing over a
#   million barcodes where valid biological signals were hidden behind a wall
#   of technical noise and experimental debris.
#
# WHAT WE HAVE ACCOMPLISHED:
#   With this final step, we have successfully completed the macro-phase of
#   "Data Purification, Filtering, and Quality Control." By progressing through
#   a sequence of high-performance filters, we evacuated empty background soup
#   (Step 4), assessed ambient RNA drift (Step 5), removed physical multi-cell
#   doublets (Step 6), eliminated dead or ruptured ghost cells (Step 7), and
#   purged uninformative gene rows and Hemoglobin transcripts (Step 8).
#   By normalizing cellular scaling depth and isolating the top 2,000 highly
#   variable biological features (Step 9), we have safely archived our master
#   processed Seurat binary file and generated human-readable audit ledgers.
#
# WHERE WE ARE HEADING (THE PROJECT SHIFT):
#   The quality control refinery is now officially locked down. Your single-
#   cell dataset is completely stable, verified, and structured for biological
#   discovery.
#
#   You are now exiting the data-cleaning stage and transitioning directly into
#   downstream "Unsupervised Mathematical Clustering and Dimensional Reduction."
#   In the next phase of your pipeline, you will pass this clean Seurat object
#   into Linear Scaling and Principal Component Analysis (PCA) to compress your
#   2,000 variable genes down into a compact space of 30 principal components.
#   This compressed matrix will then feed graph-based shared nearest neighbor
#   (SNN) community clustering and native 2D UMAP modelling, enabling you to
#   visually identify and annotate distinct biological cell lineages across your
#   cohort.
# ****************************************************************************#
