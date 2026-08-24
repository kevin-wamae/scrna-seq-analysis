# ****************************************************************************#
# STEP 3.3A: Naive merge without integration
# ****************************************************************************#


# --- THE PURPOSE OF A NAIVE BASELINE ---
# ****************************************************************************#
#   Before running Harmony or FastMNN, we first merge every sample together
#   with ZERO batch correction and run the standard clustering workflow on
#   top. This isn't wasted computation — it's a deliberate control.
#
# WHY WE BOTHER WITH AN UNCORRECTED BASELINE:
#   Integration methods can look impressive on a UMAP even when they aren't
#   doing much biologically meaningful work. Without first seeing what the
#   data looks like BEFORE correction, you have no reference point to judge
#   whether Harmony/FastMNN actually removed technical batch effects, or
#   whether there was never much batch effect to correct in the first place.
#   Every integration comparison plot in later steps will be measured against
#   this naive merge.


# --- 1. Merge All Samples ---
# ****************************************************************************#
#   `add.cell.ids` prefixes every cell barcode with its sample_id before
#   merging. This is not cosmetic — 10x barcodes (e.g. AAACCTGAGAAACCAT-1)
#   are only guaranteed unique WITHIN a single sequencing run. Without this
#   prefix, cells from different patients sharing the same raw barcode would
#   silently collide and overwrite each other during the merge.
merged_naive <- LOG_STEP("Merging all samples (naive, uncorrected)...", {
  merge(
    x = seurat_list[[1]],
    y = seurat_list[-1],
    add.cell.ids = names(seurat_list),
    project = "GSE174609_Naive_Merge"
  )
})

cat("Merged dataset:", ncol(merged_naive), "cells ×", nrow(merged_naive), "genes\n")
# Merged dataset: 72649 cells × 18861 genes


# --- 2. Run the Standard Seurat Workflow ---
# ****************************************************************************#
#   We push the merged-but-uncorrected object straight through the familiar
#   normalize -> variable features -> scale -> PCA -> UMAP -> neighbors ->
#   clusters chain, exactly as Part 2 did per-sample for SoupX calibration —
#   except now across the FULL, unintegrated cohort.
#
# PIPELINE STEP BREAKDOWN:
#   1. NormalizeData: Scales counts per cell to 10k and log-transforms them
#   2. FindVariableFeatures: Selects top 2,000 variance-driving genes
#   3. ScaleData: Shifts gene expression means to 0 and scales variance to 1
#   4. RunPCA: Compresses 2,000 variable genes down to 50 principal components
#   5. RunUMAP: Projects the top 30 PCs into 2D for visualization
#   6. FindNeighbors: Builds an SNN graph connecting cells with similar PCs
#   7. FindClusters: Uses the Louvain algorithm to partition the SNN graph

merged_naive <- LOG_STEP("Running standard workflow (normalize, PCA, UMAP, clustering)...", {
  merged_naive %>%
    NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(nfeatures = 2000, verbose = FALSE) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(npcs = 50, verbose = FALSE) %>%
    RunUMAP(dims = 1:30, reduction = "pca", verbose = FALSE) %>%
    FindNeighbors(dims = 1:30, verbose = FALSE) %>%
    FindClusters(resolution = 0.6, verbose = FALSE)
})

cat("Clustering complete:", length(unique(merged_naive$seurat_clusters)), "clusters identified\n")
# Clustering complete: 19 clusters identified


# --- 3. Checkpoint: Persist the Naive Baseline to Disk ---
# ****************************************************************************#
#   `merged_naive` is needed again much later (Step 11 comparison UMAPs and
#   Step 12 mixing metrics), so we checkpoint it now — a crashed job can
#   resume from here without redoing the merge + clustering workflow.
#
#   The serialization format is governed by CHECKPOINT_FORMAT, set once in
#   Step 3.1 (1 = qs2, 2 = rds). qs2 multithreads its compression via
#   RcppParallel, so it reuses the N_WORKERS count also resolved in Step 3.1
#   (slurm allocation - 1, or a 4-core cap on shared interactive machines);
#   rds ignores thread counts.

# save the merged object to disk in the chosen format
if (CHECKPOINT_FORMAT == 1) {
  CHECKPOINT_FILE <- file.path(
    DATA_CHECKPOINT_DIR, "01_merged_naive_clustered.qs2"
  )
  LOG_STEP(sprintf(
    "Saving naive merge checkpoint (qs2, %d threads)...", N_WORKERS
  ), {
    qs2::qs_save(merged_naive, CHECKPOINT_FILE, nthreads = N_WORKERS)
  })

} else if (CHECKPOINT_FORMAT == 2) {
  CHECKPOINT_FILE <- file.path(
    DATA_CHECKPOINT_DIR, "01_merged_naive_clustered.rds"
  )
  LOG_STEP("Saving naive merge checkpoint (rds, single-threaded)...", {
    saveRDS(merged_naive, CHECKPOINT_FILE)
  })
} else {
  stop("CHECKPOINT_FORMAT must be 1 (qs2) or 2 (rds)")
}

# print the path to the checkpoint file
cat("✓ Checkpoint written:", CHECKPOINT_FILE, "\n\n")


# ****************************************************************************#
# SUMMARY & PIPELINE MILESTONE TRANSITION
# ****************************************************************************#
# WHERE WE STARTED:
#   Before entering this step, Step 2 delivered 8 independently QC-filtered
#   Seurat objects sharing a common gene space — structurally ready to
#   combine, but still living as separate objects with no shared embedding.
#
# WHAT WE HAVE ACCOMPLISHED:
#   We merged all 8 samples into a single object with collision-safe cell
#   barcodes, then ran the full standard clustering workflow directly on top
#   — with zero batch correction applied. This gives us a concrete, uncorrected
#   baseline: 72,649 cells resolving into 19 clusters purely on raw biological
#   and technical variation combined.
#
# WHERE WE ARE HEADING (STEP 3.3B):
#   Because no batch correction was applied, some (or many) of these 19
#   clusters may be splitting cells by patient/technical origin rather than
#   true cell type — the exact confound integration exists to fix. In Step 3.3B,
#   we will run Harmony and/or FastMNN on this same merged object and compare
#   the resulting UMAPs and clustering side-by-side against this naive
#   baseline, making the effect of batch correction directly visible rather
#   than assumed.
# ****************************************************************************#