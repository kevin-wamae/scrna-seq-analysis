# ****************************************************************************#
# STEP 4.2A: CCA Integration
# ****************************************************************************#


# --- HOW CCA CORRECTS BATCH EFFECTS ---
# ****************************************************************************#
#   Canonical Correlation Analysis identifies shared correlation structures
#   across datasets — projections where the samples look most alike — and
#   uses mutually-correlated "anchor" cells across batches to pull
#   equivalent cell states into alignment. It's Seurat's classic,
#   well-tested default: robust even when samples differ in cell type
#   composition, at the cost of being the slowest and most memory-hungry
#   of the four methods here.
#
# METHOD PROFILE (CCA):
#   - Strengths:      Excellent for diverse cell types; handles datasets
#                     with different cell type compositions; robust to
#                     varying cell numbers per sample
#   - Best for:       Experiments where samples may have different cell
#                     type proportions
#   - Computational
#     cost:           Medium
#   - When to use:    Default choice for most experiments
#
# WHAT IntegrateLayers() PRODUCES:
#   The counts/data layers are left untouched — integration writes a NEW
#   corrected low-dimensional embedding ("integrated.cca") alongside the
#   existing uncorrected "pca" reduction (`orig.reduction`). All downstream
#   steps must therefore point at this corrected reduction explicitly;
#   anything still reading "pca" would silently use uncorrected values.

# Integrate using CCA
integrated_cca <- LOG_STEP("Running CCA integration...", {
  IntegrateLayers(
    object         = merged_seurat,
    method         = CCAIntegration,
    orig.reduction = "pca",              # uncorrected input embedding
    new.reduction  = "integrated.cca",   # corrected output embedding
    dims           = 1:30,
    verbose        = FALSE
  )
})


# --- STANDARD DOWNSTREAM WORKFLOW (ON THE CORRECTED EMBEDDING) ---
# ****************************************************************************#
#   The same neighbors -> clusters -> UMAP chain as the naive merge (Step
#   3.3A), with one critical difference: every step reads from
#   "integrated.cca" rather than "pca", so clusters and UMAP coordinates
#   reflect batch-corrected space. The UMAP is stored under its own name
#   ("umap.cca") so each method's embedding survives side-by-side for the
#   Step 11 comparison plots.
integrated_cca <- LOG_STEP("Clustering and UMAP on CCA embedding...", {
  integrated_cca %>%
    FindNeighbors(reduction = "integrated.cca", dims = 1:30, verbose = FALSE) %>%
    FindClusters(resolution = 0.6, verbose = FALSE) %>%
    RunUMAP(reduction = "integrated.cca", dims = 1:30,
            reduction.name = "umap.cca", verbose = FALSE)
})

cat("CCA integration complete:", length(unique(integrated_cca$seurat_clusters)), "clusters\n")
# CCA integration complete: 23 clusters


# --- 3. Checkpoint: Persist the CCA-Integrated Object to Disk ---
# ****************************************************************************#
#   `integrated_cca` is needed again later (comparison UMAPs and mixing
#   metrics), and CCA is the most expensive integration in the pipeline —
#   checkpointing here means a crashed job never repeats it.
#
#   The serialization format is governed by CHECKPOINT_FORMAT, set once in
#   Step 3.1 (1 = qs2, 2 = rds). qs2 multithreads its compression via
#   RcppParallel, reusing the N_WORKERS count also resolved in Step 3.1;
#   rds ignores thread counts.

# save the integrated object to disk in the chosen format
if (CHECKPOINT_FORMAT == 1) {
    CHECKPOINT_FILE <- file.path(
        DATA_CHECKPOINT_DIR, "02_integrated_cca.qs2"
    )
    LOG_STEP(sprintf(
        "Saving CCA checkpoint (qs2, %d threads)...", N_WORKERS
    ), {
        qs2::qs_save(integrated_cca, CHECKPOINT_FILE, nthreads = N_WORKERS)
    })

} else if (CHECKPOINT_FORMAT == 2) {
    CHECKPOINT_FILE <- file.path(
        DATA_CHECKPOINT_DIR, "02_integrated_cca.rds"
    )
    LOG_STEP("Saving CCA checkpoint (rds, single-threaded)...", {
        saveRDS(integrated_cca, CHECKPOINT_FILE)
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
#   The preceding step delivered `merged_seurat`: all samples merged with
#   split layers, normalized, variable features selected, scaled, and an
#   uncorrected PCA computed — the common launchpad every integration
#   method departs from.
#
# WHAT WE HAVE ACCOMPLISHED:
#   We ran the first of four integration methods. CCA identified anchor
#   cells across samples and wrote a batch-corrected embedding
#   ("integrated.cca") alongside the uncorrected PCA, then we re-clustered
#   and re-embedded (umap.cca) in that corrected space, yielding 23 clusters.
#   The result is checkpointed to disk, so this — the most computationally
#   expensive method in the pipeline — never needs to be repeated.
#
# WHERE WE ARE HEADING (STEP 4.2B):
#   One corrected embedding is not yet evidence of good integration — each
#   method makes different correction/preservation trade-offs. Next we run
#   RPCA, a faster and more conservative anchor-based alternative, adding
#   its own embedding to the same object family so all methods can be
#   compared head-to-head against the naive baseline.
# ****************************************************************************#