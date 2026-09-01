# ****************************************************************************#
# STEP 4.2B: RPCA Integration
# ****************************************************************************#


# --- HOW RPCA CORRECTS BATCH EFFECTS ---
# ****************************************************************************#
#   Reciprocal PCA projects each sample into every other sample's PCA space
#   (rather than deriving a shared CCA space from scratch) and finds anchor
#   cells in those reciprocal projections, identifying the variation the
#   samples share. This makes it markedly faster and lighter on memory than
#   CCA — the trade-off is a more conservative correction that assumes
#   samples share broadly similar cell type compositions.
#
# METHOD PROFILE (RPCA):
#   - Strengths:      Faster than CCA; less memory-intensive; scales well
#                     to large datasets (>100,000 cells)
#   - Best for:       Large-scale studies; similar cell type compositions
#                     across samples
#   - Computational
#     cost:           Low-Medium
#   - When to use:    When computational efficiency is critical, or
#                     datasets are very large
#
# NOTE ON PARALLEL PROCESSING & MEMORY:
#   RPCA's anchor-finding uses the `future` package internally, which is
#   why Step 3.1 configures `future.globals.maxSize` — that cap governs how
#   much data `future` may serialize for worker transfers (sized from the
#   SLURM allocation, or the conservative local default). We keep the plan
#   sequential here, so no worker copies of the multi-GB object are made;
#   if you nevertheless hit "total size of globals exceeds maximum" errors,
#   raise the cap (see Step 3.1) or reduce the dataset size.

# Integrate using RPCA
integrated_rpca <- LOG_STEP("Running RPCA integration...", {
    IntegrateLayers(
        object         = merged_seurat,
        method         = RPCAIntegration,
        orig.reduction = "pca",               # uncorrected input embedding
        new.reduction  = "integrated.rpca",   # corrected output embedding
        dims           = 1:30,
        verbose        = FALSE
    )
})


# --- 2. Standard Downstream Workflow (On the Corrected Embedding) ---
# ****************************************************************************#
#   The same neighbors -> clusters -> UMAP chain as before, reading from
#   "integrated.rpca" rather than "pca" so clusters and UMAP coordinates
#   reflect batch-corrected space. The UMAP is stored under its own name
#   ("umap.rpca") so each method's embedding survives side-by-side for the
#   integration comparison plots later in the pipeline.
integrated_rpca <- LOG_STEP("Clustering and UMAP on RPCA embedding...", {
    integrated_rpca %>%
        FindNeighbors(reduction = "integrated.rpca", dims = 1:30, verbose = FALSE) %>%
        FindClusters(resolution = 0.6, verbose = FALSE) %>%
        RunUMAP(reduction = "integrated.rpca", dims = 1:30,
            reduction.name = "umap.rpca", verbose = FALSE)
})

cat("RPCA integration complete:", length(unique(integrated_rpca$seurat_clusters)), "clusters\n")
# RPCA integration complete: 22 clusters


# --- 3. Checkpoint: Persist the RPCA-Integrated Object to Disk ---
# ****************************************************************************#
#   `integrated_rpca` is needed again later (comparison UMAPs and mixing
#   metrics) — checkpointing here means a crashed job never repeats it.
#
#   The serialization format is governed by CHECKPOINT_FORMAT, set once in
#   Step 3.1 (1 = qs2, 2 = rds). qs2 multithreads its compression via
#   RcppParallel, reusing the N_WORKERS count also resolved in Step 3.1;
#   rds ignores thread counts.

# save the integrated object to disk in the chosen format
if (CHECKPOINT_FORMAT == 1) {
    CHECKPOINT_FILE <- file.path(
        DATA_CHECKPOINT_DIR, "03_integrated_rpca.qs2"
    )
    LOG_STEP(sprintf(
        "Saving RPCA checkpoint (qs2, %d threads)...", N_WORKERS
    ), {
        qs2::qs_save(integrated_rpca, CHECKPOINT_FILE, nthreads = N_WORKERS)
    })

} else if (CHECKPOINT_FORMAT == 2) {
    CHECKPOINT_FILE <- file.path(
        DATA_CHECKPOINT_DIR, "03_integrated_rpca.rds"
    )
    LOG_STEP("Saving RPCA checkpoint (rds, single-threaded)...", {
        saveRDS(integrated_rpca, CHECKPOINT_FILE)
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
#   Step 4.2A left us with one corrected embedding (CCA) — a strong but
#   expensive anchor-based correction, and on its own no way to know
#   whether its trade-offs suit this dataset.
#
# WHAT WE HAVE ACCOMPLISHED:
#   We ran the second of four integration methods. RPCA found anchors in
#   reciprocal PCA projections — a faster, lighter, more conservative
#   correction than CCA — and wrote its own embedding ("integrated.rpca")
#   plus clustering and UMAP ("umap.rpca") into the object, yielding 22
#   clusters (vs CCA's 23). The result is checkpointed to disk, so this
#   step never needs to be repeated. Both anchor-based methods are now
#   represented side-by-side, ready for head-to-head comparison.
#
# WHERE WE ARE HEADING (STEP 4.2C — HARMONY):
#   Both methods so far are anchor-based and operate through Seurat's
#   IntegrateLayers() on split layers. Harmony takes a different route
#   entirely: it works directly on the PCA embedding, iteratively nudging
#   cell coordinates to remove batch structure — no anchors, no layer
#   machinery, and by far the fastest of the four. It therefore needs a
#   slightly different preparation (layers joined via JoinLayers() before
#   RunHarmony()), which the next script performs before adding a third
#   corrected embedding to our comparison set.
# ****************************************************************************#