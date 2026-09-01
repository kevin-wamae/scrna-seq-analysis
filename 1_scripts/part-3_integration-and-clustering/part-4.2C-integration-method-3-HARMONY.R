# ****************************************************************************#
# STEP 4.2C: Harmony Integration
# ****************************************************************************#


# --- HOW HARMONY CORRECTS BATCH EFFECTS ---
# ****************************************************************************#
#   Harmony takes a fundamentally different route from CCA and RPCA: rather
#   than finding anchor cells between samples, it iteratively adjusts cell
#   coordinates directly in PCA space until samples mix within shared
#   clusters — no anchor-finding, no layer-based machinery, and by far the
#   fastest and lightest of the four methods here.
#
# METHOD PROFILE (Harmony):
#   - Strengths:      Very fast; preserves global structure well; simple
#                     single-step process; scales to very large datasets
#   - Best for:       Large datasets, soft batch correction, preserving
#                     broad structure
#   - Computational
#     cost:           Low
#   - When to use:    First-pass integration, large datasets, when speed
#                     matters
#
# WHY THIS SCRIPT LOOKS DIFFERENT FROM CCA/RPCA:
#   Because Harmony works directly on a PCA embedding rather than per-sample
#   layers, it needs its own preparation: layers are joined into one
#   (JoinLayers()) and a fresh PCA is computed on the pooled, uncorrected
#   data before Harmony ever runs. This is a SEPARATE object
#   (`merged_harmony`) from `merged_seurat` used by CCA/RPCA — those two
#   methods need the split layers intact, so we branch off a joined copy
#   here rather than mutating the original.


# --- 1. Prepare a Joined, PCA-Reduced Object ---
# ****************************************************************************#
merged_harmony <- LOG_STEP("Joining layers and computing PCA for Harmony...", {
    merged_seurat %>%
        JoinLayers() %>%
        RunPCA(npcs = 50, verbose = FALSE)
})


# --- 2. Run Harmony Integration ---
# ****************************************************************************#
#   `RunHarmony()` corrects the PCA embedding in place, writing a new
#   "harmony" reduction into the object — the batch variable to correct
#   for ("sample_id") is passed directly rather than via `orig.reduction`/
#   `new.reduction` arguments, since Harmony doesn't share IntegrateLayers()'s
#   interface.
integrated_harmony <- LOG_STEP("Running Harmony integration...", {
    RunHarmony(merged_harmony, "sample_id")
})


# --- 3. Standard Downstream Workflow (On the Corrected Embedding) ---
# ****************************************************************************#
#   The same neighbors -> clusters -> UMAP chain as before, reading from
#   "harmony" rather than "pca" so clusters and UMAP coordinates reflect
#   batch-corrected space. The UMAP is stored under its own name
#   ("umap.harmony") so each method's embedding survives side-by-side for
#   the integration comparison plots later in the pipeline.
integrated_harmony <- LOG_STEP("Clustering and UMAP on Harmony embedding...", {
    integrated_harmony %>%
        FindNeighbors(reduction = "harmony", dims = 1:30, verbose = FALSE) %>%
        FindClusters(resolution = 0.6, verbose = FALSE) %>%
        RunUMAP(reduction = "harmony", dims = 1:30,
            reduction.name = "umap.harmony", verbose = FALSE)
})

cat("Harmony integration complete:", length(unique(integrated_harmony$seurat_clusters)), "clusters\n")
# Harmony integration complete: 21 clusters


# --- 4. Checkpoint: Persist the Harmony-Integrated Object to Disk ---
# ****************************************************************************#
#   `integrated_harmony` is needed again later (comparison UMAPs and mixing
#   metrics) — checkpointing here means a crashed job never repeats it.
#
#   The serialization format is governed by CHECKPOINT_FORMAT, set once in
#   Step 3.1 (1 = qs2, 2 = rds). qs2 multithreads its compression via
#   RcppParallel, reusing the N_WORKERS count also resolved in Step 3.1;
#   rds ignores thread counts.

# save the integrated object to disk in the chosen format
if (CHECKPOINT_FORMAT == 1) {
    CHECKPOINT_FILE <- file.path(
        DATA_CHECKPOINT_DIR, "04_integrated_harmony.qs2"
    )
    LOG_STEP(sprintf(
        "Saving Harmony checkpoint (qs2, %d threads)...", N_WORKERS
    ), {
        qs2::qs_save(integrated_harmony, CHECKPOINT_FILE, nthreads = N_WORKERS)
    })

} else if (CHECKPOINT_FORMAT == 2) {
    CHECKPOINT_FILE <- file.path(
        DATA_CHECKPOINT_DIR, "04_integrated_harmony.rds"
    )
    LOG_STEP("Saving Harmony checkpoint (rds, single-threaded)...", {
        saveRDS(integrated_harmony, CHECKPOINT_FILE)
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
#   Step 4.2B left us with two anchor-based corrected embeddings (CCA,
#   RPCA), both operating through IntegrateLayers() on the same split-layer
#   object.
#
# WHAT WE HAVE ACCOMPLISHED:
#   We ran the third of four integration methods, and the first
#   non-anchor-based one. Harmony joined the layers into a fresh object,
#   computed an uncorrected PCA, then iteratively adjusted that embedding
#   in place to remove batch structure — writing "harmony" plus its
#   clustering and UMAP ("umap.harmony") into the object, yielding 21
#   clusters (vs CCA's 23, RPCA's 22). The result is checkpointed to disk,
#   so this step never needs to be repeated. Three integration strategies —
#   two anchor-based, one iterative — are now represented side-by-side.
#
# WHERE WE ARE HEADING (STEP 4.2D — FASTMNN):
#   FastMNN returns to the anchor-based family, but with a different
#   philosophy from CCA/RPCA: rather than assuming samples share broad
#   structure, it finds mutual nearest neighbours — pairs of cells that are
#   each other's closest match across batches — and corrects conservatively
#   around them. Like CCA and RPCA, it operates through IntegrateLayers()
#   on the original `merged_seurat` object's split layers (not
#   `merged_harmony`), and is the fourth and final method added to our
#   comparison set.
# ****************************************************************************#