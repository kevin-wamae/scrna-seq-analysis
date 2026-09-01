# ****************************************************************************#
# STEP 4.2D: FastMNN Integration
# ****************************************************************************#


# --- HOW FASTMNN CORRECTS BATCH EFFECTS ---
# ****************************************************************************#
#   FastMNN finds mutual nearest neighbours — pairs of cells across two
#   batches that are each other's closest match — and uses them to compute
#   a batch-correction vector applied in low-dimensional space. Like CCA and
#   RPCA it operates on the merged object's split layers via
#   IntegrateLayers(), but unlike either it makes no assumption that
#   batches share similar overall structure, only that some cells are
#   genuinely alike across batches — making it the most conservative
#   corrector of the four: it leaves subtle biological differences intact
#   rather than risk smoothing them away.
#
# METHOD PROFILE (FastMNN):
#   - Strengths:      Preserves biological variation well; handles complex
#                     batch structures; less aggressive correction
#                     (conservative)
#   - Best for:       Datasets with subtle biological differences,
#                     hierarchical batch structures
#   - Computational
#     cost:           Medium
#   - When to use:    When biological signal is subtle, and you want
#                     conservative correction
#
# NOTE ON orig.reduction:
#   Unlike CCA/RPCA, FastMNN doesn't take an `orig.reduction` argument — it
#   corrects directly from the split layers rather than an existing PCA
#   embedding, so there's no uncorrected reduction to point it at.


# --- 1. Run FastMNN Integration ---
# ****************************************************************************#
integrated_fastmnn <- LOG_STEP("Running FastMNN integration...", {
    IntegrateLayers(
        object        = merged_seurat,
        method        = FastMNNIntegration,
        new.reduction = "integrated.mnn",   # corrected output embedding
        verbose       = FALSE
    )
})


# --- 2. Standard Downstream Workflow (On the Corrected Embedding) ---
# ****************************************************************************#
#   The same neighbors -> clusters -> UMAP chain as before, reading from
#   "integrated.mnn" so clusters and UMAP coordinates reflect
#   batch-corrected space. The UMAP is stored under its own name
#   ("umap.mnn") so each method's embedding survives side-by-side for the
#   integration comparison plots later in the pipeline — this is now the
#   fourth and final method added to that comparison set.
integrated_fastmnn <- LOG_STEP("Clustering and UMAP on FastMNN embedding...", {
    integrated_fastmnn %>%
        FindNeighbors(reduction = "integrated.mnn", dims = 1:30, verbose = FALSE) %>%
        FindClusters(resolution = 0.6, verbose = FALSE) %>%
        RunUMAP(reduction = "integrated.mnn", dims = 1:30,
                reduction.name = "umap.mnn", verbose = FALSE)
})

cat("FastMNN integration complete:", length(unique(integrated_fastmnn$seurat_clusters)), "clusters\n")
# FastMNN integration complete: 20 clusters


# --- 3. Checkpoint: Persist the FastMNN-Integrated Object to Disk ---
# ****************************************************************************#
#   `integrated_fastmnn` is needed again later (comparison UMAPs and mixing
#   metrics) — checkpointing here means a crashed job never repeats it.
#   This is the fourth and final per-method checkpoint in the pipeline.
#
#   The serialization format is governed by CHECKPOINT_FORMAT, set once in
#   Step 3.1 (1 = qs2, 2 = rds). qs2 multithreads its compression via
#   RcppParallel, reusing the N_WORKERS count also resolved in Step 3.1;
#   rds ignores thread counts.

# save the integrated object to disk in the chosen format
if (CHECKPOINT_FORMAT == 1) {
  CHECKPOINT_FILE <- file.path(
    DATA_CHECKPOINT_DIR, "05_integrated_fastmnn.qs2"
  )
  LOG_STEP(sprintf(
    "Saving FastMNN checkpoint (qs2, %d threads)...", N_WORKERS
  ), {
    qs2::qs_save(integrated_fastmnn, CHECKPOINT_FILE, nthreads = N_WORKERS)
  })

} else if (CHECKPOINT_FORMAT == 2) {
  CHECKPOINT_FILE <- file.path(
    DATA_CHECKPOINT_DIR, "05_integrated_fastmnn.rds"
  )
  LOG_STEP("Saving FastMNN checkpoint (rds, single-threaded)...", {
    saveRDS(integrated_fastmnn, CHECKPOINT_FILE)
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
#   Step 4.2C left us with three corrected embeddings (CCA, RPCA, Harmony)
#   spanning the anchor-based and iterative-adjustment approaches to batch
#   correction.
#
# WHAT WE HAVE ACCOMPLISHED:
#   We ran the fourth and final integration method in this pipeline.
#   FastMNN found mutual nearest neighbours across samples and wrote its
#   own embedding ("integrated.mnn") plus clustering and UMAP ("umap.mnn")
#   into the object, yielding 20 clusters — the most conservative result of
#   the four methods, consistent with FastMNN's bias toward preserving
#   subtle biological signal over aggressive mixing. The result is
#   checkpointed to disk, completing the set of four per-method checkpoints
#   alongside the naive baseline.
#
#   A NOTE ON WHAT WE DIDN'T RUN — scVI:
#   A fifth method, scVI, uses variational autoencoders to learn batch
#   corrections via deep learning, and can outperform linear methods on
#   very large or complex datasets. We deliberately skip it here: it
#   requires a Python/PyTorch/CUDA environment via reticulate, which is
#   slow to set up, prone to version conflicts, and often GPU-dependent —
#   overhead disproportionate to what this dataset needs. The four R-based
#   methods above are pure CRAN/Bioconductor installs, run without training
#   overhead, and are typically comparable in quality to scVI for datasets
#   of this scale. Readers who want to explore it anyway can start at the
#   scvi-tools docs (https://docs.scvi-tools.org/) or Seurat's bridge
#   integration article (https://satijalab.org/seurat/articles/seurat5_integration_bridge).
#
# WHERE WE ARE HEADING (STEP 4.3 — COMPARING INTEGRATION METHODS):
#   Four corrected embeddings now exist side-by-side with the naive,
#   uncorrected baseline — five objects in total, all checkpointed to disk:
#     1. Naive Merge (no integration — baseline)
#     2. CCA Integration
#     3. RPCA Integration
#     4. Harmony Integration
#     5. FastMNN Integration
#   None of these has yet been judged better or worse than the others.
#   Next, we compare all five side-by-side — visually, via UMAPs colored by
#   sample and condition, and quantitatively, via mixing scores — to
#   determine which correction (if any) best balances removing technical
#   batch effects against preserving real biological variation in this
#   dataset.
# ****************************************************************************#