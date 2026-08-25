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