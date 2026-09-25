# ****************************************************************************#
# STEP 6.1: Perform clustering at multiple resolutions
# ****************************************************************************#


# --- PREREQUISITES (scripts that must run before this one in this session) ---
# ****************************************************************************#
#   • part-3.0-install-packages.R — packages for library() calls.
#   • part-3.1-load-libraries-and-configuration.R — RUN_ID, LOG_STEP, and the
#     output-directory variables.
#   • part-5.3B-select-best-integration-method.R — supplies `integrated_final`
#     (the winning method's Seurat object) and `reduction_final` (its UMAP
#     reduction name), committed to as the single integration result every
#     downstream step builds on.
#   The dependency manager is always re-sourced here (base-R only, cheap) so
#   the latest DEP_TABLE is loaded on every run; ensure_dependencies() does
#   nothing if `integrated_final` and `reduction_final` are already in the
#   environment, otherwise it offers to source the prerequisites (interactive)
#   or stops with a clear message (non-interactive/batch). When everything is
#   already satisfied, an interactive session is asked whether to re-source
#   fresh or proceed as-is.
source("1_scripts/part-3_integration-and-clustering/part-3.0-dependencies.R")
ensure_dependencies(step = "part-6.1-clustering-multiple-resolutions.R")

# NOTE: Requires Step 5.3B to have already run in this session — that step is
#       what selects `integrated_final` / `reduction_final` out of
#       `methods_list`.


# --- WHAT CLUSTERING ACTUALLY DOES (FOR BEGINNERS) ---
# ****************************************************************************#
#   Integration corrected for batch effects so cells now sit in a shared
#   space based on biology rather than which sample they came from.
#   Clustering is the next question: within that shared space, which cells
#   look enough alike to call them the same population? Seurat answers this
#   with graph-based community detection (the Louvain algorithm), in three
#   conceptual steps:
#
#     1. Build a k-nearest-neighbour (KNN) graph: each cell connects to its
#        ~20 most similar neighbours, measured in the INTEGRATED PCA space
#        (not raw gene expression, which would still carry batch effects).
#     2. Convert to a shared-nearest-neighbour (SNN) graph: connections are
#        reweighted by how many neighbours two cells have in common — cells
#        with lots of shared neighbours get a stronger edge between them,
#        which makes the resulting clusters more robust to noise.
#     3. Detect communities on that graph: the algorithm looks for groups of
#        cells that are much more densely connected to each other than to
#        the rest of the graph and calls each group a cluster.
#
#   `FindNeighbors()`, called back when each integration method was run
#   (Section 4), already built that SNN graph. `FindClusters()` below does
#   only step 3 — it partitions the existing graph, it does not rebuild it —
#   which is why sweeping five resolutions here is cheap.


# --- WHY WE SWEEP RESOLUTIONS INSTEAD OF PICKING ONE ---
# ****************************************************************************#
#   The `resolution` argument to `FindClusters()` controls how easily the
#   algorithm splits cells into separate communities: higher resolution =
#   more, smaller clusters; lower resolution = fewer, broader ones. Every
#   integration method so far used a single fixed resolution (0.6) purely so
#   the five methods could be compared on equal footing in Section 5 — that
#   value was never chosen for what this dataset's biology actually needs.
#
#   ROUGH RULES OF THUMB (more cells generally supports a finer resolution
#   without over-splitting):
#     - ~3,000 cells   → resolution 0.4-0.6
#     - ~10,000 cells  → resolution 0.6-0.8
#     - ~50,000 cells  → resolution 0.8-1.2
#     - >100,000 cells → resolution 1.0-1.5
#
#   Our dataset sits at ~72,000 cells — in between the 50k and 100k+ rows
#   above — so rather than guess a single value from the table, we sweep the
#   whole 0.4-1.2 range and keep every result. The right resolution to
#   actually use gets chosen later, once we can see what each one produces
#   (guide §7.5-7.6) — not decided blind here.
#
#   A COMMON MISCONCEPTION TO AVOID: it's tempting to tune the resolution
#   until the cluster count matches the number of cell types you expect
#   (e.g. "I know there are 8 PBMC cell types, so I want 8 clusters"). This
#   doesn't work, because cell types exist at multiple levels at once — "T
#   cells" might be 1 cluster at low resolution, but split into CD4+/CD8+ at
#   medium resolution, and further into naive/memory/effector subtypes at
#   high resolution. Cell states (activated vs resting, cell-cycle phase,
#   etc.) add further splits on top of that. None of these are "wrong" —
#   getting 12-15 clusters out of 8 expected cell types is normal and often
#   correct once subtypes are accounted for. There is no resolution that
#   uniquely reproduces "the" cell types; the right choice depends on how
#   fine-grained an answer your analysis actually needs.
#
#   GUIDE §9.2 / §10.3 — HOW TO READ THE RESOLUTION SERIES:
#   The low, medium, and high labels describe progressively finer views of the
#   same integrated neighbourhood graph, not three competing truths. Lower
#   resolutions are useful for broad populations; middle resolutions often
#   expose common subtypes; higher resolutions may reveal cell states or split
#   continuous transitions into smaller communities. The expected cluster-count
#   ranges in the guide are PBMC-oriented reference points, not acceptance
#   criteria for this cohort. What matters here is whether structure appears
#   gradually and remains biologically interpretable, rather than whether the
#   count matches a preconceived list of cell types.
#
#   A useful visual question for Step 6.2A is whether a new high-resolution
#   cluster is a coherent, marker-supported population or merely a thin slice
#   of a neighbouring cloud. The code below deliberately preserves every
#   resolution so that question can be answered later; this step does not
#   declare any resolution biologically correct.


# --- 1. Define the Resolution Sweep ---
# ****************************************************************************#
#   Five evenly-spaced values spanning the low/medium/high ranges described
#   above (guide §7.2, §7.4), matched to our ~72,000-cell dataset.
resolutions <- c(0.4, 0.6, 0.8, 1.0, 1.2)


# --- 2. Cluster at Each Resolution ---
# ****************************************************************************#
#   Each pass through the loop re-partitions the SAME SNN graph at a
#   different resolution — no neighbours are recomputed, so this is fast
#   even though it runs five times.
#
#   COLUMN NAMING: `FindClusters()` stores each run's cluster assignment
#   under Seurat's internal name "RNA_snn_res.<res>", and also overwrites the
#   object's "active" cluster identity (`Idents()`) with whichever run just
#   finished. Instead of leaving cluster labels under that internal name, we
#   immediately copy each run's result into its own clearly-named column
#   ("clusters_res_<res>") via `AddMetaData()`. That way all five results
#   survive side by side, and downstream steps can read
#   `integrated_final$clusters_res_0.8`, say, directly and unambiguously.
#
#   NOTE: because `Idents()` only ever holds the LAST resolution run (1.2),
#   downstream steps must read the explicit `clusters_res_*` columns below —
#   never `integrated_final$seurat_clusters` — to get a specific resolution.
for (res in resolutions) {
    integrated_final <- FindClusters(
        integrated_final,
        resolution = res,
        verbose = FALSE
    )

    cluster_idents <- Idents(integrated_final)
    integrated_final <- AddMetaData(
        integrated_final,
        metadata = as.character(cluster_idents),
        col.name = paste0("clusters_res_", res)
    )

    cat(sprintf("Resolution %.1f: %d clusters\n", res, nlevels(cluster_idents)))
}

cat("\nCluster assignments stored in integrated_final$clusters_res_* \n")
cat("(one column per resolution; seurat_clusters now reflects res = 1.2 only)\n")


# ****************************************************************************#
# SUMMARY & PIPELINE MILESTONE TRANSITION
# ****************************************************************************#
# WHERE WE STARTED:
#   Step 5.3B committed the pipeline to one integrated object
#   (`integrated_final`), but the object carried only a single clustering
#   from a fixed, arbitrary resolution (0.6) that downstream steps would
#   have silently inherited without ever being examined against our actual
#   ~72,000-cell dataset.
#
# WHAT WE HAVE ACCOMPLISHED:
#   Without touching cell type labels or biology yet, we re-partitioned
#   `integrated_final`'s existing SNN graph at five resolutions (0.4, 0.6,
#   0.8, 1.0, 1.2), printing the community count each one recovers and
#   storing every partition under its own explicit `clusters_res_<res>`
#   metadata column. The multi-resolution landscape now lives in one
#   object, ready to be examined — and remember, a "right" number of
#   clusters doesn't exist in isolation; it depends on which resolution's
#   output actually holds up biologically.
#
# WHERE WE ARE HEADING (NEXT: VISUALIZING MULTI-RESOLUTION CLUSTERING):
#   Five columns of cluster labels aren't yet something you can evaluate by
#   eye — a resolution that looks reasonable as a raw number (e.g. "18
#   clusters") could still be splitting one real cell type into
#   near-duplicate sub-clusters, or merging two distinct ones together. The
#   next step renders each resolution's clusters as its own UMAP panel so
#   the five options can be visually compared side by side (guide §7.5),
#   before any decision is made about which resolution to carry forward.
# ****************************************************************************#
