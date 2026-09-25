# ****************************************************************************#
# STEP 6.6: Save integrated and clustered data
# ****************************************************************************#


# --- PREREQUISITES (scripts that must run before this one in this session) ---
# ****************************************************************************#
#   • part-3.0-install-packages.R — packages for library() calls.
#   • part-3.1-load-libraries-and-configuration.R — RUN_ID, CHECKPOINT_FORMAT,
#     N_WORKERS, LOG_STEP, and the output-directory variables
#     (DATA_OUT_DIR / METADATA_OUT_DIR).
#   • part-5.3B-select-best-integration-method.R — supplies `best_method` and
#     `mixing_results`, both reported in the integration summary below.
#   • part-6.3-load-clustering-checkpoints.R — supplies `integrated_final` with
#     `seurat_clusters` set, and `optimal_resolution`.
#   The dependency manager is always re-sourced here (base-R only, cheap) so
#   the latest DEP_TABLE is loaded on every run; ensure_dependencies() does
#   nothing if those objects are already in the environment, otherwise it
#   offers to source them (interactive) or stops with a clear message
#   (non-interactive/batch). When everything is already satisfied, an
#   interactive session is asked whether to re-source fresh or proceed as-is.
source("1_scripts/part-3_integration-and-clustering/part-3.0-dependencies.R")
ensure_dependencies(step = "part-6.6-save-clustered-and-integrated-data.R")

# NOTE: Requires Step 6.3 to have loaded the clustering checkpoints in this
#       session (supplying `integrated_final` / `optimal_resolution`) and Step
#       5.3B for `best_method` / `mixing_results`.


# --- WHY THIS STEP EXISTS ---
# ****************************************************************************#
#   Everything up to here lived in the R session: the objects already on disk
#   are checkpoints (intermediates), and the clustered object exists only in
#   memory after Steps 6.3-6.4. This step produces the actual DELIVERABLE of
#   Part 3 — a single portable object plus its metadata
#   ledger — so downstream cell-type annotation and differential expression can
#   start from a clean, self-contained file in a new session.
#
#   The three outputs:
#     1. `integrated_clustered_seurat.<ext>` — the final object with
#        `seurat_clusters` at `optimal_resolution`. Layers are joined first
#        (below) so the file carries one assay layer rather than the per-sample
#        layers left by integration; a joined object is far more portable and
#        is what `FindMarkers()` / annotation steps expect.
#     2. `integrated_cell_metadata.csv` — per-cell metadata keyed by barcode
#        (row.names = TRUE), the tabular mirror of the object for quick
#        inspection without loading Seurat. This is a deliberate exception to
#        the `row.names = FALSE` convention: the barcodes ARE the row
#        identifier here, and losing them would make the table ambiguous.
#     3. `integration_summary.csv` — one row describing the run: dataset,
#        sample/cell/gene counts, chosen integration method and resolution,
#        cluster count, and the winning method's mixing/separation scores.
#
#   INTERPRETATION of the headline numbers:
#     - `mixing_score` / `condition_separation` are those selected in Step
#       5.3B; they describe the integration, not the clustering.
#     - `n_clusters` should land in the interpretable 8-15 range discussed in
#       Step 6.2B; a value far outside it is a prompt to revisit the resolution
#       choice, not a failure of this save step.
#
#   WHAT THE SAVED SUMMARY MEANS:
#   These fields document the analysis decision and preserve the evidence needed
#   to revisit it. They do not certify that every cluster is a cell type or that
#   treatment biology has been fully validated. Final interpretation still
#   requires the resolution grid, sample/condition composition, marker-gene
#   validation, and downstream within-cell-type comparisons. In particular,
#   reference thresholds such as a preferred cluster-count range or a mixing-score
#   reference are study-context prompts, not universal pass/fail rules.
#
#   CHECKPOINT DECISION: no additional checkpoint. The object written here is
#   itself the preserved artifact, and Steps 6.2B/6.3 already carry the
#   `06_clustered_final` / `07_clustering_metrics` checkpoints for resuming the
#   clustering chain. Its extension still honours `CHECKPOINT_FORMAT` so a run
#   never mixes on-disk serialization formats.

cat("\n=== Saving Final Results ===\n")


# --- 1. Join Split Layers for a Portable Object ---
# ****************************************************************************#
#   Integration left per-sample layers in the assay; `JoinLayers()` collapses
#   them into a single layer. This keeps the saved file compact and avoids the
#   "multiple layers" surprises that can bite downstream functions if the
#   object is later fed to tools expecting one counts/data layer.
integrated_final <- JoinLayers(integrated_final)


# --- 2. Save the Final Integrated & Clustered Object ---
# ****************************************************************************#
#   `CHECKPOINT_FORMAT` (set once in Step 3.1: 1 = qs2, 2 = rds) decides the
#   extension and serializer, exactly as it does for every checkpoint, so the
#   deliverable matches the run's chosen format.
if (CHECKPOINT_FORMAT == 1) {
    FINAL_OBJECT <- file.path(DATA_OUT_DIR, "integrated_clustered_seurat.qs2")
    LOG_STEP("Saving final integrated object (qs2)...", {
        qs2::qs_save(integrated_final, FINAL_OBJECT, nthreads = N_WORKERS)
    })
} else if (CHECKPOINT_FORMAT == 2) {
    FINAL_OBJECT <- file.path(DATA_OUT_DIR, "integrated_clustered_seurat.rds")
    LOG_STEP("Saving final integrated object (rds, single-threaded)...", {
        saveRDS(integrated_final, FINAL_OBJECT)
    })
} else {
    stop("CHECKPOINT_FORMAT must be 1 (qs2) or 2 (rds)")
}
cat("→ Saved:", FINAL_OBJECT, "\n")


# --- 3. Save the Cell Metadata Ledger ---
# ****************************************************************************#
#   `row.names = TRUE` preserves the cell barcodes as row identifiers (see the
#   note in WHY THIS STEP EXISTS). Every per-cell covariate — sample, condition,
#   QC metrics, and the `clusters_res_*` / `seurat_clusters` columns — travels
#   with the table.
CELL_METADATA_OUT <- file.path(METADATA_OUT_DIR, "integrated_cell_metadata.csv")
write.csv(integrated_final@meta.data,
          CELL_METADATA_OUT,
          row.names = TRUE)
cat("→ Saved:", CELL_METADATA_OUT, "\n")


# --- 4. Build & Save the Integration Summary ---
# ****************************************************************************#
#   One row that captures the identity of the run: how many samples/cells/genes,
#   which integration method won (Step 5.3B), the resolution chosen (Step 6.2B),
#   the resulting cluster count, and the winning method's two headline scores.
#   This is the machine-readable counterpart to the SUMMARY footer below.
integration_summary <- data.frame(
    dataset = "GSE174609",
    n_samples = length(unique(integrated_final$sample_id)),
    n_cells_total = ncol(integrated_final),
    n_genes = nrow(integrated_final),
    integration_method = best_method,
    optimal_resolution = optimal_resolution,
    n_clusters = length(unique(integrated_final$seurat_clusters)),
    mixing_score = mixing_results$mixing_score[mixing_results$method == best_method],
    condition_separation = mixing_results$condition_separation[mixing_results$method == best_method]
)

INTEGRATION_SUMMARY_OUT <- file.path(METADATA_OUT_DIR, "integration_summary.csv")
write.csv(integration_summary, INTEGRATION_SUMMARY_OUT, row.names = FALSE)
cat("→ Saved:", INTEGRATION_SUMMARY_OUT, "\n")


# ****************************************************************************#
# SUMMARY & PIPELINE MILESTONE TRANSITION
# ****************************************************************************#
# WHERE WE STARTED:
#   The final integrated-and-clustered object existed only in memory (restored
#   by Step 6.3 and described visually by Step 6.5); nothing downstream could
#   pick it up in a fresh session without re-running the clustering chain.
#
# WHAT WE HAVE ACCOMPLISHED:
#   We joined the assay layers and wrote Part 3's deliverable: the joined
#   `integrated_clustered_seurat.qs2`/`.rds` (honouring CHECKPOINT_FORMAT) in
#   `integrated_data/`, the per-cell metadata with barcodes as row names, and a
#   one-row `integration_summary.csv` recording the dataset size, the winning
#   integration method (`best_method`), the chosen resolution, the cluster
#   count, and the method's mixing/separation scores. No extra checkpoint was
#   written — this object IS the preserved artifact.
#
# WHERE WE ARE HEADING (NEXT: DOWNSTREAM ANNOTATION):
#   Part 3 is complete. Downstream work (Part 4, scaffolded in pixi.toml) starts
#   from `integrated_clustered_seurat`, assigning cell types via marker genes
#   and running condition-level differential expression on the clusters
#   established here.
# ****************************************************************************#
