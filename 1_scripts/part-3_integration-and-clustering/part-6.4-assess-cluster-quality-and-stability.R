# ****************************************************************************#
# STEP 6.4: Assess cluster quality and stability
# ****************************************************************************#


# --- PREREQUISITES (scripts that must run before this one in this session) ---
# ****************************************************************************#
#   • part-3.0-install-packages.R — packages for library() calls (`cluster`,
#     loaded in Step 3.1, is used for the distance/silhouette concepts these
#     checks build on).
#   • part-3.1-load-libraries-and-configuration.R — RUN_ID, LOG_STEP, and the
#     output-directory variables.
#   • part-6.3-load-clustering-checkpoints.R — loads `integrated_final` with
#     `seurat_clusters` already committed to the chosen `optimal_resolution`
#     partition (from the Step 6.2B checkpoint), which is exactly the
#     clustering quality-checked below.
#   The dependency manager is always re-sourced here (base-R only, cheap) so
#   the latest DEP_TABLE is loaded on every run; ensure_dependencies() does
#   nothing if those objects are already in the environment, otherwise it
#   offers to source them (interactive) or stops with a clear message
#   (non-interactive/batch). When everything is already satisfied, an
#   interactive session is asked whether to re-source fresh or proceed as-is.
source("1_scripts/part-3_integration-and-clustering/part-3.0-dependencies.R")
ensure_dependencies(step = "part-6.4-assess-cluster-quality-and-stability.R")

# NOTE: Requires Steps 6.1 -> 6.2B to have already run once (so the clustering
#       checkpoints exist) and Step 6.3 to have loaded them in this session —
#       that is what commits `seurat_clusters` to the selected resolution.

# --- WHY THIS STEP EXISTS ---
# ****************************************************************************#
#   Step 6.1 swept five resolutions, Step 6.2A showed them on UMAP, and Step
#   6.2B picked one by silhouette score. A good silhouette score means the
#   clusters are well-separated in the integrated embedding — but it says
#   nothing about whether those clusters are biologically sensible or merely
#   artifacts of the clustering itself. This step performs two sanity checks
#   on the chosen partition (guide §7.7):
#
#     1. QUALITY (cluster sizes): are any clusters implausibly small? A tiny
#        cluster (<1% of cells) is either a genuine rare population or an
#        over-split piece of a larger group — the flag exists so you decide,
#        not so the script decides for you.
#     2. STABILITY (composition): are any clusters dominated by a single
#        sample or condition? If one sample still forms its own cluster,
#        integration has not fully removed the batch effect; if one condition
#        is over-represented that is expected and biologically relevant (we
#        kept Healthy vs Periodontitis_Post_Treatment specifically to see it).
#
#   Neither check is pass/fail. They produce ranges with interpretation, and
#   the final call is always made by eye against the Step 6.2A UMAP grid.


# --- 1. Cluster Sizes & Small-Cluster Detection ---
# ****************************************************************************#
#   The chosen partition lives in `integrated_final$seurat_clusters` (set by
#   Step 6.2B from `clusters_res_<optimal_resolution>`). First pass: how many
#   cells does each cluster contain, and does any cluster fall below 1% of
#   the dataset?
#
#   INTERPRETATION: at ~72,000 cells, 1% is roughly 720 cells. Clusters
#   smaller than that deserve a second look — they could be:
#     - a genuine rare population (e.g. a scarce immune subtype), which is
#       fine once confirmed by marker genes, or
#     - an over-split fragment of a larger population, which would argue for
#       a lower resolution.
#   The flag below is a pointer, not a verdict.
cluster_sizes <- table(integrated_final$seurat_clusters)
cluster_pct <- prop.table(cluster_sizes) * 100

cat("\nCluster sizes:\n")
print(cluster_sizes)

min_cluster_size <- 0.01 * ncol(integrated_final)
small_cluster_flag <- cluster_sizes < min_cluster_size

small_clusters <- cluster_sizes[small_cluster_flag]

if (length(small_clusters) > 0) {
    cat("\n\u26a0  Warning: Small clusters detected (<1% of cells):\n")
    print(small_clusters)
    cat("   These may represent rare cell types or over-clustering\n")
}

# Persist the size summary (cluster sizes as percentages + the flag) so the
# chosen partition's cell counts sit in the metadata ledger alongside the
# resolution-comparison table written by Step 6.2B.
cluster_size_summary <- data.frame(
    cluster = names(cluster_sizes),
    n_cells = as.integer(cluster_sizes),
    pct_cells = as.numeric(cluster_pct),
    small_cluster = as.vector(small_cluster_flag)
)
write.csv(cluster_size_summary,
    file.path(METADATA_OUT_DIR, "cluster_size_summary.csv"),
    row.names = FALSE
)
cat("→ Saved:", file.path(METADATA_OUT_DIR, "cluster_size_summary.csv"), "\n")


# --- 2. Sample Composition Across Clusters ---
# ****************************************************************************#
#   A cluster that is overwhelmingly one sample (>70% of its cells from a
#   single sample) suggests residual batch structure: integration should
#   have scattered that sample's cells across multiple clusters rather than
#   letting them pile up on their own. This is the direct follow-up to the
#   batch-effect work in Section 3/4.
#
#   INTERPRETATION: with 8 samples, a perfectly mixed cluster would have each
#   sample contribute roughly 12.5%; anything above ~70% single-sample is a
#   strong sign of incomplete correction.
sample_cluster_table <- table(integrated_final$seurat_clusters, integrated_final$sample_id)
sample_cluster_pct <- prop.table(sample_cluster_table, margin = 1) * 100

sample_specific <- apply(sample_cluster_pct, 1, max) > 70
if (any(sample_specific)) {
    cat("\n\u26a0  Warning: Sample-dominated clusters detected (>70% from single sample):\n")
    print(names(which(sample_specific)))
    cat("   These may indicate incomplete batch correction\n")
}

# Persist the sample composition matrix (clusters × samples, percentages).
sample_composition_out <- as.data.frame.matrix(sample_cluster_pct)
sample_composition_out$cluster <- rownames(sample_composition_out)
write.csv(sample_composition_out,
    file.path(METADATA_OUT_DIR, "cluster_sample_composition.csv"),
    row.names = FALSE
)
cat("→ Saved:", file.path(METADATA_OUT_DIR, "cluster_sample_composition.csv"), "\n")


# --- 3. Visualize Sample Distribution ---
# ****************************************************************************#
#   Stacked bars turn the percentages above into something readable at a
#   glance: one bar per cluster, coloured by sample. A healthy mix looks like
#   a stack with 8 roughly equal bands; a dominated cluster reads as one tall
#   band swallowing its neighbours to the right of the "mixed" threshold.
#   Colours come from the shared `sample_colors` defined in Step 3.1, so a
#   given sample is the same colour here as in every earlier figure.
sample_dist_data <- as.data.frame.matrix(sample_cluster_pct)
sample_dist_data$cluster <- rownames(sample_dist_data)
sample_dist_long <- reshape2::melt(sample_dist_data,
    id.vars = "cluster",
    variable.name = "sample", value.name = "percentage"
)

p_sample_dist <- ggplot(sample_dist_long, aes(x = cluster, y = percentage, fill = sample)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(
        title = "Sample Distribution Across Clusters",
        subtitle = "Check for sample-dominated clusters (poor integration)",
        x = "Cluster", y = "Percentage of Cells"
    ) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right"
    ) +
    scale_fill_manual(values = sample_colors)

ggsave(
    file.path(PLOTS_CLUSTERING_DIR, "07_sample_distribution_clusters.png"),
    plot = p_sample_dist, width = 12, height = 6, dpi = 300
)
cat("→ Saved:", file.path(PLOTS_CLUSTERING_DIR, "07_sample_distribution_clusters.png"), "\n")


# --- 4. Condition Composition Across Clusters ---
# ****************************************************************************#
#   The cohort keeps two conditions (Healthy, Periodontitis_Post_Treatment).
#   Unlike the sample check above, a condition-skewed cluster is NOT a
#   warning by itself — condition is exactly the biological axis we preserved
#   through integration. What we look for here is the opposite problem: have
#   the conditions been so aggressively corrected that genuinely different
#   disease-state populations got fused into a single cluster?
condition_dist <- table(integrated_final$seurat_clusters, integrated_final$condition)
condition_dist_pct <- prop.table(condition_dist, margin = 1) * 100

cat("\nCondition distribution across clusters:\n")
print(round(condition_dist_pct, 1))

# Persist the condition composition matrix alongside the sample matrix.
condition_composition_out <- as.data.frame.matrix(condition_dist_pct)
condition_composition_out$cluster <- rownames(condition_composition_out)
write.csv(condition_composition_out,
    file.path(METADATA_OUT_DIR, "cluster_condition_composition.csv"),
    row.names = FALSE
)
cat("→ Saved:", file.path(METADATA_OUT_DIR, "cluster_condition_composition.csv"), "\n")


# --- 5. Visualize Condition Distribution ---
# ****************************************************************************#
#   Side-by-side (dodged) bars per cluster, coloured by condition with the
#   shared `condition_colors` defined in Step 3.1 (same palette as Steps 3.3B
#   and 5.2A).
condition_dist_data <- as.data.frame.matrix(condition_dist_pct)
condition_dist_data$cluster <- rownames(condition_dist_data)
condition_dist_long <- reshape2::melt(condition_dist_data,
    id.vars = "cluster",
    variable.name = "condition",
    value.name = "percentage"
)

p_condition_dist <- ggplot(
    condition_dist_long,
    aes(x = cluster, y = percentage, fill = condition)
) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(
        title = "Condition Distribution Across Clusters",
        subtitle = "Check if biological conditions are preserved",
        x = "Cluster", y = "Percentage of Cells"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = condition_colors)

ggsave(
    file.path(PLOTS_CLUSTERING_DIR, "08_condition_distribution_clusters.png"),
    plot = p_condition_dist, width = 12, height = 6, dpi = 300
)
cat("→ Saved:", file.path(PLOTS_CLUSTERING_DIR, "08_condition_distribution_clusters.png"), "\n")


# ****************************************************************************#
# SUMMARY & PIPELINE MILESTONE TRANSITION
# ****************************************************************************#
# WHERE WE STARTED:
#   Step 6.2B picked the optimal resolution on silhouette score and wrote that
#   partition into `integrated_final$seurat_clusters` — well-separated in the
#   integrated embedding, but never examined as compositions of cells.
#
# WHAT WE HAVE ACCOMPLISHED:
#   We sanity-checked the chosen partition two ways: cluster sizes (flagging
#   any <1% cluster for a rare-population-vs-over-split decision) and its
#   sample/condition composition (flagging >70% single-sample clusters as
#   possible residual batch, while reading condition skew as the preserved
#   biological signal). Both the size summary and the two composition
#   matrices are persisted to the metadata ledger
#   (`cluster_size_summary.csv`, `cluster_sample_composition.csv`,
#   `cluster_condition_composition.csv`), and the Stacked/Dodged compositions
#   are saved as plots 07 and 08.
#
# WHERE WE ARE HEADING (NEXT: FINAL VISUALIZATION + SAVE):
#   With a vetted resolution and no qualifying sample-driven clusters, the
#   pipeline moves to the closing steps of this section: a single UMAP figure
#   presenting the final integrated-and-clustered object (guide §7.8, STEP
#   19), then saving that object and its metadata as the deliverable for
#   downstream cell-type annotation and differential expression (guide §8,
#   STEP 20).
# ****************************************************************************#
