# ****************************************************************************#
# STEP 5.1B: Calculate integration quality metrics
# ****************************************************************************#

# --- WHY THIS STEP EXISTS ---
# ****************************************************************************#
#   Step 3.3C's per-cluster composition statistic only ever looked at the
#   naive merge. Now that we have four independently integrated objects
#   (CCA, RPCA, Harmony, FastMNN) sitting alongside that naive baseline, we
#   need one comparable number per method to objectively say which
#   correction actually worked best — a visual UMAP impression (Step 5.1A)
#   isn't citable evidence on its own.
#
#   QUANTITATIVE COMPARISON: MIXING AND SEPARATION METRICS
#   To objectively compare integration methods, we calculate mixing scores
#   that quantify how well cells from different samples mix within
#   biological clusters. Higher scores indicate better integration.
#
#   HOW THE MIXING METRIC WORKS: for each cell, we find its k = 50 nearest
#   neighbors in the UMAP space and calculate the diversity of samples
#   among those neighbors using the inverse Simpson's Index. If all
#   neighbors are from different samples (perfect mixing), the score is
#   high. If all neighbors are from the same sample (poor integration),
#   the score is low. Averaging this over all cells gives one mixing score
#   per method.


# --- 1. Helper: Local Inverse Simpson's Index for One Cell ---
# ****************************************************************************#
#   Small enough to test in isolation and to vectorize with purrr::map_dbl
#   over neighborhoods, rather than burying the per-cell math inside the
#   outer sapply() loop.
inverse_simpson <- function(group_labels) {
    props <- table(group_labels) / length(group_labels)
    1 / sum(props^2)
}


# --- 2. Helper: Mixing Metric for One Seurat Object ---
# ****************************************************************************#
#   Measures sample diversity in each cell's k-nearest-neighborhood on a
#   given 2D embedding, then averages across all cells to give one score
#   for the whole object. `FNN::get.knn()` does the neighbor search;
#   everything after that is purrr piping over the resulting index matrix.
calculate_mixing_metric <- function(seurat_obj, group_by = "sample_id",
                                    reduction = "umap", k = 50) {
    # Get the 2D embedding coordinates for the requested reduction
    embedding <- Embeddings(seurat_obj, reduction = reduction)[, 1:2]

    # Find each cell's k nearest neighbors
    nn_index <- FNN::get.knn(embedding, k = k)$nn.index

    # Pull the grouping variable (e.g. sample_id) once, up front
    group_labels <- seurat_obj@meta.data[[group_by]]

    # For each cell (each row of nn_index), look up its neighbors' sample
    # labels and score how mixed that neighborhood is, then average across
    # all cells for one summary score
    nn_index %>%
        asplit(1) %>%
        map_dbl(~ inverse_simpson(group_labels[.x])) %>%
        mean()
}


# --- 3. Define the Methods to Compare ---
# ****************************************************************************#
#   One entry per checkpointed object from Steps 3.3A / 4.2A-D, paired with
#   the UMAP reduction each method produced. Keeping this as a single named
#   list is what lets the rest of the step iterate with purrr instead of
#   repeating the same calculation five times by hand.
methods_list <- list(
    "Naive"   = list(obj = merged_naive,       reduction = "umap"),
    "CCA"     = list(obj = integrated_cca,     reduction = "umap.cca"),
    "RPCA"    = list(obj = integrated_rpca,    reduction = "umap.rpca"),
    "Harmony" = list(obj = integrated_harmony, reduction = "umap.harmony"),
    "FastMNN" = list(obj = integrated_fastmnn, reduction = "umap.mnn")
)


# --- 4. Calculate Mixing Score and Cluster Count for Each Method ---
# ****************************************************************************#
#   `map_dfr()` iterates the helper over every method and row-binds the
#   results straight into a tidy data frame, replacing the two parallel
#   sapply() calls (one for mixing_score, one for n_clusters) with a single
#   pass that keeps each method's two metrics together.
mixing_results <- LOG_STEP("Calculating integration mixing scores across all methods...", {
    methods_list %>%
        imap_dfr(function(method_info, method_name) {
            tibble(
                method = method_name,
                mixing_score = calculate_mixing_metric(
                    method_info$obj,
                    group_by  = "sample_id",
                    reduction = method_info$reduction
                ),
                n_clusters = method_info$obj$seurat_clusters %>%
                    unique() %>%
                    length()
            )
        })
})

# Display results
cat("\nIntegration Quality Metrics:\n")
cat("(Higher mixing score = better sample mixing)\n\n")
print(mixing_results)


# --- 5. Save Metrics ---
# ****************************************************************************#
#   Uses METADATA_OUT_DIR, registered back in Step 3.1's dynamic output-path
#   setup, rather than hardcoding the date-stamped job run path here.
write.csv(mixing_results,
    file.path(METADATA_OUT_DIR, "integration_comparison_metrics.csv"),
    row.names = FALSE
)
cat("→ Saved:", file.path(METADATA_OUT_DIR, "integration_comparison_metrics.csv"), "\n")


# --- 6. Visualize Mixing Scores ---
# ****************************************************************************#
#   One bar per method, ordered worst-to-best left-to-right, with the exact
#   score labeled above each bar so the figure can be read without cross-
#   referencing the CSV.
p_mixing <- mixing_results %>%
    ggplot(aes(x = reorder(method, -mixing_score), y = mixing_score, fill = method)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = round(mixing_score, 2)), vjust = -0.5) +
    labs(
        title = "Integration Method Comparison: Sample Mixing",
        subtitle = "Higher score = better integration (samples mix within cell types)",
        x = "Integration Method",
        y = "Mixing Score (Local Inverse Simpson's Index)"
    ) +
    theme(legend.position = "none") +
    scale_fill_brewer(palette = "Set2")

# Save mixing scores image
ggsave(file.path(PLOTS_COMPARISON_DIR, "04_mixing_scores.png"),
    p_mixing, width = 10, height = 6, dpi = 300
)
cat("→ Saved:", file.path(PLOTS_COMPARISON_DIR, "04_mixing_scores.png"), "\n\n")


# ****************************************************************************#
# SUMMARY & PIPELINE MILESTONE TRANSITION
# ****************************************************************************#
# WHERE WE STARTED:
#   Step 3.3C's per-cluster composition statistic only ever looked at the
#   naive merge in isolation — no way to compare it against the four
#   integration methods on the same numeric scale.
#
# WHAT WE HAVE ACCOMPLISHED:
#   We computed a local inverse Simpson's Index mixing score for all five
#   objects (naive baseline + CCA, RPCA, Harmony, FastMNN), giving each
#   method one comparable number for how well it interleaves samples within
#   each cell's k = 50 nearest neighbors. Results are saved to
#   integration_comparison_metrics.csv and visualized in
#   04_mixing_scores.png, letting us rank methods by mixing quality instead
#   of eyeballing UMAP panels.
#
# WHERE WE ARE HEADING (STEP 13 — ASSESSING BIOLOGICAL SIGNAL PRESERVATION):
#   A high mixing score alone isn't the whole story: integration is a
#   balancing act between removing technical variation (samples mixing)
#   and preserving real biological variation (Healthy vs Post-Treatment
#   remaining distinguishable). A method that mixes everything together
#   indiscriminately may have over-integrated, blending away genuine
#   disease signal along with the batch noise. Step 13 complements this
#   step's mixing score with a condition-separation metric — the distance
#   between each condition's mean UMAP position — so we can look for the
#   combination that actually indicates good integration: high mixing
#   *and* moderate-to-high condition separation, not just one or the
#   other.
# ****************************************************************************#