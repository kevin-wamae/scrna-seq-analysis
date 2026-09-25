# ****************************************************************************#
# STEP 3.3B: Visualize batch effects in naive merge
# ****************************************************************************#


# --- PREREQUISITES (scripts that must run before this one in this session) ---
# ****************************************************************************#
#   • part-3.0-install-packages.R — packages for library() calls.
#   • part-3.1-load-libraries-and-configuration.R — RUN_ID, LOG_STEP,
#     output-directory variables, and the shared `sample_colors` /
#     `condition_colors` palettes.
#   • part-3.2-load-10x-quality-controlled-data.R — supplies `sample_metadata`
#     and the sample objects merged in Step 3.3A.
#   • part-3.3A-merge-naive.R — supplies `merged_naive`, the object plotted
#     here.
#   The dependency manager is always re-sourced here (base-R only, cheap) so
#   the latest DEP_TABLE is loaded on every run; ensure_dependencies() does
#   nothing if those objects are already in the environment, otherwise it
#   offers to source them (interactive) or stops with a clear message
#   (non-interactive/batch). When everything is already satisfied, an
#   interactive session is asked whether to re-source fresh or proceed as-is.
source("1_scripts/part-3_integration-and-clustering/part-3.0-dependencies.R")
ensure_dependencies(step = "part-3.3B-visualize-batch-effects.R")


# --- THE PURPOSE OF THIS DIAGNOSTIC ---
# ****************************************************************************#
#   The naive merge gave us cluster counts and cell totals, but a number can't tell
#   you WHERE batch effects live in the data — it can't distinguish "19
#   clusters because there are 19 real cell types" from "19 clusters because
#   each sample partly formed its own island." A UMAP colored by sample makes
#   that distinction visible: if patients form isolated islands rather than
#   mixing within shared biological cell-type clusters, that's the batch
#   effect Harmony/FastMNN are meant to correct.
#
#   IMPORTANT CAVEAT: with 72,649 cells, this UMAP is extremely dense —
#   points overplot each other, and whichever color/group gets drawn last
#   visually sits on top regardless of the true mixing underneath. A region
#   that "looks" well-mixed can just be overplotting. That's exactly why this
#   diagnostic uses four panels instead of one: no single view is trustworthy
#   alone at this scale, and even together they only give a visual
#   impression — the quantitative step replaces that impression with an actual per-cluster
#   composition statistic. Treat this figure as our fixed "before" picture
#   for comparison, not as proof of anything on its own.


# --- 1. Color Palettes (Shared) ---
# ****************************************************************************#
#   `sample_colors` and `condition_colors` are defined once in Step 3.1 and
#   reused here, so this "before" figure and every later figure colour the
#   same sample/condition identically:
#     - `sample_colors`: one distinct color per sample (8 total), used in p1
#       to check whether any single sample clusters apart from the rest
#       (batch effect) versus blending into shared clusters (no batch effect).
#     - `condition_colors`: one color per biological condition, used in p2/p4
#       to check whether disease state drives real separation — the signal we
#       actually want to see once batch effects are corrected.
#
#   Panel p3 colours by cluster and keeps the ggplot defaults, so no palette
#   is needed for it.


# --- 2. Build Diagnostic UMAP Panels ---
# ****************************************************************************#
#   Four complementary views of the same embedding. Each one is only a
#   partial answer on its own — together they cover for each other's blind
#   spots (especially the overplotting problem noted above):
#     - p1: Are samples separating out on their own? This is the actual
#       batch-effect detector — look for a sample forming its own isolated
#       island rather than blending with the other 7.
#     - p2: Is biological condition driving separation instead? Collapsing
#       8 samples down to 2 conditions can visually flatten a real batch
#       effect (the majority group just paints over the minority group), so
#       don't read this panel as "no batch effect" on its own — cross-check
#       against p1 and p4.
#     - p3: What do the naive (uncorrected) cluster boundaries look like?
#       A plain reference for "these are the 19 groups we're trying to
#       explain," with no sample/condition coloring yet.
#     - p4: Faceted per-sample, one small subplot per sample on shared axes.
#       This is the fix for overplotting: a sample-specific blob that's
#       invisible in p1's crowded overlay (buried under other colors) will
#       stand out clearly here as present in one facet and absent from the
#       rest.

# UMAP colored by sample - shows batch effects
p1_naive <- DimPlot(merged_naive, reduction = "umap", group.by = "sample_id",
                    pt.size = 0.05, cols = sample_colors) +
  ggtitle("Naive Merge: Colored by Sample") +
  theme(legend.position = "right", legend.text = element_text(size = 8))

# UMAP colored by condition
p2_naive <- DimPlot(merged_naive, reduction = "umap", group.by = "condition",
                    pt.size = 0.05, cols = condition_colors) +
  ggtitle("Naive Merge: Colored by Condition")

# UMAP colored by clusters
p3_naive <- DimPlot(merged_naive, reduction = "umap", group.by = "seurat_clusters",
                    pt.size = 0.05, label = TRUE, label.size = 5) +
  ggtitle("Naive Merge: Clusters") +
  NoLegend()

# Split by sample to see separation
p4_naive <- DimPlot(merged_naive, reduction = "umap", group.by = "condition",
                    split.by = "sample_id", pt.size = 0.05, ncol = 4,
                    cols = condition_colors) +
  ggtitle("Naive Merge: Split by Sample") +
  theme(strip.text = element_text(size = 9, face = "bold"))


# --- 3. Assemble & Save Diagnostic Panel ---
# ****************************************************************************#
#   `patchwork`'s `|` (side-by-side) and `/` (stacked) operators combine all
#   four panels into one 2x2 figure, so sample-level, condition-level,
#   cluster-level, and per-sample views can all be compared at a glance in a
#   single saved file rather than four separate ones.
combined_naive <- (p1_naive | p2_naive) / (p3_naive | p4_naive)

ggsave(
  file.path(PLOTS_OUT_DIR, "01_naive_merge_batch_effects.png"),
  plot = combined_naive, width = 16, height = 12, dpi = 300
)

cat("→ EXAMINE file:", file.path(PLOTS_OUT_DIR, "01_naive_merge_batch_effects.png"), "\n")
cat("   Look for samples forming isolated islands (batch effect) vs\n")
cat("   blending into shared clusters (no strong batch effect).\n")
cat("   Remember: with this many cells, a visual read is only a first\n")
cat("   pass — treat p1/p4 as the primary evidence and p2 with caution.\n\n")


# ****************************************************************************#
# SUMMARY & PIPELINE MILESTONE TRANSITION
# ****************************************************************************#
# WHERE WE STARTED:
#   Step 3.3A gave us a merged-but-uncorrected object and a bare number: 19
#   clusters. That number alone can't distinguish real biology from batch
#   effects — a cluster count says nothing about WHY cells grouped the way
#   they did.
#
# WHAT WE HAVE ACCOMPLISHED:
#   We rendered the naive merge as a four-panel diagnostic UMAP — colored by
#   sample, by condition, by cluster identity, and split per-sample — and
#   saved it to 01_naive_merge_batch_effects.png. This turns "batch effects
#   might be present" into something directly visible: wherever a specific
#   sample or patient forms an isolated island rather than blending into a
#   shared cluster, that's technical variation masquerading as biological
#   signal. This figure is now our fixed "before" reference — every
#   integration method's UMAP in later steps gets compared back against it.
#
# WHERE WE ARE HEADING (STEP 3.3C):
#   A visual impression of batch effects is persuasive but not quantitative
#   — "that cluster looks mostly red" isn't a number you can report or set a
#   threshold against, and at this cell count it can even be an artifact of
#   overplotting rather than a real pattern. In the next step, we will calculate
#   the actual per-cluster sample composition and flag any cluster dominated
#   by a single sample, turning what Step 3.3A showed us into a concrete,
#   citable statistic.
# ****************************************************************************#