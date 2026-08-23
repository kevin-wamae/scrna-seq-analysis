# ****************************************************************************#
# STEP 4: Visualize batch effects in naive merge
# ****************************************************************************#


# --- THE PURPOSE OF THIS DIAGNOSTIC ---
# ****************************************************************************#
#   Step 3 gave us cluster counts and cell totals, but a number can't tell
#   you WHERE batch effects live in the data — it can't distinguish "19
#   clusters because there are 19 real cell types" from "19 clusters because
#   each sample partly formed its own island." A UMAP colored by sample makes
#   that distinction visible: if patients form isolated islands rather than
#   mixing within shared biological cell-type clusters, that's the batch
#   effect Harmony/FastMNN are meant to correct in Step 5.
#
#   IMPORTANT CAVEAT: with 72,649 cells, this UMAP is extremely dense —
#   points overplot each other, and whichever color/group gets drawn last
#   visually sits on top regardless of the true mixing underneath. A region
#   that "looks" well-mixed can just be overplotting. That's exactly why this
#   diagnostic uses four panels instead of one: no single view is trustworthy
#   alone at this scale, and even together they only give a visual
#   impression — Step 5 replaces that impression with an actual per-cluster
#   composition statistic. Treat this figure as our fixed "before" picture
#   for comparison, not as proof of anything on its own.


# --- 1. Define Color Palettes ---
# ****************************************************************************#
#   Two separate palettes, matched to two separate questions:
#     - `sample_colors`: one distinct color per patient/sample (8 total).
#       Used in p1 to check whether any single sample clusters apart from
#       the rest (batch effect) versus blending into shared clusters
#       (no batch effect).
#     - `condition_colors`: one color per biological condition (Healthy vs
#       Post-Treatment). Used in p2/p4 to check whether disease state drives
#       real separation — the signal we actually want to see once batch
#       effects are corrected — independent of which specific patient a
#       cell came from.
#
#   8 samples need 8 distinct, visible colors (Set1-style categorical
#   palette), manually assigned rather than auto-generated so colors stay
#   stable across re-runs regardless of factor level ordering.
sample_colors <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", # Healthy 1-4: Red, Blue, Green, Purple
  "#FF7F00", "#A65628", "#F781BF", "#999999"  # Post_Patient 1-4: Orange, Brown, Pink, Gray
)
names(sample_colors) <- sample_metadata$sample_id

# Colors for conditions, set to match the colors used in the source paper for
# the different conditions. If you want to change the colors, change them
# here. NOTE: the names below must exactly match the values that appear in
# `merged_naive$condition` (check with `unique(merged_naive$condition)`) — a
# mismatched name doesn't error, DimPlot just silently falls back to its
# default palette for that condition, which would quietly break the
# paper-matched coloring without any warning.
condition_colors <- c(
  "Healthy" = "#2E86AB",                     # Blue
  "Periodontitis_Post_Treatment" = "#F18F01" # Orange
)


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