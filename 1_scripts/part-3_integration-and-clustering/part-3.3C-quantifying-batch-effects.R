# ****************************************************************************#
# STEP 3.3C: Quantify batch effects with mixing metrics
# ****************************************************************************#

# TODO: Fix code headings to match file name
# TODO: Rework the transition notes to align previous and upcoming sections
#       to make a cohesive narrative for the whole analysis pipeline

# WHY THIS STEP EXISTS:
#   Step 3.3B gave us a visual impression of batch effects from four UMAP
#   panels — but at 72,649 cells, "that cluster looks mostly one color" is
#   an eyeball call, not evidence, and dense overplotting can hide or fake
#   mixing either way. This step replaces that impression with an actual
#   number per cluster: for each of the 19 clusters, what fraction of its
#   cells came from each of the 8 samples? A cluster made of real, shared
#   biology should draw fairly evenly from multiple samples; a cluster that
#   is disproportionately one sample is a signature of technical batch
#   effect rather than a distinct cell type.

# Calculate per-cluster sample composition
# `table()` cross-tabulates: rows = cluster, columns = sample_id, cell =
# raw cell count. `prop.table(..., margin = 1)` normalizes each ROW (each
# cluster) to sum to 100%, so we can compare composition across clusters of
# very different sizes on the same 0-100 scale.
cluster_composition <- table(merged_naive$seurat_clusters, merged_naive$sample_id)
cluster_composition_pct <- prop.table(cluster_composition, margin = 1) * 100

# Find clusters dominated by single samples (>50% from one sample = batch-driven)
# THRESHOLD NOTE: 50% is a simple, interpretable cutoff — "more than half of
# this cluster's cells came from just one of our 8 samples" — not a
# statistically derived value. With 8 roughly-equal-sized samples, an evenly
# mixed cluster would show ~12.5% per sample, so 50% is already a generous
# 4x-over-even bar before we flag anything. Tighten this (e.g. to 30-40%)
# for a stricter check, or loosen it if your samples are known to be very
# uneven in cell count.
dominant_sample_clusters <- apply(cluster_composition_pct, 1, max) > 50
n_batch_clusters <- sum(dominant_sample_clusters)

cat("\nBatch effect assessment:\n")
cat("  Clusters dominated by single sample (>50%):", n_batch_clusters,
    "/", nrow(cluster_composition), "\n")
# Clusters dominated by single sample (>50%): 1 / 19

if (n_batch_clusters > 0) {
  cat("  ⚠  Naive merge shows significant batch effects\n")
  cat("  → Integration is necessary\n")
}
# → Integration is necessary


# ****************************************************************************#
# SUMMARY & PIPELINE MILESTONE TRANSITION
# ****************************************************************************#
# WHERE WE STARTED:
#   Step 3.3B left us with a visual, qualitative impression from four UMAP
#   panels — persuasive, but not something we could report as a number or
#   use to make a go/no-go decision.
#
# WHAT WE HAVE ACCOMPLISHED:
#   We turned that impression into a concrete, citable statistic: 1 of 19
#   naive-merge clusters (cluster 18, per the UMAP) is dominated by a single
#   sample at >50% composition. Even a single flagged cluster confirms what
#   Step 3.3B's panel 1 suggested — cells are grouping partly by which
#   sample they came from, not purely by biology — so this naive merge is
#   not safe to analyze as-is.
#
# WHERE WE ARE HEADING (SECTION 5):
#   Because the naive merge shows measurable batch-driven clustering,
#   Section 5 introduces integration: algorithms that re-align samples in a
#   shared space based on common biological variation while suppressing the
#   technical variation this metric just caught. We won't rely on just one
#   method — Section 5 walks through and compares four different R-based
#   integration approaches (CCA, RPCA, Harmony, FastMNN) plus a note on a
#   fifth (scVI) that this tutorial doesn't cover hands-on, and this
#   3.3C composition metric becomes our benchmark: after each method we'll
#   recompute it and check whether n_batch_clusters drops from today's 1/19,
#   giving us an apples-to-apples way to judge which integration method
#   actually fixed the problem — not just a "the UMAP looks nicer" call.
# ****************************************************************************#