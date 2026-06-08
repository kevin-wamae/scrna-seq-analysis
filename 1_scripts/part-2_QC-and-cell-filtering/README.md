# scRNA-seq: Quality Control & Cell Filtering

A step-by-step, heavily-commented quality control (QC) workflow for single-cell
RNA-seq data from the **10x Genomics** platform, implemented in **R** with
**Seurat v5** and complementary Bioconductor packages.

The pipeline takes a raw Cell Ranger count matrix and produces a clean Seurat
object ready for integration, clustering, cell-type annotation, and downstream
biology. It is written to be *read* as much as run: every script explains the
biology and the statistics behind each decision, and frames its outputs as
**ranges with interpretation** rather than fixed pass/fail numbers — so the
reasoning transfers to any tissue, not just the dataset you start with.

> This is the QC stage of a larger series. An upstream step covers processing
> raw FASTQ files through Cell Ranger to generate the count matrices used here.

---

## Philosophy: ranges, not magic numbers

There is no universal "good" UMI count or "correct" mitochondrial cutoff. A
resting lymphocyte and an activated macrophage are both healthy cells with
wildly different profiles; a clean PBMC prep and a frozen tumor biopsy have
different baselines for almost every metric. So this workflow does two things
everywhere:

1. **Reports a metric, then tells you the range it usually falls in and what
   the edges of that range mean** — so you can judge *your* sample on its own
   terms.
2. **Prefers visual inspection and data-driven thresholds over hard-coded
   constants.** You look at the distribution, find the natural break, and set
   the cutoff there.

Treat every number in this README as a *guide to interpretation*, not a rule to
copy.

---

## What you'll learn

- Load and inspect a raw 10x matrix, and reason about sparsity and UMI distributions
- Separate real cells from empty droplets statistically with **EmptyDrops** (`DropletUtils`)
- Estimate and *conditionally* remove ambient RNA with **SoupX**
- Detect and remove doublets with **scDblFinder**
- Set biologically-informed, visually-justified cell-level QC thresholds
- Filter low-information and contaminant genes
- Normalize and select highly variable features
- Produce a reproducible, auditable QC-filtered dataset

---

## Pipeline at a glance

```text
Raw 10x matrix (per sample)
        │
   0  Install packages
   1  Load libraries & configure environment
   2  Load 10x Genomics data (raw vs filtered decision)
   3  Initial data exploration (sparsity, UMI distribution)
   4  Empty droplet detection ........ DropletUtils / EmptyDrops
   5  Ambient RNA correction ......... SoupX        (conditional, skip if <5%)
   6  Doublet detection .............. scDblFinder
   7  Cell-level QC (manual, visual threshold setting)
   8  Gene-level QC (rare genes + hemoglobin removal)
   9  Normalization & variable feature selection
  10  Save clean dataset + QC summary
        │
   Clean, QC-filtered sample → ready for integration & clustering
```

Steps 4 and 5 only apply when you load the **raw** matrix. If you load Cell
Ranger's **filtered** matrix instead, skip 4 and 5 and start at step 6 — Cell
Ranger has already done basic empty-droplet removal and ambient correction.

> **Run QC independently per sample, then integrate.** Different samples have
> different quality profiles and need their own inspection and thresholds. Keep
> the *matrix choice* (raw vs filtered) and the *general approach* consistent
> across samples — but do not blindly reuse numeric thresholds between them.

---

## Repository structure

```text
QC-and-cell-filtering/
│
├── 0-install-packages.R
├── 1-load-libraries-and-configuration.R
├── 2-load-10x-genomics-data.R
├── 3-initial-data-exploration.R
├── 4-correction-empty-droplet.R
├── 5-correction-ambient-RNA.R
├── 6-correction-doublets.R
├── 7-correction-qc_visualization_filtering.R
├── 8-correction-gene_level_QC.R
├── 9-normalisation-and-variable-feature-selection.R
├── 10-save-clean-dataset.R
│
├── plots/          # diagnostic figures (created on first run)
├── qc_metrics/     # QC_summary.csv
└── filtered_data/  # clean .rds + cell_metadata.csv
```

> **Output folders.** Script 1 leaves `setwd()` and the `dir.create()` calls
> commented out and recommends working inside an **RStudio Project** (`.Rproj`)
> so paths stay relative and portable. If you are *not* using an RStudio
> Project, create `plots/`, `qc_metrics/`, and `filtered_data/` yourself (or
> uncomment those lines) before running, or `ggsave()` / `write.csv()` will fail.

---

## Prerequisites

| Requirement | Recommended             |
| ----------- | ----------------------- |
| R           | ≥ 4.3.0                |
| RAM         | ≥ 16 GB (8 GB minimum) |
| Storage     | ≥ 2 GB per sample      |
| OS          | Linux, macOS, Windows   |

### Packages

| Package                       | Source       | Role                                             |
| ----------------------------- | ------------ | ------------------------------------------------ |
| `Seurat` / `SeuratObject` | CRAN         | Core scRNA-seq framework & object model          |
| `DropletUtils`              | Bioconductor | Empty droplet detection (EmptyDrops)             |
| `SoupX`                     | CRAN         | Ambient RNA correction                           |
| `scDblFinder`               | Bioconductor | Doublet detection                                |
| `SingleCellExperiment`      | Bioconductor | Object model used by the Bioc tools              |
| `ggplot2` / `patchwork`   | CRAN         | Plotting and panel composition                   |
| `dplyr`                     | CRAN         | Fast metadata manipulation                       |
| `colorout`                  | r-multiverse | Colorized console errors/warnings (optional QoL) |

Install everything (idempotent — only fetches what's missing):

```r
source("0-install-packages.R")
```

`colorout` is a quality-of-life package from a community repository; if that
repo is unreachable, remove it from the install script and the `library()` call
in script 1 — nothing else depends on it.

---

## Pipeline reference

Each step below pairs *what the script does* with *how to read its output*.

### Step 0 — Install packages

Pulls dependencies from three ecosystems (CRAN, Bioconductor, community) behind
`requireNamespace()` guards so nothing is re-downloaded on repeat runs.

### Step 1 — Libraries & configuration

Loads libraries, pins reproducibility with `set.seed(100)`, and verifies a
Seurat ≥ 5.0.0 install. Explains the Seurat 5 *layer architecture* (raw and
corrected counts can coexist in separate layers of one object) and recommends
RStudio Projects over `setwd()` for portable paths.

### Step 2 — Load 10x data

Reads Cell Ranger output with `Read10X()` and builds the object with
`min.cells = 0, min.features = 0` — **no early filtering**, because EmptyDrops
and SoupX need the full ambient background intact. Sample metadata is attached
in a single `dplyr::mutate()`.

> **Edit before running:** set the Cell Ranger output path, and choose your
> matrix. The script defaults to the **raw** matrix (Option A); the **filtered**
> matrix (Option B) is provided commented-out.

| Matrix             | Contents                                       | When to use                                |
| ------------------ | ---------------------------------------------- | ------------------------------------------ |
| **raw**      | All droplets (cells + empty + soup + doublets) | Full control; run steps 4–5 (the default) |
| **filtered** | Only Cell Ranger-called cells                  | Faster;**skip steps 4–5**           |

### Step 3 — Initial exploration

Reports total droplets/genes, matrix sparsity, median UMI per droplet, and a
count of high-UMI droplets.

| Metric                              | Typical range | What it tells you                                                                                                                                |
| ----------------------------------- | ------------- | ------------------------------------------------------------------------------------------------------------------------------------------------ |
| Matrix sparsity                     | 95–99.9%     | The hallmark of scRNA-seq (biological + dropout zeros).**Below ~95%** warns of ambient-soup flooding or large doublet clumps.              |
| Median UMI per droplet (raw matrix) | ~1            | Expected and*not* a failure — empty droplets vastly outnumber cells and drag the median down. The real cells live in a second, high-UMI mode. |
| High-UMI droplet count              | varies        | A quick peek above the noise floor at your likely cellular yield.**Sanity check only — never a hard filter.**                             |

### Step 4 — Empty droplet detection

Runs `emptyDrops()` (`lower = 100`, `niters = 10000`, `test.ambient = TRUE`),
calls cells at **FDR < 0.01**, and treats `NA` as empty. The key idea is that it
tests each droplet's *expression profile* against the ambient soup, not just its
total count — so a low-UMI droplet with a cell-type-specific gene mix is rescued
as a real cell, while a low-UMI droplet of generic housekeeping genes is dropped.
Saves `plots/01_empty_droplets.png`.

| Reading the result                                                    | Interpretation                                                                             |
| --------------------------------------------------------------------- | ------------------------------------------------------------------------------------------ |
| Cells called as % of droplets                                         | Commonly ~0.5–2% of all raw droplets — most captures are empty.                          |
| Two clear modes on the `log10(UMI)` histogram with a valley between | High-quality prep: background soup cleanly separated from real cells.                      |
| No valley / heavy overlap                                             | Possible degradation or heavy ambient contamination — inspect before trusting downstream. |

> **Critical:** the full raw matrix is cached as `raw_counts_all_droplets`
> *before* subsetting — SoupX needs it in step 5.

### Step 5 — Ambient RNA correction (conditional)

Builds a `SoupChannel` from the table-of-droplets (raw) and table-of-cells
(validated), does a quick clustering to give SoupX cell identities, then
estimates contamination with `autoEstCont()` inside a `tryCatch`. Correction is
applied **only when warranted**.

| Contamination fraction | Interpretation                                     | Action                                                                         |
| ---------------------- | -------------------------------------------------- | ------------------------------------------------------------------------------ |
| `NULL` or < 5%       | Clean prep, or no clear ambient signature to model | **Skip** — forcing a correction risks eroding real low-abundance signal |
| 5–15%                 | Standard run                                       | Correct to sharpen markers                                                     |
| > 20%                  | Heavy cell lysis (fibrous tumor, frozen brain)     | Correction essential                                                           |

A `NULL` result is common and usually *good news* — it means the estimator
couldn't find enough ambient signal to model, typical of gentle dissociation or
samples lacking strong anchor markers (hemoglobin, myelin).

### Step 6 — Doublet detection

`JoinLayers()`, convert to SCE, run `scDblFinder(dbr = NULL)` (harmless warnings
suppressed), write `doublet_class` / `doublet_score` to metadata, build a quick
PCA/UMAP (native `uwot`, no Python needed), save `plots/02_doublets_umap.png`,
and keep singlets.

Doublet rate scales with cell loading:

| Cells recovered | Expected doublet rate |
| --------------- | --------------------- |
| ~1,000          | ~1%                   |
| ~5,000          | 3–5%                 |
| ~10,000         | 7–10%                |

So a ~10% rate is acceptable on a large run but a red flag on a small one. On
the UMAP, watch for: a separate island of doublets (a phantom "cell type"), an
entire real cluster flagged red (over-sensitive detection — lower `dbr`), or a
red halo around every cluster (uncorrected ambient contamination from step 5).

### Step 7 — Cell-level QC (manual, visual)

Computes `percent.mt` (`^MT-`), `percent.ribo` (`^RP[SL]`), and
`log10GenesPerUMI`; prints distribution summaries; and saves `03_qc_violins.png`,
`04_qc_scatter.png`, and **both** `05_filtering_thresholds_log.png` and
`05_filtering_thresholds_linear.png`. You set thresholds *after* reading the
plots.

| Metric                   | Healthy range                                     | What the extremes mean                                                                                                                                              |
| ------------------------ | ------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `nCount_RNA` (UMIs)    | median > ~3,000 is robust; < ~1,000 is concerning | Sequencing depth. Very low = under-sequencing, loading failure, or degradation.                                                                                     |
| `nFeature_RNA` (genes) | median > ~1,200 is robust; < ~500 is concerning   | Transcript diversity. Very low = debris or empty-ish droplets.                                                                                                      |
| `percent.mt`           | typically 1–5% in healthy cells                  | High = ruptured/dying cells (cytoplasmic RNA leaks, mitochondria stay trapped). Use the**95th percentile** to anchor a cap that trims only the necrotic tail. |
| `percent.ribo`         | ~10–30% in active cells                          | Cross-validation signal:`ribo → 0` while `mt` spikes marks a dead "ghost cell."                                                                                |
| `log10GenesPerUMI`     | higher = more complex                             | Low values flag debris/low-complexity droplets.                                                                                                                     |

How to set each cutoff visually: put the **minimums** just above the low-complexity
"shelf" at the bottom-left of the UMI-vs-genes scatter; put the **MT cap** at the
elbow where the dot cloud turns vertical toward 100% (≈10–15% for fluid preps,
up to ~20% for tough tumors); put the **maximums** where the main cloud thins
into sparse high outliers (residual doublets).

A built-in guard reports the removal rate and warns at the edges:

| Removal rate | Verdict                                             |
| ------------ | --------------------------------------------------- |
| 10–25%      | Healthy sweet spot                                  |
| < 5%         | Too lenient — likely retaining debris/doublets     |
| 30–40%      | Too strict — usually the MT% cap; relax it         |
| > 50%        | Catastrophic — sample quality or threshold problem |

> **Edit before running:** the five threshold values
> (`nfeature_min/max`, `ncount_min/max`, `mt_thresh`). They are sample-specific.

### Step 8 — Gene-level QC

Builds a per-gene detection table, plots `06_gene_detection.png`, then keeps
genes detected in **≥ `min_pct_cells`** of cells **and** drops **hemoglobin
genes** (`^HB[AB]`, red-blood-cell contamination). MT and ribosomal genes are
*flagged but retained*. The percentage-based floor scales with dataset size,
unlike a fixed "≥3 cells" rule.

| Detection floor   | Strategy     | Best for                                                        |
| ----------------- | ------------ | --------------------------------------------------------------- |
| 0.1% of cells     | Lenient      | Heterogeneous tissue / rare cell types — protects rare markers |
| 1% of cells       | Standard     | Uniform populations — strips background, speeds compute        |
| ≥3 cells (fixed) | Conservative | Removing singleton artifacts / mapping errors                   |

### Step 9 — Normalization & variable features

`NormalizeData(method = "LogNormalize", scale.factor = 10000)` then
`FindVariableFeatures(method = "vst", nfeatures = 2000)`. Prints the top variable
genes and saves `07_variable_features.png`. On that plot, the dense black shelf
is the non-variable background (noise + housekeeping); the red points are the
~2,000 features whose *standardized* variance drives cell-type structure.

### Step 10 — Save clean dataset

Writes the final object plus a single comprehensive audit trail:

```text
filtered_data/<sample>_qc_filtered.rds   # raw + normalized layers, metadata, HVGs
filtered_data/cell_metadata.csv          # per-cell metrics
qc_metrics/QC_summary.csv                # one comprehensive ledger
```

`QC_summary.csv` tracks cell survival at every gate
(`initial_droplets → after_emptydrops → after_soupx → after_doublets → after_cell_qc → final_cells`), gene attrition
(`initial_genes → after_frequency_filter → final_genes`), per-stage removal
percentages, final median UMI/genes/MT%, contamination and doublet rates, and
the exact thresholds used — everything a methods section needs.

Expected attrition, for context: gene removal of **30–50%** at the frequency
filter is normal (many features are unexpressed or dropouts), while the
hemoglobin step usually removes **<1%**.

---

## Outputs

```text
plots/
├── 01_empty_droplets.png
├── 02_doublets_umap.png
├── 03_qc_violins.png
├── 04_qc_scatter.png
├── 05_filtering_thresholds_log.png
├── 05_filtering_thresholds_linear.png
├── 06_gene_detection.png
└── 07_variable_features.png

qc_metrics/QC_summary.csv
filtered_data/<sample>_qc_filtered.rds
filtered_data/cell_metadata.csv
```

---

## Reference ranges (guides, not rules)

| Metric                | Typical range | At the edges                                       |
| --------------------- | ------------- | -------------------------------------------------- |
| Final cells           | 3,000–10,000 | Depends entirely on loading and tissue             |
| Median UMI            | 2,000–10,000 | Low = shallow depth; very high = possible doublets |
| Median genes / cell   | 1,000–3,000  | Low = debris-heavy; tissue-dependent               |
| Median MT%            | < 3–5%       | High = dissociation stress / cell death            |
| Ambient contamination | < 5%          | High in lysis-prone tissue (tumor, brain)          |
| Doublet rate          | 1–10%        | Scales with cell loading                           |
| Cell removal rate     | 10–25%       | < 5% too lenient, > 30% too strict                 |

---

## Troubleshooting

| Symptom                            | Likely cause / fix                                                                                                                            |
| ---------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------- |
| Very few cells after filtering     | Thresholds too strict, or poor sample quality. Relax the MT cap for your tissue; check Cell Ranger `web_summary.html`.                      |
| > 30% cells removed                | Usually the MT% cap. Inspect the removal breakdown; if one gate accounts for >90% of losses, that threshold is wrong.                         |
| Doublet rate high for a small prep | Over-loading or an over-sensitive detector. Optionally threshold on `scDblFinder.score`.                                                    |
| A cluster driven by QC metrics     | Incomplete filtering*or* a real stressed population (check HSP genes). Tighten if technical; annotate if biological.                        |
| SoupX returns `NULL`             | Often**good** — the sample is clean. The `tryCatch` proceeds on original counts. Only set contamination manually with strong reason. |

---

## Next steps

Per-sample QC complete → integration → scaling → PCA → UMAP → clustering →
cell-type annotation → differential expression → pathway enrichment.

---

## Guiding principle

> QC is not about maximizing the number of cells retained — it's about
> maximizing biological signal while minimizing technical noise. Every filtering
> decision should be driven by evidence from *your* data, not by thresholds
> copied from another dataset.

---

## References

Hao et al. 2024, *Nat Biotechnol* (Seurat 5) · Lun et al. 2019, *Genome Biol*
(EmptyDrops) · Young & Behjati 2020, *GigaScience* (SoupX) · Germain et al.
2021, *F1000Res* (scDblFinder) · Luecken & Theis 2019, *Mol Syst Biol*
(best practices).

Docs: [Seurat](https://satijalab.org/seurat/) ·
[DropletUtils](https://github.com/MarioniLab/DropletUtils) ·
[SoupX](https://github.com/constantAmateur/SoupX) ·
[scDblFinder](https://github.com/plger/scDblFinder)
