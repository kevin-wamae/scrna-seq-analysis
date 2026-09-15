# scrna-seq-analysis

Single-cell RNA-seq analysis pipeline (10x Genomics, Seurat v5) for a
periodontal-disease cohort (GSE174609). From raw FASTQ (NCBI SRA) through
Cell Ranger counting, per-sample QC, batch-effect integration, best-method
selection, and clustering at multiple resolutions. Written to be *read as
much as run*: every script explains the biology and statistics, frames
outputs as **ranges with interpretation** rather than pass/fail numbers,
and ends with a milestone-transition summary.

## Repository layout
- `pixi.toml` — conda env source of truth. Features: `part1-feature`,
  `part2-feature`, `part3-feature`, all on `core-feature` (Seurat, dplyr,
  ggplot2, dev tooling). Envs: `part1`, `part2`, `part3`. Task `srun` opens
  an interactive HPC session: `srun --ntasks=1 --cpus-per-task=16 --mem=50G --time=02:00:00 --pty bash`.
- `1_scripts/part-0_pixi-setup/` — `validate_environment.sh` checks every env's tools.
- `1_scripts/part-1_fastq-to-counts/` — bash: SRA download → FastQC/MultiQC → cellranger count → copy.
- `1_scripts/part-2_QC-and-cell-filtering/` — R, per-sample: empty droplets →
  ambient RNA (SoupX) → doublets (scDblFinder) → QC viz/filter → gene-level
  QC → normalize/HVG → save clean `.rds` (`part-2.0`…`part-2.10`).
- `1_scripts/part-3_integration-and-clustering/` — R: dependency manager +
  install → merge naive baseline → 4 integration methods (CCA, RPCA, Harmony,
  FastMNN) → compare visually + quantitatively → pick best (`part-5.3B`) →
  cluster at multiple resolutions (`part-6.1`). Scripts: `part-3.0-install-packages.R`,
  `part-3.0-dependencies.R`, `part-3.1`…`part-5.3B`, `part-6.1-clustering-multiple-resolutions.R`.
- `2_input/` — `cellranger-matrix-counts/`, and `sample-metadata/sample_names.tsv`
  (the **single source of truth** for sample_id/condition/patient_id/QC thresholds).
- `3_output/<RUN_ID>/` — `qc_and_filtering/` (plots, filtered_data, metrics) and
  `integration_and_clustering/` (plots, integrated_data, metadata, checkpoints).
- `logs/` — HPC job logs.

## Cohort & analysis design
- 12 samples total; downstream analysis keeps **8**: Healthy vs
  Periodontitis_Post_Treatment (`Periodontitis_Pre_Treatment` filtered out in
  part-3.2 and elsewhere).
- `sample_names.tsv` columns: sample_id, sra_id, condition, patient_id,
  nfeature_min/max, ncount_min/max, mt_thresh.

## How to run
- Validate env: `bash 1_scripts/part-0_pixi-setup/validate_environment.sh`
- Part 1: `pixi run -e part1 <tool>`; shell scripts in part-order.
- Parts 2–3: `pixi run srun` → interactive R on the HPC node, then `source()`
  scripts in order. On macOS/locally, use the relevant pixi env's R
  (RStudio / `pixi run -e part3 Rscript ...`).
- Working directory must be the **repo root** (scripts use relative paths like
  `2_input/...`, `3_output/...`).
- The dependency manager (below) makes "source in any order" safe for every
  step that has a `DEP_TABLE` entry + top-of-script `ensure_dependencies()` call.

## Pipeline conventions & config
- `RUN_ID <- "2026_06_09_brown_job_3058993"` — date-stamped job id shared by
  parts 2 and 3; every output nests under `3_output/<RUN_ID>/`. Never re-hardcode.
- `CHECKPOINT_FORMAT` — 1 = qs2 (fast, multithreaded; recommended), 2 = rds.
  Applies to every checkpoint read/write.
- `set.seed(100)` early, once, before any stochastic step (Harmony, UMAP, etc.).
- Worker count: `N_WORKERS` from `SLURM_CPUS_PER_TASK` (minus 1) else
  `min(4, cores-1)`. `future.globals.maxSize` sized from SLURM mem allocation,
  else interactive prompt (default 4 GB), else silent 4 GB.

## Coding conventions (follow these)
**1. Naming** — scripts: `part-<part>.<step><letter>-<slug>.R`; objects
snake_case; steps referenced as e.g. "Step 5.2B". (Note: `part-6.1`'s internal
banner says "STEP 15" — historically inconsistent, left as-is for now.)

**2. Comment style**
- File banner: `# ***` line, `# STEP <N>: Title`, closing `# ***`.
- Section headings: `# --- UPPERCASE HEADING ---` immediately followed by a
  `# ****` divider.
- Prose blocks: explain WHY a method exists, its strengths/costs, what the
  output means (with interpretation ranges), and caveats. Indent continuation
  lines with `#   `.
- `# NOTE: Requires Step X ...` at the top when the script needs objects from
  earlier in-session steps.
- `# TODO:` markers for known improvements.
- End every script with `SUMMARY & PIPELINE MILESTONE TRANSITION`:
  WHERE WE STARTED / WHAT WE HAVE ACCOMPLISHED / WHERE WE ARE HEADING.

**3. R idioms**
- Tidyverse `%>%` chains; Seurat 5 pipeable verbs (`NormalizeData`,
  `FindVariableFeatures`, `ScaleData`, `RunPCA`, `RunUMAP`, `FindNeighbors`,
  `FindClusters`, `IntegrateLayers`, `JoinLayers`, `RunHarmony`,
  `AddMetaData`, `Embeddings`).
- purrr/furrr: `map`, `map_dbl`, `imap_dfr`, `iwalk`, `set_names`, `reduce`,
  `future_map` (with `.options = furrr_options(seed = TRUE)`).
- Base-R-style helpers are preferred and expected: `file.path`, `read.delim`,
  `write.csv`, `dir.create`, plus ggplot2 modern aliases `labs`, `geom_bar`,
  `ggsave`, `scale_fill_brewer`, `geom_text_repel`. Do NOT rewrite these to
  "more standard" equivalents — keep the codebase uniform.

**4. `LOG_STEP`** — wrap anything slow or verbose (defined in part-3.1):
`LOG_STEP("Running Harmony integration...", { ... })`. Prints start, times the
block, prints "✓ Done in Xs", flushes stdout. Use for all multi-second steps.

**5. Parallelism** — `plan(multisession, workers = ...)` then `plan(sequential)`
once finished. Size worker counts from `N_WORKERS`/`availableCores()`; be wary
of serializing multi-GB globals to workers.

**6. Checkpoints** — persist heavy intermediates (merged/integrated objects) to
`DATA_CHECKPOINT_DIR` honoring `CHECKPOINT_FORMAT` (`qs2::qs_save`/`qs_read`,
`nthreads = N_WORKERS`, or `saveRDS`/`readRDS`). Loader pattern lives in
part-5.1 via a named `CHECKPOINT_NAMES` vector.

**7. Metadata** — always read `2_input/sample-metadata/sample_names.tsv` rather
than hardcoding sample lists; that file is authoritative (single source of truth).

**8. Packages** — never add `install.packages()` for a package pixi manages.
Only `part-2.0` / `part-3.0` install non-conda packages (colorout,
SeuratWrappers), each guarded by `if (!requireNamespace(...))`.

**9. Dependency manager (part-3.0-dependencies.R)** — base-R only.
`DEP_TABLE` maps each step → `requires` + `sentinel` objects;
`ensure_dependencies(step = "<filename>")` at the top of each script skips if
sentinel objects exist, otherwise offers to source prerequisites (interactive)
or errors with the missing list (batch). **When adding a new script**: add one
row to `DEP_TABLE` (see "HOW TO ADD A NEW STEP" in that file) AND call
`ensure_dependencies()` at the top — both, not one. `part-6.1` does not yet
have either; it currently relies on `integrated_final` already being in the
session from `part-5.3B`.

**10. Output file naming conventions** — everything nests under
`3_output/<RUN_ID>/`, split by stage:
- Part 2 (`qc_and_filtering/`): `plots/`, `filtered_data/`, `metrics/` — each
  with a per-sample subdirectory `<sample_id>/`.
- Part 3 (`integration_and_clustering/`): `plots/` (subdirs
  `integration_comparison/`, `clustering/`), `integrated_data/`, `metadata/`,
  `checkpoints/`.

File types follow a strict, readable scheme (lowercase snake_case everywhere):
- **Plots** — `NN_<slug>.png`, where `NN` is a zero-padded 2-digit run-order
  (chronological position of the figure in the pipeline, so files sort and
  read in narrative order) and `<slug>` the snake_case content, e.g.
  `01_naive_merge_batch_effects.png`, `02_integration_by_sample.png`,
  `03_integration_by_condition.png`, `04_mixing_scores.png`,
  `05_mixing_vs_preservation.png`; part-2 QC uses e.g.
  `03_qc_violins.png`, `04_qc_scatter.png`, `05_filtering_thresholds_log/linear.png`.
  Two views of the same figure share a number + a `_<variant>` suffix
  (`_log`/`_linear`). Always save via
  `ggsave(file.path(PLOTS_OUT_DIR, "NN_slug.png"), plot = p_<slug>, width, height, dpi = 300)`;
  plot objects are named `p_<slug>` to match.
- **Tables / CSVs** — descriptive snake_case with **no** numeric prefix:
  `integration_comparison_metrics.csv`, `cell_metadata.csv`, per-sample
  `<sample_id>_QC_summary.csv`. Written with `write.csv(..., row.names = FALSE)`
  into `metadata/` (part-3) or `metrics/` (part-2).
- **Seurat objects / checkpoints** — `NN_<slug>.<ext>` with the same 2-digit
  run-order prefix: `01_merged_naive`, `02_integrated_cca`,
  `03_integrated_rpca`, `04_integrated_harmony`, `05_integrated_fastmnn`.
  Extension is `.qs2` or `.rds` depending on `CHECKPOINT_FORMAT`; slugs are
  mapped by the named `CHECKPOINT_NAMES` vector in part-5.1. Part-2 clean
  sample objects are `<sample_id>_qc_filtered.rds` in `filtered_data/<sample_id>/`.
- After every save, log the path with `cat("→ Saved:", <path>, "\n")`.

## Verification
- No formal test framework. Syntax-check changed `.R` files:
  `Rscript -e 'invisible(parse(file = "path/to/script.R"))'`
- Functional verification: run scripts in order inside the pixi part3 env.

## Status & next steps
- Sep 2026 sessions: per-script dependency manager + PREREQUISITES headers
  (`42d4ec8`); helper renamed to `part-3.0-dependencies.R` (`f4f836f`);
  clustering at multiple resolutions added, `part-6.1` (`ec97805`); this
  AGENTS.md added; output file naming conventions (plots/tables/checkpoints)
  documented.
- Pending: wire `part-6.1` into `DEP_TABLE` + add its `ensure_dependencies()`
  block; align its "STEP 15" banner/filename if we touch it later; downstream
  cell-type annotation (Part 4 scaffolded in pixi.toml, commented out).
- After each work session, update this section so the next session resumes.