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
  cluster at multiple resolutions (`part-6.1`) → visualize + evaluate
  resolutions (`part-6.2A`, `part-6.2B`, checkpointed) → load clustering
  checkpoints (`part-6.3`) → cluster quality (`part-6.4`). Scripts:
  `part-3.0-install-packages.R`, `part-3.0-dependencies.R`,
  `part-3.1`…`part-5.3B`, `part-6.1`, `part-6.2A`, `part-6.2B`,
  `part-6.3`, `part-6.4`.
- `2_input/` — `cellranger-matrix-counts/`, and `sample-metadata/sample_names.tsv`
  (the **single source of truth** for sample_id/condition/patient_id/QC thresholds).
- `3_output/<RUN_ID>/` — `qc_and_filtering/` (plots, filtered_data, metrics) and
  `integration_and_clustering/` (plots, integrated_data, metadata, checkpoints).
- `3_guide/` — tutorial clippings from ngs101.com that the scripts' commentary
  is written against: `Part 1 - From FASTQ to Count Matrix.md`,
  `Part 2A/2B - Quality Control and Cell Filtering.md`,
  `Part 3 - Integration and Clustering.md`. Comments cite them as
  "guide §<section>" / "STEP <N>", and each script's SUMMARY & PIPELINE
  MILESTONE TRANSITION mirrors the guide's narrative — keep new steps aligned
  with it.
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
part-5.1 via a named `CHECKPOINT_NAMES` vector, and in part-6.3 for the Part 6
objects. **When adding or modernizing a step, always decide whether it needs a
checkpoint, and state that decision in the script's commentary.** Add one when
the step is expensive (minutes+) or produces a heavy, reusable intermediate
that downstream steps would otherwise recompute; skip it for cheap steps and
say why. Any step that writes a checkpoint must have a matching loader
(mirroring part-5.1/part-6.3), with both registered in `DEP_TABLE` so consumers
depend on the loader rather than the heavy producers.

**7. Metadata** — always read `2_input/sample-metadata/sample_names.tsv` rather
than hardcoding sample lists; that file is authoritative (single source of truth).
It is read once in `part-3.1`, which defines `sample_metadata` plus the shared
`sample_colors` / `condition_colors` palettes every plotting step reuses — do
not re-read the TSV or redefine those colours locally.

**8. Packages** — never add `install.packages()` for a package pixi manages.
Only `part-2.0` / `part-3.0` install non-conda packages (colorout,
SeuratWrappers), each guarded by `if (!requireNamespace(...))`.

**9. Dependency manager (part-3.0-dependencies.R)** — base-R only.
`DEP_TABLE` maps each step → `requires` + `sentinel` objects;
`ensure_dependencies(step = "<filename>")` at the top of each script skips if
sentinel objects exist, otherwise offers to source prerequisites (interactive)
or errors with the missing list (batch). Every script **unconditionally
re-sources** the manager at the top (base-R only, cheap), so edits to
`DEP_TABLE` are picked up on every run without restarting the session or
clearing stale objects. When all prerequisites are already satisfied, a
top-level interactive run is asked whether to re-source them fresh (re-runs
slow steps) or proceed as-is (default); batch and nested auto-sourcing runs
skip silently (as do steps with no prerequisites at all).
**When adding a new script**: add one
row to `DEP_TABLE` (see "HOW TO ADD A NEW STEP" in that file) AND call
`ensure_dependencies()` at the top — both, not one. `part-6.1` follows this
(requires `part-3.1` + `part-5.3B`, sentinels `integrated_final` /
`reduction_final`). Note `FindClusters()` overwrites `seurat_clusters` each
run, so after `part-6.1` downstream code must read the `clusters_res_*`
columns it writes.

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
  `03_integrated_rpca`, `04_integrated_harmony`, `05_integrated_fastmnn`,
  `06_clustered_final` (final clustered object from part-6.3),
  `07_clustering_metrics` (small list: reduction_final, resolutions,
  optimal_resolution, resolution_comparison).
  Extension is `.qs2` or `.rds` depending on `CHECKPOINT_FORMAT`; Part-5
  slugs are mapped by the named `CHECKPOINT_NAMES` vector in part-5.1, and
  the Part-6 pair (06/07) is read by part-6.3. Part-2 clean sample objects
  are `<sample_id>_qc_filtered.rds` in `filtered_data/<sample_id>/`.
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
- `part-6.1` modernized: wired into `DEP_TABLE` + `ensure_dependencies()` block,
  "STEP 15" banner aligned to "STEP 6.1", and `clusters_res_*` columns via
  `AddMetaData` (no colname gsub). Kept lean in scope (no added outputs/plots)
  but fully commented, with an interpretation guide and a SUMMARY & PIPELINE
  MILESTONE TRANSITION matching the 3_guide narrative.
- `part-6.4` modernized (house-style only): wired into `DEP_TABLE` +
  `ensure_dependencies()` block, "STEP 18" banner aligned to "STEP 6.4",
  `Post_Treatment`→`Periodontitis_Post_Treatment` palette fix, and the two
  composition/size summaries now persist as CSVs
  (`cluster_size_summary.csv`, `cluster_sample_composition.csv`,
  `cluster_condition_composition.csv`). Same analyses and plots as before.
- Skill `.opencode/skills/scrna-pipeline-modernize/` encodes the
  modernize-vs-improve distinction and the house-style checklist; load it when
  asked to "modernize"/align a pipeline script to conventions without changing
  its analysis.
- Shared cohort metadata + colour palettes centralized in `part-3.1`
  (`sample_metadata`, `sample_colors`, `condition_colors`); `3.2`, `3.3B`,
  `5.2A`, and `6.4` now consume them instead of re-reading the TSV or
  redefining palettes locally.
- Part 6 checkpointing: `part-6.2B` now writes `06_clustered_final` (the final
  clustered object) + `07_clustering_metrics` (`reduction_final`,
  `resolutions`, `optimal_resolution`, `resolution_comparison`) to
  `DATA_CHECKPOINT_DIR`; new `part-6.3-load-clustering-checkpoints.R` loads
  them back (mirrors `part-5.1`). `part-6.4` now requires the `6.3` loader
  instead of `6.1`/`6.2B`, so downstream sessions skip the heavy
  `5.3B`−`6.2B` chain once the checkpoints exist. `part-6.2A` still requires
  `6.1` (it is upstream of the checkpoints). First build still runs the
  full chain; `6.2B` itself does not self-short-circuit.
- `part-6.5` and `part-6.6` modernized (house-style + consistency only): wired
  into `DEP_TABLE` + `ensure_dependencies()`, banners `6.5`/`6.6` (was STEP
  19/20), PREREQUISITES + WHY/section prose + SUMMARY transitions, saves logged
  with `cat("→ Saved:")`, and `6.5`'s stale `best_method` TODO removed. Fixed
  the underlying dependency bug: `part-5.3B`'s `ready()` now also requires
  `best_method`, so `6.5`/`6.6` (which both use it) get 5.3B properly
  re-sourced rather than silently skipped. `6.6` now honours
  `CHECKPOINT_FORMAT` for the final object extension (`integrated_clustered_seurat.qs2`/`.rds`)
  and wraps the save in `LOG_STEP`; no new analysis, plots, or filenames.
- Pending: downstream cell-type annotation (Part 4 scaffolded in pixi.toml,
  commented out).
- After each work session, update this section so the next session resumes.