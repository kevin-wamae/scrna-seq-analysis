---
name: scrna-pipeline-modernize
description: Use when asked to "modernize", "align to conventions", "bring a script up to house/house-style", "make it match the other scripts", "add commenting and transition notes", or otherwise update an R script (usually a part-*.R) in the scrna-seq-analysis pipeline WITHOUT changing its analysis. Also load when working on any pipeline script that lacks the dependency-manager block, banner, commentary, output conventions, or SUMMARY & PIPELINE MILESTONE TRANSITION. Distinguishes "modernize" (comment/wire/docs only — never change analysis scope) from "improve" (scope may expand; confirm first). Use ONLY for scripts in this repository's pipeline folders (1_scripts/part-2_*, part-3_*).
---

# Modernizing scrna-seq-analysis pipeline scripts

The repo's scripts are written to be *read as much as run*: every script
explains the biology/statistics, frames outputs as ranges with interpretation,
and ends with a milestone-transition summary. "Modernize" means bringing a
script up to that house style **without touching what it computes or plots**.

## The core distinction

- **Modernize / align to conventions** → add commentary, wiring, and docs
  only. The analysis, thresholds, outputs, and plots stay exactly as they are.
- **Improve** → changes the analysis (new metrics, new plots, new outputs).
  Never guess this is wanted. If the user says "improve", ask whether they
  mean house-style-only or actual analysis changes.

Read `AGENTS.md` (`§Coding conventions`) first — it is the source of truth
and this skill summarises it for the modernization workflow.

## Modernize checklist

Apply each item that the target script is missing. Do not invent work beyond
the list.

1. **Banner** — `# ***` line, `# STEP <N>: Title`, closing `# ***`. The
   STEP number follows the *running step* (e.g. `STEP 6.4`), not the guide's
   historical numbering. Keep the existing title unless it is wrong.

2. **PREREQUISITES block + dependency-manager wiring** — at the top:
   - `# --- PREREQUISITES (scripts that must run before this one in this
     session) ---` with a `# ****` divider, bullet-listing each prerequisite,
   - the standard "dependency manager is always re-sourced here…" paragraph
     (copy verbatim from a neighbouring script, e.g. `part-6.4`),
   - `source("1_scripts/part-3_integration-and-clustering/part-3.0-dependencies.R")`
   - `ensure_dependencies(step = "<this script's exact filename>")`
   - `# NOTE: Requires Step X …` when it needs objects from a specific
     earlier in-session step.

3. **DEP_TABLE row** — add one entry for the script in
   `part-3.0-dependencies.R`:
   - `requires`: the scripts listed in the PREREQUISITES block,
   - `sentinel`: an object the script actually creates (check it via
     `exists()`), which downstream steps use to detect that this step ran.
   The row key MUST exactly match the filename passed to
   `ensure_dependencies()` and the script's filename (including `.R`).

4. **Section prose** — add `# --- UPPERCASE HEADING ---` headers with `# ****`
   dividers, and prose blocks explaining WHY each check exists, what its
   output means with interpretation ranges, and caveats. Indent continuation
   lines with `#   `. Reference the tutorial guide where relevant (e.g.
   `guide §7.7`).

5. **Output conventions** — align file names to `AGENTS.md` §10 while keeping
   the *same* outputs:
   - Plots: `NN_<slug>.png` in the stage's `plots/` subdir, saved via
     `ggsave(..., width, height, dpi = 300)`, plot objects named `p_<slug>`.
     Keep the existing NN run-order numbers.
   - Tables/CSVs: descriptive snake_case, `write.csv(..., row.names = FALSE)`,
     into `metadata/` (part-3) or `metrics/` (part-2).
   - After every save: `cat("→ Saved:", <path>, "\n")`.
   - If the original script only `print()`s a summary that the conventions
     expect as a table, persisting it as a CSV is a *delivery* change, not an
     analysis change — but confirm with the user before adding it.

6. **SUMMARY & PIPELINE MILESTONE TRANSITION** — end the script with the
   standard block: WHERE WE STARTED / WHAT WE HAVE ACCOMPLISHED / WHERE WE
   ARE HEADING, mirroring the target script's neighbour (e.g. `part-6.3`) and
   the guide narrative. Start the `# ***` footer line like the header.

## Do-not rules (violations break repo conventions)

- **No scope creep.** Never add analytically new metrics, thresholds, plots,
  or computations under the guise of "modernizing".
- No `install.packages()` for packages pixi manages (see AGENTS.md §8).
- Don't rewrite base-R-style helpers (`file.path`, `read.delim`, `write.csv`,
  `dir.create`, ggplot2 aliases) into "more standard" equivalents — keep the
  codebase uniform.
- Keep existing package choices (e.g. `reshape2::melt` is loaded in
  `part-3.1` and used deliberately; don't swap it for tidyr).
- Match the exact wording of the shared "dependency manager is always
  re-sourced here…" paragraph and colours (e.g. condition palette in
  `part-3.3B`) so scripts stay consistent.

## Verification

- Syntax-check the changed `.R` file(s):
  `Rscript -e 'invisible(parse(file = "path/to/script.R"))'`
- Confirm the DEP_TABLE entry resolves: source
  `part-3.0-dependencies.R` and run `dep_resolve("<script filename>")`.
- If outputs were the point, the CSVs land under
  `3_output/<RUN_ID>/integration_and_clustering/metadata/` and plots under
  `.../plots/clustering/`.