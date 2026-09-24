# ****************************************************************************#
# STEP 3.0: Dependency manager for the integration & clustering pipeline
# ****************************************************************************#
#   Every script in this folder passes in-memory objects to the next one
#   (e.g. `seurat_list`, `merged_seurat`, `mixing_results`), so a step run
#   on its own in a fresh session errors the moment it touches an object an
#   earlier script should have created. This manager makes that safe:
#
#   - A per-script dependency table maps each step to the scripts that must
#     run before it, and the objects those scripts are expected to leave in
#     the session (its "sentinel" objects).
#   - `ensure_dependencies()` checks whether those sentinels already exist.
#       * Some missing -> prints the exact list that needs to run, then
#         (interactively) asks whether to source them now, or (in a
#         non-interactive/batch run) stops with a clear error listing them.
#       * All present  -> in a non-interactive or nested (auto-sourcing) run,
#         proceeds with the current session's objects silently. In a
#         top-level interactive session it instead asks whether to re-source
#         the prerequisites fresh (default: no — proceed as-is). Re-sourcing
#         guarantees a fresh state but re-runs any slow steps.
#
#   BASE-R ONLY: no packages are loaded here, so this can run before Step
#   3.1's `library()` calls. Sourcing does not change the working directory,
#   so the existing relative paths ("2_input/...", "3_output/...") still
#   resolve from the project root, exactly as they do when the scripts are
#   run by hand.


# --- LOCATE THIS FOLDER ---
# ****************************************************************************#
#   The folder containing the pipeline scripts. When this file is sourced at
#   top level, the caller's frame carries the sourced file's path in `ofile`,
#   so sibling scripts can be resolved without hardcoding a path. That frame
#   trick is unreliable in a NESTED source — e.g. when `ensure_dependencies()`
#   sources a prerequisite script whose own top-of-file block re-sources this
#   manager from inside a function, `sys.frame(1)` no longer carries `ofile`
#   and the naive fallback wrongly resolves to the working directory (the repo
#   root) — so every candidate below is validated against the files it must
#   actually contain before being accepted.
DEP_DIR <- local({
    dep_dir_ok <- function(d) {
        !is.null(d) && length(d) == 1L && !is.na(d) && nzchar(d) &&
            file.exists(file.path(d, "part-3.0-dependencies.R")) &&
            file.exists(file.path(d, "part-3.1-load-libraries-and-configuration.R"))
    }
    frame_ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
    cand <- if (
        is.character(frame_ofile) && length(frame_ofile) == 1L && nzchar(frame_ofile)
    ) {
        normalizePath(dirname(frame_ofile), mustWork = FALSE)
    } else {
        NA_character_
    }
    if (dep_dir_ok(cand)) return(cand)
    # Fall back to the documented convention: the working directory is the
    # repo root, so the scripts live at the fixed relative location below.
    fallback <- normalizePath(
        file.path(getwd(), "1_scripts", "part-3_integration-and-clustering"),
        mustWork = FALSE
    )
    if (dep_dir_ok(fallback)) return(fallback)
    # The caller may already be running from inside the scripts folder.
    wd <- normalizePath(getwd(), mustWork = FALSE)
    if (dep_dir_ok(wd)) return(wd)
    stop(
        "part-3.0-dependencies.R: could not locate the pipeline script folder.\n",
        "Run with the repository root (or the part-3 folder) as the working directory."
    )
})


# --- PER-SCRIPT DEPENDENCY TABLE ---
# ****************************************************************************#
#   Each entry maps a step (by filename) to:
#     - `requires` : prerequisite scripts, in run order (topological order is
#                    resolved automatically from this list).
#     - `sentinel` : object names to check with `exists()` to decide whether
#                    the step has already been run in this session.
#     - `ready`    : optional function() -> TRUE/FALSE for steps whose
#                    "already done" state can't be read from a single object
#                    (e.g. packages installed, or a column appended to a
#                    data frame). If provided, the step counts as satisfied
#                    when `ready()` is TRUE or all `sentinel` objects exist.
#
# HOW TO ADD A NEW STEP (e.g. part-10.1):
#   Add ONE row here, listing only that step's direct prerequisites and the
#   objects it needs:
#
#       "part-10.1-<slug>.R" = list(
#           requires = c("part-5.2A-<slug>.R", "part-9.5B-<slug>.R"),
#           sentinel = c("<object1>", "<object2>")   # what 10.1 needs
#       )
#
#   `ensure_dependencies()` computes ONLY that step's transitive closure in
#   topological order and checks only those sentinels — unrelated steps
#   (3.3B, 4.2C, ...) are never sourced or touched. No other script needs
#   editing; each new script just calls `ensure_dependencies(step = "<its
#   own filename>")` in its own top block.
DEP_TABLE <- list(
    "part-3.0-install-packages.R" = list(
        requires = character(),
        sentinel = character(),
        ready = function()
            requireNamespace("colorout", quietly = TRUE) &&
            requireNamespace("SeuratWrappers", quietly = TRUE)
    ),
    "part-3.1-load-libraries-and-configuration.R" = list(
        requires = c("part-3.0-install-packages.R"),
        sentinel = c("LOG_STEP", "DATA_CHECKPOINT_DIR")
    ),
    "part-3.2-load-10x-quality-controlled-data.R" = list(
        requires = c("part-3.1-load-libraries-and-configuration.R"),
        sentinel = c("seurat_list")
    ),
    "part-3.3A-merge-naive.R" = list(
        requires = c(
            "part-3.1-load-libraries-and-configuration.R",
            "part-3.2-load-10x-quality-controlled-data.R"
        ),
        sentinel = c("merged_naive")
    ),
    "part-3.3B-visualize-batch-effects.R" = list(
        requires = c(
            "part-3.1-load-libraries-and-configuration.R",
            "part-3.2-load-10x-quality-controlled-data.R",
            "part-3.3A-merge-naive.R"
        ),
        sentinel = c("merged_naive", "sample_metadata")
    ),
    "part-3.3C-quantifying-batch-effects.R" = list(
        requires = c(
            "part-3.1-load-libraries-and-configuration.R",
            "part-3.2-load-10x-quality-controlled-data.R",
            "part-3.3A-merge-naive.R"
        ),
        sentinel = c("cluster_composition")
    ),
    "part-4.1-prepare-data-for-integration.R" = list(
        requires = c(
            "part-3.1-load-libraries-and-configuration.R",
            "part-3.2-load-10x-quality-controlled-data.R"
        ),
        sentinel = c("merged_seurat")
    ),
    "part-4.2A-integration-method-1-CCA.R" = list(
        requires = c(
            "part-3.1-load-libraries-and-configuration.R",
            "part-4.1-prepare-data-for-integration.R"
        ),
        sentinel = c("integrated_cca")
    ),
    "part-4.2B-integration-method-2-RCPA.R" = list(
        requires = c(
            "part-3.1-load-libraries-and-configuration.R",
            "part-4.1-prepare-data-for-integration.R"
        ),
        sentinel = c("integrated_rpca")
    ),
    "part-4.2C-integration-method-3-HARMONY.R" = list(
        requires = c(
            "part-3.1-load-libraries-and-configuration.R",
            "part-4.1-prepare-data-for-integration.R"
        ),
        sentinel = c("integrated_harmony")
    ),
    "part-4.2D-integration-method-3-FastMNN.R" = list(
        requires = c(
            "part-3.1-load-libraries-and-configuration.R",
            "part-4.1-prepare-data-for-integration.R"
        ),
        sentinel = c("integrated_fastmnn")
    ),
    "part-5.1-load-seurat-object-checkpoints.R" = list(
        requires = c("part-3.1-load-libraries-and-configuration.R"),
        sentinel = c(
            "merged_naive", "integrated_cca", "integrated_rpca",
            "integrated_harmony", "integrated_fastmnn"
        )
    ),
    "part-5.2A-integration-comparison-visual.R" = list(
        requires = c(
            "part-3.1-load-libraries-and-configuration.R",
            "part-5.1-load-seurat-object-checkpoints.R"
        ),
        sentinel = c(
            "merged_naive", "integrated_cca", "integrated_rpca",
            "integrated_harmony", "integrated_fastmnn"
        )
    ),
    "part-5.2B-integration-comparison-quantitatively.R" = list(
        requires = c(
            "part-3.1-load-libraries-and-configuration.R",
            "part-5.1-load-seurat-object-checkpoints.R"
        ),
        sentinel = c("methods_list", "mixing_results")
    ),
    "part-5.3A-assess-biological-preservation.R" = list(
        requires = c(
            "part-3.1-load-libraries-and-configuration.R",
            "part-5.1-load-seurat-object-checkpoints.R",
            "part-5.2B-integration-comparison-quantitatively.R"
        ),
        sentinel = c("separation_results")
    ),
    "part-5.3B-select-best-integration-method.R" = list(
        requires = c(
            "part-3.1-load-libraries-and-configuration.R",
            "part-5.1-load-seurat-object-checkpoints.R",
            "part-5.2B-integration-comparison-quantitatively.R",
            "part-5.3A-assess-biological-preservation.R"
        ),
        sentinel = character(),
        ready = function()
            exists("methods_list",  envir = .GlobalEnv, inherits = FALSE) &&
            exists("mixing_results", envir = .GlobalEnv, inherits = FALSE) &&
            "condition_separation" %in% names(mixing_results)
    ),
    "part-6.1-clustering-multiple-resolutions.R" = list(
        requires = c(
            "part-3.1-load-libraries-and-configuration.R",
            "part-5.3B-select-best-integration-method.R"
        ),
        sentinel = c("integrated_final", "reduction_final")
    ),
    "part-6.2-view-multiple-resolutions-clustering.R" = list(
        requires = c(
            "part-3.1-load-libraries-and-configuration.R",
            "part-6.1-clustering-multiple-resolutions.R"
        ),
        sentinel = c("integrated_final", "resolutions")
    ),
    "part-6.3-evaluate-optimal-clustering-resolution.R" = list(
        requires = c(
            "part-3.1-load-libraries-and-configuration.R",
            "part-5.3B-select-best-integration-method.R",
            "part-6.1-clustering-multiple-resolutions.R"
        ),
        sentinel = c("resolution_comparison", "optimal_resolution")
    ),
    "part-6.4-assess-cluster-quality-and-stability.R" = list(
        requires = c(
            "part-3.1-load-libraries-and-configuration.R",
            "part-6.3-evaluate-optimal-clustering-resolution.R"
        ),
        sentinel = c("cluster_sizes")
    )
)


# --- INTERNALS ---
# ****************************************************************************#

# Is one step already satisfied in the current session?
dep_satisfied <- function(entry) {
    if (!is.null(entry$ready)) {
        ok <- tryCatch(isTRUE(entry$ready()), error = function(e) FALSE)
        if (ok) return(TRUE)
        if (length(entry$sentinel) == 0L) return(FALSE)
    }
    all(vapply(
        entry$sentinel,
        function(nm) exists(nm, envir = .GlobalEnv, inherits = FALSE),
        logical(1)
    ))
}

# Return the transitive closure of a step's prerequisites in topological
# order (every prerequisite comes before the steps that depend on it).
dep_resolve <- function(step) {
    if (!step %in% names(DEP_TABLE)) {
        stop("ensure_dependencies: unknown step '", step,
            "'. Please add it to DEP_TABLE in part-3.0-dependencies.R.")
    }
    seen <- character()
    order <- character()
    visit <- function(s) {
        if (s %in% seen) return(invisible())
        seen <<- c(seen, s)
        for (r in DEP_TABLE[[s]]$requires) visit(r)
        order <<- c(order, s)
    }
    visit(step)
    order
}


# --- PUBLIC API ---
# ****************************************************************************#

# Ensure the current script's prerequisites have been run in this session.
# Call this at the top of each pipeline script, e.g.:
#   ensure_dependencies(step = "part-5.3B-select-best-integration-method.R")
#
# Returns invisible(TRUE) if prerequisites were (re)sourced, invisible(FALSE)
# if they were already present and kept as-is (or the user declined to act).
ensure_dependencies <- function(step = NULL) {
    if (is.null(step) || !nzchar(step)) {
        stop("ensure_dependencies: 'step' must be the current script's filename.")
    }

    order  <- dep_resolve(step)
    to_check <- setdiff(order, step)          # the current step itself is last
    missing <- vapply(
        to_check,
        function(s) !dep_satisfied(DEP_TABLE[[s]]),
        logical(1)
    )
    need <- to_check[missing]

    if (length(need) == 0L) {
        # --- ALL PREREQUISITES ALREADY SATISFIED ---
        # Silent fast path for batch, nested auto-sourcing runs, and steps
        # that have no prerequisites at all (e.g. part-3.0-install-packages,
        # where there is nothing to re-source). A top-level interactive run
        # with an actual prerequisite chain offers the choice between a
        # guaranteed-fresh state (re-source the whole chain, re-running slow
        # steps) and keeping the current session's objects (default).
        if (!interactive() || isTRUE(getOption("dep_auto_source")) ||
                length(to_check) == 0L) {
            cat(sprintf(
                "\n[DEP] %s: all prerequisites already satisfied — proceeding with current session objects.\n",
                step
            ))
            return(invisible(FALSE))
        }
        ans <- readline(
            "\n[DEP] All prerequisites already satisfied.\n     Re-source them fresh? (re-runs any slow steps) [y/N]: "
        )
        if (!(tolower(trimws(ans)) %in% c("y", "yes"))) {
            cat(sprintf(
                "\n[DEP] %s: proceeding with current session objects (declined fresh re-source).\n",
                step
            ))
            return(invisible(FALSE))
        }
        to_source <- to_check
        cat("[DEP] Re-sourcing satisfied prerequisites fresh...\n")

    } else {
        # --- SOME PREREQUISITES MISSING ---
        need_names <- basename(need)
        cat(sprintf(
            "\n[DEP] %s requires these prerequisite script(s) not yet run in this session:\n",
            step
        ))
        for (n in need_names) cat("   - ", n, "\n", sep = "")

        run <- if (isTRUE(getOption("dep_auto_source"))) {
            TRUE                                     # already mid auto-resolution
        } else if (interactive()) {
            ans <- readline("  Run them now? [y/N]: ")
            tolower(trimws(ans)) %in% c("y", "yes")
        } else {
            FALSE
        }

        if (run) {
            to_source <- need
            cat("[DEP] Sourcing prerequisite scripts in order...\n")
        } else if (!interactive()) {
            stop(sprintf(
                "Missing prerequisite(s) for %s in a non-interactive session: %s.\nRun them first — see the PREREQUISITES note at the top of this script.",
                step, paste(need_names, collapse = ", ")
            ))
        } else {
            cat(sprintf(
                "\n[DEP] Declined to run prerequisites for %s — continuing (may error if the required objects are absent).\n",
                step
            ))
            return(invisible(FALSE))
        }
    }

    # Shared source loop (used by both the missing-prereqs and fresh re-source
    # paths). `dep_auto_source` suppresses nested prompts while sourcing.
    old <- getOption("dep_auto_source")
    options(dep_auto_source = TRUE)          # suppress nested prompts
    on.exit(options(dep_auto_source = old), add = TRUE)
    for (s in to_source) {
        p <- file.path(DEP_DIR, s)
        if (!file.exists(p)) {
            stop("Prerequisite script not found: ", p)
        }
        cat("  → ", basename(s), "\n", sep = "")
        source(p, local = FALSE, echo = FALSE)
    }
    invisible(TRUE)
}
