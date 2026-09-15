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
#       * All present  -> prints a skip message and does nothing.
#       * Some missing -> prints the exact list that needs to run, then
#         (interactively) asks whether to source them now, or (in a
#         non-interactive/batch run) stops with a clear error listing them.
#
#   BASE-R ONLY: no packages are loaded here, so this can run before Step
#   3.1's `library()` calls. Sourcing does not change the working directory,
#   so the existing relative paths ("2_input/...", "3_output/...") still
#   resolve from the project root, exactly as they do when the scripts are
#   run by hand.


# --- LOCATE THIS FOLDER ---
# ****************************************************************************#
#   When sourced via `source()`, the caller's frame carries the sourced
#   file's path in `ofile`. Use that to resolve sibling scripts without
#   hardcoding a path. Falls back to the working directory if that frame is
#   unavailable (e.g. the helper is run directly).
DEP_FRAME <- tryCatch(sys.frame(1), error = function(e) NULL)
DEP_DIR <- if (
    !is.null(DEP_FRAME) &&
    exists("ofile", envir = DEP_FRAME, inherits = FALSE) &&
    is.character(DEP_FRAME$ofile)
) {
    normalizePath(dirname(DEP_FRAME$ofile))
} else {
    normalizePath(getwd())
}


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
# Returns invisible(TRUE) if prerequisites were sourced, invisible(FALSE) if
# they were already present (or the user declined to run them).
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
        cat(sprintf(
            "\n[DEP] %s: all prerequisites already available in this session — nothing to run.\n",
            step
        ))
        return(invisible(FALSE))
    }

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
        cat("[DEP] Sourcing prerequisite scripts in order...\n")
        old <- getOption("dep_auto_source")
        options(dep_auto_source = TRUE)          # suppress nested prompts
        on.exit(options(dep_auto_source = old), add = TRUE)
        for (s in need) {
            p <- file.path(DEP_DIR, s)
            if (!file.exists(p)) {
                stop("Prerequisite script not found: ", p)
            }
            cat("  → ", basename(s), "\n", sep = "")
            source(p, local = FALSE, echo = FALSE)
        }
        return(invisible(TRUE))
    }

    if (!interactive()) {
        stop(sprintf(
            "Missing prerequisite(s) for %s in a non-interactive session: %s.\nRun them first — see the PREREQUISITES note at the top of this script.",
            step, paste(need_names, collapse = ", ")
        ))
    }

    cat(sprintf(
        "\n[DEP] Declined to run prerequisites for %s — continuing (may error if the required objects are absent).\n",
        step
    ))
    invisible(FALSE)
}
