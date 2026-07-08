# ==============================================================================
# zzz.R
# Package loader / dependency declarations
# ==============================================================================

# Required packages
.required_packages <- c(
    "dplyr",
    "tibble",
    "tidyr",
    "ggplot2",
    "scales",
    "stringr",
    "glue",
    "BiocParallel",
    "lme4",
    "glmmTMB",
    "broom.mixed",
    "Rcpp",
    "RcppEigen"
)

# Check and load dependencies
.check_deps <- function() {
    missing <- .required_packages[!vapply(.required_packages, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing)) {
        stop("Missing required packages: ", paste(missing, collapse = ", "),
             "\nInstall with: install.packages(c('", paste(missing, collapse = "', '"), "'))",
             call. = FALSE)
    }
    invisible(TRUE)
}

# Source all R files in order
.load_all <- function(path = NULL) {
    .check_deps()

    # If no path supplied, derive from sourced file location or fall back to "."
    if (is.null(path)) {
        this_file <- sys.frame(1)$ofile
        path <- if (!is.null(this_file)) dirname(this_file) else "."
    }

    # Attach commonly used packages
    suppressPackageStartupMessages({
        library(dplyr)
        library(tibble)
    })

    # Core numbered library files live in core/ (fall back to R/ then path/ for older layouts)
    core_dir <- file.path(path, "core")
    if (!dir.exists(core_dir)) core_dir <- file.path(path, "R")
    if (!dir.exists(core_dir)) core_dir <- path
    files <- list.files(core_dir, pattern = "^[0-9]+-.*\\.R$", full.names = TRUE)
    files <- files[order(files)]  # ensure numeric ordering
    
    for (f in files) {
        source(f, local = FALSE)
        message("Loaded: ", basename(f))
    }

    # Source helper files if helpers/ directory exists
    helpers_dir <- file.path(path, "helpers")
    if (dir.exists(helpers_dir)) {
        helper_files <- list.files(helpers_dir, pattern = "\\.R$", full.names = TRUE)
        for (f in helper_files) {
            source(f, local = FALSE)
            message("Loaded helper: ", basename(f))
        }
        ## Compile any C++ helpers via Rcpp::sourceCpp. Compiled binaries
        ## are cached under tools::R_user_dir or the system cache so this is
        ## fast on second load.
        cpp_files <- list.files(helpers_dir, pattern = "\\.cpp$", full.names = TRUE)
        ## Skip TMB templates: detected by the literal `#include <TMB.hpp>` token.
        ## TMB templates are compiled lazily by their R wrappers via TMB::compile()
        ## (which knows about the TMB include path).
        if (length(cpp_files)) {
            keep <- vapply(cpp_files, function(f) {
                hd <- tryCatch(readLines(f, n = 50, warn = FALSE),
                                error = function(e) character(0))
                !any(grepl("TMB\\.hpp", hd))
            }, logical(1))
            cpp_files <- cpp_files[keep]
        }
        if (length(cpp_files)) {
            if (!requireNamespace("Rcpp", quietly = TRUE)) {
                warning("Rcpp not installed; skipping compilation of: ",
                        paste(basename(cpp_files), collapse = ", "))
            } else {
                for (f in cpp_files) {
                    tryCatch({
                        Rcpp::sourceCpp(f, verbose = FALSE)
                        message("Compiled helper: ", basename(f))
                    }, error = function(e) {
                        warning(sprintf("Failed to compile %s: %s",
                                        basename(f), conditionMessage(e)))
                    })
                }
            }
        }
    }

    invisible(TRUE)
}

# Auto-load if this file is sourced directly
if (sys.nframe() == 0 || identical(environment(), globalenv())) {
    ofile <- sys.frame(1)$ofile
    .load_all(if (!is.null(ofile)) dirname(ofile) else ".")
}
