# ==============================================================================
# Near/Far Spillover Decomposition for PACE
# ==============================================================================
#
# Key insight from the literature:
#   - Transcript misassignment operates at short range (~10-15µm centroid-centroid)
#   - Biological cell-cell influence operates at longer range (up to 25+µm)
#   - Marco Salas et al. (Nature Methods, 2025): >10.71µm from centroid,
#     transcripts correlate more with background than cell-type signature
#
# Approach:
#   - Random slopes on TOTAL neighbourhood counts (25µm) → biology
#   - Fixed effects on NEAR neighbourhood counts (15µm) → absorb contamination
#   - Near fixed effects are gene-agnostic (same for all cell types)
#   - Random slopes capture cell-type-specific spatial variation
#
# Model per gene j:
#   gene_j ~ 1 + near_type1 + ... + near_typeK + (1 + type1 + ... + typeK || celltype)
#
# The near FEs absorb short-range contamination. The random slopes on total
# counts capture biology net of contamination.
# ==============================================================================

library(spatstat.geom)

#' Add near-range neighbourhood counts to data frame
#'
#' Computes the count of each cell type within a near radius for each cell.
#' These are used as fixed-effect spillover covariates.
#'
#' @param df Data frame with cell-level data (must have 'celltype' column)
#' @param cell_data Data frame with cell_id, cell_type, x, y
#' @param types Character vector of cell type names (matching column names in df)
#' @param near_radius Near radius in microns (default: 15)
#' @param suffix Suffix for near columns (default: "_near")
#' @return Modified df with added near-range columns
add_near_counts <- function(df, cell_data, types, near_radius = 15, suffix = "_near",
                            edge = "isotropic",
                            kernel = c("none", "exp", "gaussian"),
                            kernel_h = 5, kernel_max_dist = 20) {

    kernel <- match.arg(kernel)

    # Ensure imageID exists (single sample case)
    if (!"imageID" %in% names(cell_data)) {
        cell_data$imageID <- "sample1"
    }

    # Build the near matrix: either binary count (kernel="none") or
    # kernel-weighted (kernel="exp" / "gaussian"). Kernel mode lets the
    # spillover term absorb only very-short-range lateral transcript bleed
    # (Salas 2025) without competing with the broader-range neighbour
    # kernel features for biological signal.
    near_mat <- if (kernel == "none") {
        spillover_by_image(cell_data, radius = near_radius, edge = edge)
    } else {
        if ("imageID" %in% names(cell_data) &&
            length(unique(cell_data$imageID)) > 1L) {
            kernel_spillover_by_image(cell_data, max_dist = kernel_max_dist,
                                       h = kernel_h, kernel = kernel,
                                       edge = edge)
        } else {
            ## Single-image case: use compute_kernel_abundance directly
            mat <- compute_kernel_abundance(cell_data,
                                             max_dist = kernel_max_dist,
                                             h = kernel_h,
                                             kernel = kernel, edge = edge)
            rownames(mat) <- cell_data$cell_id
            mat
        }
    }
    
    # Zero out homotypic counts (a cell's own type is not a contamination source)
    ct <- cell_data$cell_type
    for (s in types) {
        if (s %in% colnames(near_mat)) {
            near_mat[ct == s, s] <- 0
        }
    }
    
    # Match rows: near_mat rownames are cell_id, df rows align with cell_data
    for (s in types) {
        safe_s <- gsub(" ", "_", s)
        col_name <- paste0(safe_s, suffix)
        
        if (s %in% colnames(near_mat)) {
            df[[col_name]] <- as.numeric(near_mat[, s])
        } else {
            df[[col_name]] <- 0
        }
    }
    
    df
}


#' Run near/far spillover decomposition model
#'
#' Wrapper around extractRandomEffectsNew that:
#' 1. Adds near-range counts to df
#' 2. Runs the model with near counts as fixed effects (vars_spill)
#'    and total counts as random slopes (vars)
#'
#' @param genes Gene names
#' @param df Data frame with gene columns and total neighbourhood counts
#' @param vars Cell type covariate names (for 25µm random slopes)
#' @param cell_data Data frame with cell_id, cell_type, x, y
#' @param types Cell type names
#' @param near_radius Near radius for contamination (default: 15)
#' @param family glmmTMB family
#' @param BPPARAM BiocParallel backend
#' @return Results from extractRandomEffectsNew
run_near_far_model <- function(
        genes, df, vars, cell_data, types,
        near_radius = 15,
        family = glmmTMB::nbinom1(),
        include_spillover = TRUE,
        BPPARAM = BiocParallel::MulticoreParam(workers = 4),
        use_offset = TRUE,
        area_var = NULL,
        ...
) {

    # Step 0: Compute per-cell library size if needed and missing
    if (use_offset && !"nCount" %in% names(df)) {
        df <- add_nCount(df, genes)
    }

    # Step 1: Add near-range counts
    cat("Computing near-range neighbourhood counts (radius =", near_radius, "µm)...\n")
    df <- add_near_counts(df, cell_data, types, near_radius = near_radius)

    # Step 2: Define near column names as spillover covariates
    vars_spill <- paste0(gsub(" ", "_", types), "_near")

    cat("Model structure:\n")
    cat("  Library-size offset:           ", if (use_offset) "offset(log(nCount))" else "(none)", "\n")
    cat("  Cell-area covariate:           ", if (!is.null(area_var)) sprintf("scale(log(%s))", area_var) else "(none)", "\n")
    cat("  Fixed effects (contamination): ", paste(vars_spill, collapse = ", "), "\n")
    cat("  Random slopes (biology):       ", paste(vars, collapse = ", "), "\n")

    # Step 3: Run extraction
    results <- extractRandomEffectsNew(
        genes = genes,
        df = df,
        vars = vars,
        vars_spill = vars_spill,
        include_spillover = include_spillover,
        family = family,
        BPPARAM = BPPARAM,
        use_offset = use_offset,
        area_var = area_var,
        ...
    )

    results
}


#' Run near/far spillover decomposition model (Responder-aware)
#'
#' Wrapper around extractRandomEffectsResponder that:
#' 1. Adds near-range counts to df
#' 2. Runs the responder model with near counts as fixed effects (vars_spill)
#'    and total counts as random slopes with Responder interactions
#'
#' Model per gene j:
#'   gene_j ~ 1 + near_type1 + ... + near_typeK +
#'            (1 + Responder * (type1 + ... + typeK) || celltype) +
#'            (1 || imageID)
#'
#' @param genes Gene names
#' @param df Data frame with gene columns, total neighbourhood counts,
#'   Responder column, and imageID column
#' @param vars Cell type covariate names (for 25µm random slopes)
#' @param cell_data Data frame with cell_id, cell_type, x, y
#' @param types Cell type names
#' @param near_radius Near radius for contamination (default: 15)
#' @param resp_var Name of responder column (default: "Responder")
#' @param resp_level Active responder level (default: "PD")
#' @param family glmmTMB family
#' @param BPPARAM BiocParallel backend
#' @return Results from extractRandomEffectsResponder
run_near_far_model_responder <- function(
        genes, df, vars, cell_data, types,
        near_radius = 15,
        resp_var = "Responder",
        resp_level = "PD",
        family = glmmTMB::nbinom1(),
        include_spillover = TRUE,
        BPPARAM = BiocParallel::MulticoreParam(workers = 4),
        use_offset = TRUE,
        area_var = NULL,
        image_slope_vars = NULL,
        ...
) {

    # Step 0: Compute per-cell library size if needed and missing
    if (use_offset && !"nCount" %in% names(df)) {
        df <- add_nCount(df, genes)
    }

    # Step 1: Add near-range counts
    cat("Computing near-range neighbourhood counts (radius =", near_radius, "µm)...\n")
    df <- add_near_counts(df, cell_data, types, near_radius = near_radius)

    # Step 2: Define near column names as spillover covariates
    vars_spill <- paste0(gsub(" ", "_", types), "_near")

    # Resolve TRUE -> vars / FALSE -> NULL up front so the displayed formula
    # matches what the inner fitter actually uses
    if (isTRUE(image_slope_vars)) {
        image_slope_vars <- vars
    } else if (isFALSE(image_slope_vars)) {
        image_slope_vars <- NULL
    }

    img_re_str <- if (!is.null(image_slope_vars) && length(image_slope_vars)) {
        paste0("(1 + ", paste(image_slope_vars, collapse = " + "), " || imageID)")
    } else {
        "(1 || imageID)"
    }

    cat("Model structure (responder-aware):\n")
    cat("  Library-size offset:           ", if (use_offset) "offset(log(nCount))" else "(none)", "\n")
    cat("  Cell-area covariate:           ", if (!is.null(area_var)) sprintf("scale(log(%s))", area_var) else "(none)", "\n")
    cat("  Fixed effects (contamination): ", paste(vars_spill, collapse = ", "), "\n")
    cat("  Fixed Responder (stabiliser):  ", resp_var, "\n")
    cat("  Random slopes (biology):       ", paste(vars, collapse = ", "), "\n")
    cat("  Responder variable:            ", resp_var, " (active level: ", resp_level, ")\n")
    cat("  Image random effects:          ", img_re_str, "\n")
    cat("  Formula per gene:\n")
    offset_str <- if (use_offset) "offset(log(nCount)) + " else ""
    cat("    gene ~ 1 +", offset_str, resp_var, "+", paste(vars_spill, collapse = " + "), "+\n")
    cat("           (1 +", resp_var, "* (", paste(vars, collapse = " + "), ") || celltype) +\n")
    cat("          ", img_re_str, "\n")

    # Step 3: Run responder extraction
    results <- extractRandomEffectsResponder(
        genes = genes,
        df = df,
        vars = vars,
        vars_spill = vars_spill,
        include_spillover = include_spillover,
        resp_var = resp_var,
        resp_level = resp_level,
        family = family,
        BPPARAM = BPPARAM,
        use_offset = use_offset,
        area_var = area_var,
        image_slope_vars = image_slope_vars,
        ...
    )

    results
}


#' Add per-cell total transcript count column for library-size offset
#'
#' Computes nCount as the row-sum of all gene columns in df. Required when
#' use_offset = TRUE in extractRandomEffects. Filters out cells with zero
#' total counts to avoid log(0) in the offset.
#'
#' @param df Data frame with gene columns
#' @param genes Character vector of gene column names to sum
#' @return df with nCount column added, zero-count cells removed
add_nCount <- function(df, genes) {
    gene_cols <- intersect(genes, names(df))
    if (length(gene_cols) == 0L) {
        stop("None of the requested gene columns are in df. Cannot compute nCount.")
    }
    n_count <- rowSums(as.matrix(df[, gene_cols, drop = FALSE]))
    df$nCount <- as.numeric(n_count)
    n_zero <- sum(df$nCount <= 0, na.rm = TRUE)
    if (n_zero > 0L) {
        message(sprintf(
            "  Dropping %d cells with zero total transcript count (cannot take log).",
            n_zero
        ))
        df <- df[df$nCount > 0, , drop = FALSE]
    }
    df
}

