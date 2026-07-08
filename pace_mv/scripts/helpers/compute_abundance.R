#' Compute neighborhood abundance for single image (count-based)
#' 
#' Used for multi-image/responder analysis.
#' 
#' @param df_img Data frame for one image (must have cell_id, cell_type, x, y)
#' @param max_dist Maximum distance
#' @param edge Edge correction method
#' @return Matrix of neighbor counts (cells x cell types)
#' @export
compute_neighborhood_abundance <- function(df_img, max_dist, 
                                           edge = c("none", "border", "isotropic")) {
    edge <- match.arg(edge)
    
    n <- nrow(df_img)
    ct <- factor(df_img$cell_type)
    ct_levels <- levels(ct)
    
    ow <- spatstat.geom::owin(xrange = range(df_img$x), yrange = range(df_img$y))
    ppp <- spatstat.geom::ppp(x = df_img$x, y = df_img$y, window = ow)
    cp <- spatstat.geom::closepairs(ppp, rmax = max_dist, what = "ijd", distinct = TRUE)
    
    if (length(cp$i) == 0) {
        mat <- matrix(0, nrow = n, ncol = length(ct_levels),
                      dimnames = list(df_img$cell_id, ct_levels))
        return(mat)
    }
    
    ct_id <- as.integer(ct)
    src_type <- ct_id[cp$j]
    
    mat <- matrix(0, nrow = n, ncol = length(ct_levels),
                  dimnames = list(df_img$cell_id, ct_levels))
    
    edges_by_type <- split(seq_along(src_type), src_type)
    for (k in names(edges_by_type)) {
        idx <- edges_by_type[[k]]
        tmp <- rowsum(rep(1, length(idx)), group = cp$i[idx], reorder = FALSE)
        mat[as.integer(rownames(tmp)), as.integer(k)] <- tmp[, 1]
    }
    
    if (edge == "isotropic") {
        p <- borderEdge(ppp, max_dist)
        mat <- sweep(mat, 1, p, "/")
    } else if (edge == "border") {
        inside <- spatstat.geom::bdist.points(ppp) >= max_dist
        mat[!inside, ] <- NA
    }
    
    mat
}

#' Compute kernel-weighted neighbourhood abundance for single image
#'
#' Per-cell sum of kernel-weighted contributions from each celltype within
#' max_dist. Replaces the binary indicator I(d < max_dist) used by
#' compute_neighborhood_abundance with a smooth weight w(d) = K(d/h).
#' Direct contact (small d) contributes ~1; distant cells contribute less.
#'
#' @param df_img cell-level df for one image (cell_id, cell_type, x, y)
#' @param max_dist truncation radius (cells beyond this contribute 0)
#' @param h kernel bandwidth (in same units as x, y)
#' @param kernel "exp" -> exp(-d/h); "gaussian" -> exp(-(d/h)^2)
#' @param edge edge-correction method (passed through)
#' @return cells x celltypes matrix of kernel-weighted abundance
#' @export
compute_kernel_abundance <- function(df_img, max_dist = 50, h = 15,
                                       kernel = c("exp", "gaussian"),
                                       edge = c("none", "border",
                                                 "isotropic")) {
  kernel <- match.arg(kernel)
  edge   <- match.arg(edge)

  n <- nrow(df_img)
  ct <- factor(df_img$cell_type)
  ct_levels <- levels(ct)

  ow  <- spatstat.geom::owin(xrange = range(df_img$x),
                              yrange = range(df_img$y))
  ppp <- spatstat.geom::ppp(x = df_img$x, y = df_img$y, window = ow)
  cp  <- spatstat.geom::closepairs(ppp, rmax = max_dist,
                                    what = "ijd", distinct = TRUE)
  if (length(cp$i) == 0L) {
    mat <- matrix(0, nrow = n, ncol = length(ct_levels),
                  dimnames = list(df_img$cell_id, ct_levels))
    return(mat)
  }

  ## Kernel weight per pair
  w <- if (kernel == "exp") exp(-cp$d / h) else exp(-(cp$d / h)^2)

  ct_id    <- as.integer(ct)
  src_type <- ct_id[cp$j]

  mat <- matrix(0, nrow = n, ncol = length(ct_levels),
                dimnames = list(df_img$cell_id, ct_levels))

  edges_by_type <- split(seq_along(src_type), src_type)
  for (k in names(edges_by_type)) {
    idx <- edges_by_type[[k]]
    tmp <- rowsum(w[idx], group = cp$i[idx], reorder = FALSE)
    mat[as.integer(rownames(tmp)), as.integer(k)] <- tmp[, 1]
  }

  if (edge == "isotropic") {
    p   <- borderEdge(ppp, max_dist)
    mat <- sweep(mat, 1, p, "/")
  } else if (edge == "border") {
    inside <- spatstat.geom::bdist.points(ppp) >= max_dist
    mat[!inside, ] <- NA
  }
  mat
}


#' Kernel-weighted neighbourhood abundance over multiple images
#'
#' Mirror of spillover_by_image() but with smooth kernel weighting.
#' @export
kernel_spillover_by_image <- function(cell_data, max_dist = 50, h = 15,
                                       kernel = c("exp", "gaussian"),
                                       edge = c("none", "border",
                                                 "isotropic")) {
  kernel <- match.arg(kernel)
  edge   <- match.arg(edge)

  image_ids <- unique(cell_data$imageID)
  all_cell_types <- sort(unique(cell_data$cell_type))
  results <- list()

  for (img_id in image_ids) {
    df_img <- cell_data[cell_data$imageID == img_id, , drop = FALSE]
    if (nrow(df_img) > 1L) {
      mat <- compute_kernel_abundance(df_img, max_dist = max_dist,
                                       h = h, kernel = kernel, edge = edge)
      rownames(mat) <- df_img$cell_id
      missing_cols <- setdiff(all_cell_types, colnames(mat))
      if (length(missing_cols) > 0L) {
        zero_mat <- matrix(0, nrow = nrow(df_img),
                           ncol = length(missing_cols),
                           dimnames = list(df_img$cell_id, missing_cols))
        mat <- cbind(mat, zero_mat)
      }
      mat <- mat[, all_cell_types, drop = FALSE]
      results[[img_id]] <- mat
    }
  }
  if (length(results) > 0L) do.call(rbind, results)
  else matrix(0, 0, length(all_cell_types),
              dimnames = list(NULL, all_cell_types))
}


#' Border edge correction
#' 
#' Compute edge correction weights based on fraction of disk inside window.
#' 
#' @param X Point pattern (ppp object)
#' @param maxD Maximum distance (radius)
#' @return Vector of edge correction weights
#' @export
borderEdge <- function(X, maxD) {
    W <- spatstat.geom::Window(X)
    near <- spatstat.geom::bdist.points(X) < maxD
    e <- rep(1, spatstat.geom::npoints(X))
    if (any(near)) {
        circs <- spatstat.geom::discs(X[near], maxD, separate = TRUE)
        circs <- spatstat.geom::solapply(circs, spatstat.geom::intersect.owin, W)
        areas <- vapply(circs, spatstat.geom::area, numeric(1)) / (pi * maxD^2)
        e[near] <- areas
    }
    e
}

#' Compute spillover per image (count-based, for multi-image/responder)
#' 
#' @param cell_data Data frame with cell_id, imageID, cell_type, x, y
#' @param radius Maximum distance for neighbors
#' @param edge Edge correction method
#' @return Matrix of neighbor counts (cells x cell types)
#' @export
#' Compute neighbourhood abundance for one image, ANNULAR (inner < d <= outer).
#'
#' Same as compute_neighborhood_abundance but counts only neighbours whose
#' distance falls in the (inner, outer] shell. Used together with a
#' separate fixed-effect spillover term at radius = inner so that the
#' two designs are orthogonal (every neighbour cell appears in exactly one
#' bin: the inner ball or the inner-to-outer shell).
#'
#' @param df_img Data frame for one image (must have cell_id, cell_type, x, y)
#' @param inner Inner radius (excluded)
#' @param outer Outer radius (included)
#' @param edge Edge correction method
#' @return Matrix of neighbour counts (cells x cell types)
#' @export
compute_neighborhood_annular <- function(df_img, inner, outer,
                                          edge = c("none", "border", "isotropic")) {
    edge <- match.arg(edge)
    n <- nrow(df_img)
    ct <- factor(df_img$cell_type)
    ct_levels <- levels(ct)

    ow <- spatstat.geom::owin(xrange = range(df_img$x), yrange = range(df_img$y))
    ppp <- spatstat.geom::ppp(x = df_img$x, y = df_img$y, window = ow)
    cp <- spatstat.geom::closepairs(ppp, rmax = outer, what = "ijd", distinct = TRUE)

    keep <- cp$d > inner   # exclude inner ball; (inner, outer] shell
    if (sum(keep) == 0) {
        return(matrix(0, nrow = n, ncol = length(ct_levels),
                      dimnames = list(df_img$cell_id, ct_levels)))
    }
    cp$i <- cp$i[keep]; cp$j <- cp$j[keep]; cp$d <- cp$d[keep]
    ct_id <- as.integer(ct)
    src_type <- ct_id[cp$j]
    mat <- matrix(0, nrow = n, ncol = length(ct_levels),
                  dimnames = list(df_img$cell_id, ct_levels))
    edges_by_type <- split(seq_along(src_type), src_type)
    for (k in names(edges_by_type)) {
        idx <- edges_by_type[[k]]
        tmp <- rowsum(rep(1, length(idx)), group = cp$i[idx], reorder = FALSE)
        mat[as.integer(rownames(tmp)), as.integer(k)] <- tmp[, 1]
    }
    if (edge == "isotropic") {
        ## Approximate: scale by border-edge correction at the OUTER radius
        ## (the larger window penalty is the binding constraint for shells).
        p <- borderEdge(ppp, outer)
        mat <- sweep(mat, 1, p, "/")
    } else if (edge == "border") {
        inside <- spatstat.geom::bdist.points(ppp) >= outer
        mat[!inside, ] <- NA
    }
    mat
}

#' Annular spillover across images (mirrors spillover_by_image).
#'
#' @export
spillover_annular_by_image <- function(cell_data, inner, outer,
                                        edge = c("none", "border", "isotropic")) {
    edge <- match.arg(edge)
    image_ids <- unique(cell_data$imageID)
    all_cell_types <- sort(unique(cell_data$cell_type))
    results <- list()
    for (img_id in image_ids) {
        df_img <- cell_data[cell_data$imageID == img_id, , drop = FALSE]
        if (nrow(df_img) > 1) {
            mat <- compute_neighborhood_annular(df_img, inner = inner,
                                                outer = outer, edge = edge)
            rownames(mat) <- df_img$cell_id
            missing_cols <- setdiff(all_cell_types, colnames(mat))
            if (length(missing_cols) > 0) {
                for (col in missing_cols) {
                    mat <- cbind(mat, rep(0, nrow(mat)))
                    colnames(mat)[ncol(mat)] <- col
                }
            }
            mat <- mat[, all_cell_types, drop = FALSE]
        } else {
            mat <- matrix(0, nrow = nrow(df_img), ncol = length(all_cell_types),
                          dimnames = list(df_img$cell_id, all_cell_types))
        }
        results[[img_id]] <- mat
    }
    do.call(rbind, results)
}

spillover_by_image <- function(cell_data, radius,
                               edge = c("none", "border", "isotropic")) {
    edge <- match.arg(edge)
    
    image_ids <- unique(cell_data$imageID)
    all_cell_types <- sort(unique(cell_data$cell_type))
    results <- list()
    
    for (img_id in image_ids) {
        df_img <- cell_data[cell_data$imageID == img_id, , drop = FALSE]
        
        if (nrow(df_img) > 1) {
            mat <- compute_neighborhood_abundance(df_img, max_dist = radius, edge = edge)
            rownames(mat) <- df_img$cell_id
            
            if (is.null(dim(mat)) || ncol(mat) == 0) {
                mat <- matrix(0, nrow = nrow(df_img), ncol = length(all_cell_types),
                              dimnames = list(df_img$cell_id, all_cell_types))
            } else {
                missing_cols <- setdiff(all_cell_types, colnames(mat))
                if (length(missing_cols) > 0) {
                    for (col in missing_cols) {
                        mat <- cbind(mat, rep(0, nrow(mat)))
                        colnames(mat)[ncol(mat)] <- col
                    }
                }
                mat <- mat[, all_cell_types, drop = FALSE]
            }
            
            results[[img_id]] <- mat
        }
    }
    
    combined <- do.call(rbind, results)
    combined
}