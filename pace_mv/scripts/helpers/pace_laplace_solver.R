## pace_laplace_solver.R — TMB Laplace inner solve for PACE PQL.
##
## Per-gene Laplace approximation that replaces PACE's IRLS+ridge WLS solve.
## Keeps the OUTER iteration (tau EM, alpha update, bleed EM) unchanged --
## only the per-gene (beta, u) optimisation switches from first-order Taylor
## (IRLS) to proper Laplace via TMB.

.LAPLACE_DLL <- "pace_laplace_nb1"

## Lazy compile + load the TMB template.  Idempotent.
## Locates helpers/ by searching ./, ../, and the path of this source file.
.find_helpers_dir <- function() {
  candidates <- c(
    file.path(getwd(), "scripts", "helpers"),
    file.path(getwd(), "helpers"),
    file.path(dirname(getwd()), "helpers"),
    ## Path of THIS source file's directory (works for any cwd)
    tryCatch({
      sf <- normalizePath(sys.frames()[[1]]$ofile, mustWork = FALSE)
      if (!is.null(sf) && nzchar(sf)) dirname(sf) else NA_character_
    }, error = function(e) NA_character_)
  )
  for (d in candidates)
    if (!is.na(d) && file.exists(file.path(d, paste0(.LAPLACE_DLL, ".cpp"))))
      return(d)
  NA_character_
}

ensure_laplace_dll <- function() {
  here <- .find_helpers_dir()
  if (is.na(here))
    stop("Cannot locate ", .LAPLACE_DLL, ".cpp in any of: ",
         "./helpers, ../helpers, or sys.frame path. CWD=", getwd())
  src  <- file.path(here, paste0(.LAPLACE_DLL, ".cpp"))
  so   <- file.path(here, paste0(.LAPLACE_DLL, TMB::dynlib("")))
  ## Compile if missing or stale
  if (!file.exists(so) ||
      file.info(src)$mtime > file.info(so)$mtime) {
    if (getOption("pace.laplace.verbose", FALSE))
      cat("  [laplace] compiling TMB template ...\n")
    TMB::compile(src, flags = "-O2")
  }
  ## Load if not already
  if (!(.LAPLACE_DLL %in% sapply(getLoadedDLLs(), `[[`, "name")))
    dyn.load(TMB::dynlib(file.path(here, .LAPLACE_DLL)))
  invisible(TRUE)
}

## Solve ONE gene's (beta, u) via TMB Laplace.
##
##   y          observed counts (n)
##   X          fixed-effect design (n x p)  (dense)
##   Z          sparse random-effect design (n x q)
##   offset     log offset (n)
##   tau_inv    ridge precision per column (q)  (= 1 / tau_g)
##   alpha      NB1 dispersion (scalar, fixed)
##   beta_init  warm start for beta (default = 0)
##   u_init     warm start for u    (default = 0)
##
## Returns list(beta = p, u = q, se_beta = p, se_u = q, conv = bool).
pace_laplace_one_gene <- function(y, X, Z, offset, tau_inv, alpha,
                                    beta_init = NULL, u_init = NULL,
                                    return_var = TRUE,
                                    maxit = 200L) {
  ensure_laplace_dll()
  if (!inherits(Z, "dgTMatrix") && !inherits(Z, "dgCMatrix"))
    Z <- methods::as(Z, "CsparseMatrix")
  p <- ncol(X); q <- ncol(Z)
  if (is.null(beta_init)) beta_init <- rep(0, p)
  if (is.null(u_init))    u_init    <- rep(0, q)
  alpha_safe <- max(alpha, 1e-4)
  tau_inv_safe <- pmin(pmax(tau_inv, 0), 1e8)

  obj <- TMB::MakeADFun(
    data = list(y = as.numeric(y), X = as.matrix(X), Z = Z,
                 offset = as.numeric(offset),
                 tau_inv = as.numeric(tau_inv_safe),
                 alpha = alpha_safe),
    parameters = list(beta = as.numeric(beta_init),
                       u    = as.numeric(u_init)),
    random = "u",
    DLL = .LAPLACE_DLL,
    silent = TRUE)

  opt <- tryCatch(
    stats::nlminb(obj$par, obj$fn, obj$gr,
                   control = list(iter.max = maxit, eval.max = maxit * 2,
                                   rel.tol = 1e-7)),
    error = function(e) list(convergence = 99L))

  ## Extract converged values BEFORE we discard the obj
  pl <- obj$env$parList(obj$env$last.par)
  beta_hat <- as.numeric(pl$beta)
  u_hat    <- as.numeric(pl$u)
  ## Release obj's internal state; obj will go out of scope at function exit
  ## but explicit nuke makes the per-gene gc() actually reclaim the tape.
  on.exit({ rm(obj); invisible(gc(verbose = FALSE, full = FALSE)) },
           add = TRUE)

  ## Compute diag(var(u)) and diag(var(beta)) ANALYTICALLY (cheaper than sdreport):
  ##   At the mode, marginal Laplace variance of joint params is inv(H_joint).
  ##   H = [H_bb H_bu; H_ub H_uu]
  ##   For NB1 with log link: W_i = mu_i / (1 + alpha)
  ##   H_bb = X' diag(W) X
  ##   H_bu = X' diag(W) Z
  ##   H_uu = Z' diag(W) Z + diag(tau_inv)
  ##
  ## Var(beta) = inv(H_bb - H_bu H_uu^-1 H_ub)
  ## Var(u)    = inv(H_uu - H_ub H_bb^-1 H_bu)  (Schur complement)
  ## For diagonal-only, we use the cheaper computation below.
  re_var <- rep(0, q)
  beta_var <- rep(0, p)
  if (return_var) {
    eta <- as.numeric(X %*% beta_hat + Z %*% u_hat + offset)
    mu  <- exp(eta)
    W   <- mu / (1 + alpha_safe)
    sqW <- sqrt(W)
    ## H_uu (sparse) + diag(tau_inv).  Row-scale Z by sqrt(W) explicitly.
    Z_w  <- Matrix::Diagonal(x = sqW) %*% Z
    H_uu <- methods::as(Matrix::crossprod(Z_w), "CsparseMatrix")
    Matrix::diag(H_uu) <- Matrix::diag(H_uu) + tau_inv_safe
    ## H_bb (dense, p × p; canonical PACE has p = 1 or 2 so this is trivial)
    X_w  <- X * sqW
    H_bb <- crossprod(X_w)
    H_bu <- crossprod(X_w, Z_w)   ## p × q
    rm(Z_w, X_w, sqW, eta, mu)
    ## Schur for diag var(u): diag of inv( H_uu - H_ub H_bb^-1 H_bu )
    H_bb_inv <- tryCatch(solve(H_bb), error = function(e) MASS::ginv(H_bb))
    Schur_u  <- H_uu - Matrix::crossprod(H_bu, H_bb_inv %*% H_bu)
    rm(H_uu, H_bu)
    diag_inv <- tryCatch(Matrix::diag(Matrix::solve(Schur_u)),
                          error = function(e) rep(NA_real_, q))
    re_var <- pmax(as.numeric(diag_inv), 0)
    beta_var <- pmax(diag(H_bb_inv), 0)
    rm(Schur_u, H_bb_inv, H_bb, diag_inv, W)
  }

  list(beta = beta_hat, u = u_hat,
       re_var = re_var, beta_var = beta_var,
       conv = identical(opt$convergence, 0L))
}

## Batch solver: process a chunk of genes with shared X, Z, offset matrix.
##
## Memory-conscious: serial by default (one TMB obj at a time, freed after use).
## Set R_LAPLACE_WORKERS > 1 to enable MulticoreParam parallelism (forks the
## parent — each worker holds its own copy of TMB state ~200 MB so cap workers
## carefully on memory-constrained machines).
pace_laplace_chunk <- function(Y_chunk, X, Z, offset_mat_chunk,
                                tau_inv_chunk, alpha_chunk,
                                beta_init_mat = NULL, u_init_mat = NULL,
                                return_var = TRUE,
                                BPPARAM = NULL) {
  ensure_laplace_dll()
  n_g <- ncol(Y_chunk)
  p   <- ncol(X); q <- ncol(Z)
  if (is.null(beta_init_mat)) beta_init_mat <- matrix(0, p, n_g)
  if (is.null(u_init_mat))    u_init_mat    <- matrix(0, q, n_g)

  wkrs <- as.integer(Sys.getenv("R_LAPLACE_WORKERS", unset = "1"))
  wkrs <- max(1L, wkrs)

  solve_one <- function(gi) {
    out <- pace_laplace_one_gene(
      y         = Y_chunk[, gi],
      X         = X,
      Z         = Z,
      offset    = offset_mat_chunk[, gi],
      tau_inv   = tau_inv_chunk[, gi],
      alpha     = alpha_chunk[gi],
      beta_init = beta_init_mat[, gi],
      u_init    = u_init_mat[, gi],
      return_var = return_var)
    out
  }

  if (wkrs == 1L) {
    ## Serial loop with **per-gene GC** to release TMB autodiff tapes.
    ## Without per-gene gc, TMB objects accumulate (~500 MB tape each) before
    ## R's automatic GC fires, blowing memory through 50 GB on long runs.
    ## Per-gene gc costs ~0.1s but caps peak memory to a single tape (~1 GB).
    res_list <- vector("list", n_g)
    for (gi in seq_len(n_g)) {
      res_list[[gi]] <- solve_one(gi)
      ## explicit release of the TMB autodiff tape
      invisible(gc(verbose = FALSE, full = FALSE))
    }
    return(res_list)
  }

  ## Parallel path: explicit MulticoreParam with capped workers.
  bp <- BiocParallel::MulticoreParam(workers = wkrs, RNGseed = 1L)
  BiocParallel::bplapply(seq_len(n_g), solve_one, BPPARAM = bp)
}
