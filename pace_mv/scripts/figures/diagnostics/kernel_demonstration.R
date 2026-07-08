## kernel_demonstration.R
## Pedagogical figures for the primer: show how the bio (Gaussian h=30) and
## tech (exponential h=5) kernels behave on real BC and Mel data.
##
## Three panels per cohort:
##   A. Kernel weight vs distance, overlaid on the empirical pairwise-distance
##      density from the cohort (within-image, sub-50 µm).
##   B. Per-cell K^bio distribution for the flagship focal->neighbour pair
##      (BC: Stromal -> Tumour; Mel: Macrophage -> Tumour), with K^tech for the
##      same pair on a second panel below.
##   C. K^bio vs hard 50-µm count for the same pair (kernel-weighted vs
##      step-function count), demonstrating that the kernel is a smoothed count.
##
## Output: pace_mv/plots/kernel_demo_bc.pdf, kernel_demo_mel.pdf

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(dplyr); library(FNN)
  library(SingleCellExperiment); library(zellkonverter)
})

setwd("/Users/elijahwillie/Documents/Academic/PhD/lateral_spillover/pace_mv")

H_BIO  <- 30
H_TECH <- 5
PALETTE <- c(bio = "#1F77B4", tech = "#D62728", hard = "#666666")

theme_primer <- theme_classic(base_size = 10) +
  theme(plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 9, colour = "grey30"),
        legend.position = "right",
        panel.grid.major.y = element_line(colour = "grey92"))

## -------- Panel A: kernel curves + empirical distance density ---------
##  NN distances are query=sample, data=FULL within each image, so the
##  reported distance is the true nearest-neighbour distance, not a
##  sparse-subsample inflation.
make_panel_A <- function(coords_per_image, cohort) {
  set.seed(1)
  dists <- unlist(lapply(coords_per_image, function(co) {
    n <- nrow(co); if (n < 50) return(NULL)
    s <- sample.int(n, min(2000, n))
    nn <- FNN::get.knnx(data = co, query = co[s, , drop = FALSE], k = 21L)
    ## drop column 1 (self) so we have k=20 true neighbours
    as.numeric(nn$nn.dist[, -1, drop = FALSE])
  }))
  dists <- dists[dists <= 100 & dists > 0]
  cat(sprintf("[%s] panel A: %s NN distances (k=20, true within-image) <=100µm\n",
              cohort, format(length(dists), big.mark = ",")))

  curves <- data.frame(d = seq(0, 100, length.out = 401)) |>
    mutate(`K^bio (Gauss h=30)`  = exp(-d^2 / H_BIO^2),
           `K^tech (Exp h=5)`    = exp(-d   / H_TECH),
           `Hard 50 µm`          = as.numeric(d < 50)) |>
    tidyr::pivot_longer(-d, names_to = "kernel", values_to = "w") |>
    mutate(kernel = factor(kernel, levels = c("K^bio (Gauss h=30)",
                                              "K^tech (Exp h=5)",
                                              "Hard 50 µm")))

  cols <- c("K^bio (Gauss h=30)" = unname(PALETTE["bio"]),
            "K^tech (Exp h=5)"   = unname(PALETTE["tech"]),
            "Hard 50 µm"         = unname(PALETTE["hard"]))
  ltys <- c("K^bio (Gauss h=30)" = "solid",
            "K^tech (Exp h=5)"   = "solid",
            "Hard 50 µm"         = "dashed")

  ggplot() +
    geom_histogram(data = data.frame(d = dists),
                   aes(x = d, y = after_stat(density) * 30),
                   bins = 60, fill = "grey85", colour = "grey55", linewidth = 0.15) +
    geom_line(data = curves, aes(d, w, colour = kernel, linetype = kernel),
              linewidth = 0.95) +
    geom_vline(xintercept = c(H_TECH, H_BIO), linetype = "dotted", colour = "grey30") +
    annotate("text", x = H_TECH + 1, y = 0.40, label = "h_tech=5",
             hjust = 0, vjust = -0.4, size = 3, colour = unname(PALETTE["tech"])) +
    annotate("text", x = H_BIO + 1,  y = 0.40, label = "h_bio=30",
             hjust = 0, vjust = -0.4, size = 3, colour = unname(PALETTE["bio"])) +
    scale_colour_manual(values = cols, name = NULL) +
    scale_linetype_manual(values = ltys, name = NULL) +
    scale_x_continuous(breaks = seq(0, 100, 20)) +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 1.05)) +
    labs(title = paste0("A. Kernel weight vs distance - ", cohort),
         subtitle = "grey histogram: true NN distances (k=20, sampled queries against full image, density rescaled to fit)",
         x = "distance d (µm) from focal cell",
         y = "kernel weight w(d)") +
    theme_primer
}

## -------- Panel B: K^bio + K^tech distributions for one focal celltype ----
make_panel_B <- function(K_bio, K_tech, celltype_vec, focal_ct, cohort) {
  idx <- which(celltype_vec == focal_ct)
  if (length(idx) == 0) stop("focal celltype not found: ", focal_ct)
  cat(sprintf("[%s] panel B: focal=%s, n_cells=%s\n",
              cohort, focal_ct, format(length(idx), big.mark = ",")))

  nbr_types <- colnames(K_bio)
  df_long <- data.frame(
    K = c(as.vector(K_bio[idx, ]), as.vector(K_tech[idx, ])),
    nbr = rep(rep(nbr_types, each = length(idx)), 2),
    kernel = rep(c("K^bio", "K^tech"), each = length(idx) * length(nbr_types))
  )
  df_long$nbr <- factor(df_long$nbr, levels = nbr_types)

  ggplot(df_long, aes(x = K + 1e-3, fill = kernel)) +
    geom_histogram(bins = 40, position = "identity", alpha = 0.7, colour = NA) +
    facet_grid(kernel ~ nbr, scales = "free", switch = "y") +
    scale_x_continuous(trans = "log10",
                       breaks = c(0.01, 0.1, 1, 10),
                       labels = c("0", "0.1", "1", "10")) +
    scale_fill_manual(values = c(`K^bio` = unname(PALETTE["bio"]),
                                  `K^tech` = unname(PALETTE["tech"]))) +
    labs(title = paste0("B. Per-cell kernel exposure for focal = ", focal_ct,
                         " (", cohort, ")"),
         subtitle = sprintf("each panel = %s cells; x = kernel-weighted exposure to a neighbour celltype",
                             format(length(idx), big.mark = ",")),
         x = "kernel-weighted neighbour count (log10 scale, +1e-3)",
         y = "cell count") +
    theme_primer +
    theme(legend.position = "none",
          strip.text.x = element_text(size = 8),
          strip.text.y = element_text(size = 9, face = "bold"),
          strip.background.y = element_blank())
}

## -------- Panel C: K^bio vs hard 50-µm count for flagship pair ----------
##  Recomputes K^bio FRESH without the model's KNN-50 truncation, using all
##  cells out to where the bio kernel weight is negligible (~3 * h_bio = 90 µm).
##  This shows the conceptual smoothed-count relationship; the saved K_bio in
##  the RDS uses KNN_K=50 truncation for tractability, which can suppress K^bio
##  in densely-packed regions (e.g., Mel tumour cores where the 50th nearest
##  cell sits at ~10 µm and tumours beyond it don't contribute to the saved K).
make_panel_C <- function(coords, celltype_vec, image_vec, focal_ct, nbr_ct, cohort,
                          k_search = 400L, hard_radius = 50,
                          bio_cutoff = 90, h_bio = H_BIO) {
  idx <- which(celltype_vec == focal_ct)
  N50    <- integer(length(idx))
  K_true <- numeric(length(idx))
  censored <- logical(length(idx))
  for (im in unique(image_vec)) {
    rows_all <- which(image_vec == im)
    if (length(rows_all) < 5) next
    rows_focal <- intersect(idx, rows_all)
    if (length(rows_focal) == 0) next
    co_all <- coords[rows_all, , drop = FALSE]
    ct_all <- celltype_vec[rows_all]
    co_focal <- coords[rows_focal, , drop = FALSE]
    k_eff <- min(k_search, length(rows_all) - 1L)
    nn <- FNN::get.knnx(data = co_all, query = co_focal, k = k_eff + 1L)
    nbrs_idx <- nn$nn.index[, -1, drop = FALSE]
    nbrs_d   <- nn$nn.dist[,  -1, drop = FALSE]
    for (rr in seq_along(rows_focal)) {
      d_row <- nbrs_d[rr, ]
      ct_row <- ct_all[nbrs_idx[rr, ]]
      sel_n50 <- (d_row < hard_radius) & (ct_row == nbr_ct)
      sel_bio <- (d_row < bio_cutoff)  & (ct_row == nbr_ct)
      pos <- which(idx == rows_focal[rr])
      N50[pos]      <- sum(sel_n50)
      K_true[pos]   <- sum(exp(-(d_row[sel_bio]^2) / (h_bio^2)))
      censored[pos] <- (d_row[length(d_row)] < bio_cutoff)
    }
  }
  d <- data.frame(K = K_true, N50 = N50, censored = censored)
  cat(sprintf("[%s] panel C: n_focal=%d, K_true=[%.2f, %.2f], N50=[%d, %d], %d censored\n",
              cohort, length(idx), min(K_true), max(K_true), min(N50), max(N50), sum(censored)))

  ggplot(d, aes(N50, K)) +
    geom_jitter(width = 0.15, height = 0, alpha = 0.25, size = 0.7,
                colour = unname(PALETTE["bio"])) +
    geom_smooth(method = "loess", se = FALSE, colour = "grey20", linewidth = 0.6) +
    labs(title = sprintf("C. K^bio vs hard 50-µm count - %s -> %s (%s)",
                          focal_ct, nbr_ct, cohort),
         subtitle = sprintf("K^bio recomputed without the model's 50-NN truncation (sum over all %s neighbours within %.0f µm); strictly monotone in hard count",
                             nbr_ct, bio_cutoff),
         x = paste0("hard count: # ", nbr_ct, " cells within 50 µm"),
         y = paste0("K^bio toward ", nbr_ct, " (kernel-weighted)")) +
    theme_primer
}

## ========================= BC ===========================================
cat("\n========== BC ==========\n")
bc_rds <- readRDS("data/breast_cancer/mvpql_canonical.rds")
bc_src <- readRDS("data/breast_cancer/Y_df_for_mcsd.rds")
bc_df  <- bc_src$df
bc_K_bio  <- bc_rds$K_bio
bc_K_tech <- bc_rds$K_tech

## verify alignment between K matrices and df
stopifnot(nrow(bc_K_bio) == nrow(bc_df))

bc_coords <- as.matrix(bc_df[, c("x", "y")])
bc_celltype <- as.character(bc_df$cellType)
bc_image <- rep("BC_single_image", nrow(bc_df))  # BC is one image

pA_bc <- make_panel_A(list(BC = bc_coords), "BC")
pB_bc <- make_panel_B(bc_K_bio, bc_K_tech, bc_celltype, focal_ct = "Stromal", "BC")
pC_bc <- make_panel_C(bc_coords, bc_celltype, bc_image,
                      focal_ct = "Stromal", nbr_ct = "Tumour", "BC")

pdf("plots/diagnostics/kernel_demo_bc.pdf", width = 11, height = 11)
print(pA_bc / pB_bc / pC_bc + plot_layout(heights = c(1, 1.2, 1)))
dev.off()
cat("Saved plots/diagnostics/kernel_demo_bc.pdf\n")

## ========================= Mel ==========================================
cat("\n========== Mel ==========\n")
mel_rds <- readRDS("data/simvi_melanoma/mvpql_canonical.rds")
mel_K_bio  <- mel_rds$K_bio
mel_K_tech <- mel_rds$K_tech

## reload Mel from h5ad to get coords + celltype in same order as K matrices
sce <- readH5AD(file = "data/simvi_melanoma/Melanoma_5612.h5ad", reader = "R")
assayNames(sce)[assayNames(sce) == "X"] <- "counts"
px_size <- 0.12028
sce <- sce[, which(sce$BEST_RESPONSE_BY_SCAN %in% c("PD", "SD"))]
TYPES_MEL <- c("Endothelial", "Fibroblast", "Macrophage", "B_Cell", "Tumour", "T_Cell")
mel_df <- cbind(colData(sce), reducedDim(sce, "spatial")) |>
  as.data.frame() |>
  dplyr::mutate(imageID = SPID, x = V1 * px_size, y = V2 * px_size) |>
  dplyr::mutate(celltype = make.names(celltype))
mel_df$celltype <- dplyr::recode(mel_df$celltype,
  'Tumor_a' = "Tumour", 'Tumor_b' = "Tumour", 'Tumor_c' = "Tumour",
  'Tumor_d' = "Tumour", 'Tumor_e' = "Tumour", 'Tumor_f' = "Tumour",
  'B.cell' = "B_Cell", 'T.CD8.memory' = "T_Cell",
  'endothelial' = "Endothelial", 'fibroblast' = "Fibroblast",
  'macrophage' = "Macrophage", 'T.CD8.naive' = "T_Cell", 'Treg' = "T_Cell")
mel_df <- mel_df |> dplyr::filter(celltype %in% TYPES_MEL)

stopifnot(nrow(mel_K_bio) == nrow(mel_df))

mel_coords <- as.matrix(mel_df[, c("x", "y")])
mel_celltype <- as.character(mel_df$celltype)
mel_image <- as.character(mel_df$imageID)

coords_by_image <- split(seq_len(nrow(mel_df)), mel_image) |>
  lapply(function(rows) mel_coords[rows, ])

pA_mel <- make_panel_A(coords_by_image, "Mel")
pB_mel <- make_panel_B(mel_K_bio, mel_K_tech, mel_celltype, focal_ct = "Macrophage", "Mel")
pC_mel <- make_panel_C(mel_coords, mel_celltype, mel_image,
                       focal_ct = "Macrophage", nbr_ct = "Tumour", "Mel")

pdf("plots/diagnostics/kernel_demo_mel.pdf", width = 11, height = 11)
print(pA_mel / pB_mel / pC_mel + plot_layout(heights = c(1, 1.2, 1)))
dev.off()
cat("Saved plots/diagnostics/kernel_demo_mel.pdf\n")
