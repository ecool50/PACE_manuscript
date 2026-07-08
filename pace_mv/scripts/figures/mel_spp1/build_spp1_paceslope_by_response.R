# =============================================================================
# Mel SPP1 macrophage <- Tumour: binned-means figure with PACE slopes overlaid
# -----------------------------------------------------------------------------
# Population view of SPP1 in macrophages vs tumour-neighbour density, split by
# ICI response (PD vs non-PD, the SIMVI contrast, all 26 images).
#
#   solid line + ribbon : binned macrophage means +/- SE (the observed data)
#   dashed line         : the PACE slope, read straight from the fitted model
#                         (mv$fit$U) -- NOT an OLS / geom_smooth fit.
#
# The PACE slope is a slope per unit tumour density on the latent log-mean
# scale. log1p(CP10K) ~= that scale for non-tiny values, so the dashed line is
# drawn with the PACE slope value, anchored at each group's data centroid.
# Output: plots/mel_spp1/spp1_macrophage_by_response_pdnonpd.pdf
# =============================================================================

setwd("/Users/elijahwillie/Documents/Academic/PhD/lateral_spillover/pace_mv")

suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(zellkonverter)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
})

canonical_fit_path <- "data/simvi_melanoma/sweeps/mvpql_percell_hc.rds"
h5ad_path          <- "data/simvi_melanoma/Melanoma_5612.h5ad"
output_pdf         <- "plots/mel_spp1/spp1_macrophage_by_response_pdnonpd.pdf"

cell_types <- c("Tumour", "Endothelial", "Fibroblast", "Macrophage",
                "B_Cell", "T_CD8_memory", "Treg")
tumour_subtypes <- c("Tumor_a", "Tumor_b", "Tumor_c", "Tumor_d", "Tumor_e", "Tumor_f")

mv <- readRDS(canonical_fit_path)

# -----------------------------------------------------------------------------
# Reconstruct the exact cell table the canonical fit was built from, so the
# kernel density (mv$K_bio) and the fit metadata align row-for-row.
# -----------------------------------------------------------------------------
sce <- readH5AD(h5ad_path, reader = "R")
assayNames(sce)[assayNames(sce) == "X"] <- "counts"

# keep only the four RECIST categories, then recode PD vs non-PD (SIMVI contrast)
sce <- sce[, which(sce$BEST_RESPONSE_BY_SCAN %in% c("PD", "SD", "CR", "PR"))]
sce$BEST_RESPONSE_BY_SCAN <- ifelse(as.character(sce$BEST_RESPONSE_BY_SCAN) == "PD",
                                    "PD", "nonPD")

cell_table <- cbind(colData(sce),
                    t(as.matrix(assay(sce, "counts"))),
                    reducedDim(sce, "spatial")) |>
  as.data.frame() |>
  mutate(imageID  = paste(SPID, fov, sep = "_"),
         celltype = make.names(celltype))

cell_table$celltype <- dplyr::recode(cell_table$celltype,
                                     `B.cell`        = "B_Cell",
                                     `T.CD8.memory`  = "T_CD8_memory",
                                     endothelial     = "Endothelial",
                                     fibroblast      = "Fibroblast",
                                     macrophage      = "Macrophage")
cell_table$celltype <- ifelse(cell_table$celltype %in% tumour_subtypes,
                              "Tumour", cell_table$celltype)
cell_table <- cell_table |> filter(celltype %in% cell_types)
cell_table <- cell_table[order(cell_table$imageID), ]

# gene filter: max-per-celltype detection >= 5%, then library size + drop empties
panel_genes <- make.names(rownames(sce))
gene_cols   <- intersect(panel_genes, names(cell_table))
ct_vector   <- as.character(cell_table$celltype)
ct_index    <- lapply(cell_types, function(ct) which(ct_vector == ct))

max_detection <- vapply(gene_cols, function(g) {
  is_expressed <- as.numeric(cell_table[[g]]) > 0
  max(vapply(ct_index,
             function(idx) if (!length(idx)) 0 else mean(is_expressed[idx]),
             numeric(1)))
}, numeric(1))

gene_cols <- gene_cols[max_detection >= 0.05]
cell_table$nCount <- rowSums(as.matrix(cell_table[, gene_cols, drop = FALSE]))
cell_table <- cell_table[cell_table$nCount > 0, ]

# alignment gate: the reconstructed table must match the fit's cell metadata
stopifnot(nrow(cell_table) == nrow(mv$cell_meta),
          max(abs(cell_table$nCount - mv$cell_meta$nCount)) == 0)

# -----------------------------------------------------------------------------
# Macrophage SPP1 expression vs tumour-neighbour density
# -----------------------------------------------------------------------------
macrophage_rows <- which(cell_table$celltype == "Macrophage")
mac <- data.frame(
  tumour_density = mv$K_bio[macrophage_rows, "Tumour"],
  spp1_cp10k     = log1p(1e4 * cell_table[["SPP1"]][macrophage_rows] /
                           cell_table$nCount[macrophage_rows]),
  Response       = factor(mv$cell_meta$Condition[macrophage_rows],
                          levels = c("PD", "nonPD"))
)

# solid line: binned means +/- SE (10 density quantile bins, drop sparse bins)
density_breaks <- unique(quantile(mac$tumour_density, seq(0, 1, length.out = 11),
                                  na.rm = TRUE))
mac$density_bin <- cut(mac$tumour_density, breaks = density_breaks,
                       include.lowest = TRUE)
binned_means <- mac |>
  group_by(Response, density_bin) |>
  summarise(x        = mean(tumour_density),
            mean_spp1 = mean(spp1_cp10k),
            se        = sd(spp1_cp10k) / sqrt(n()),
            n         = n(),
            .groups   = "drop") |>
  filter(n >= 10)

# dashed line: PACE slope straight from the fitted model (no re-fit, no OLS)
#   Macrophage::Tumour                 = SPP1 macrophage-focal slope (non-PD reference)
#   Macrophage::ResponderPD:Tumour     = PD modulation of that slope
slope_blups <- mv$fit$U[, "SPP1"]
names(slope_blups) <- rownames(mv$fit$U)
slope_nonpd <- unname(slope_blups["Macrophage::Tumour"])
slope_pd    <- slope_nonpd + unname(slope_blups["Macrophage::ResponderPD:Tumour"])
pace_slope  <- c(nonPD = slope_nonpd, PD = slope_pd)

# anchor each PACE slope at its group's data centroid, draw across the panel
panel_x_range <- c(0, 65)
pace_lines <- mac |>
  group_by(Response) |>
  summarise(density_mean = mean(tumour_density),
            spp1_mean     = mean(spp1_cp10k),
            .groups       = "drop") |>
  rowwise() |>
  mutate(segment = list(data.frame(
    x = panel_x_range,
    y = spp1_mean + pace_slope[[as.character(Response)]] * (panel_x_range - density_mean)
  ))) |>
  tidyr::unnest(segment) |>
  ungroup()

# -----------------------------------------------------------------------------
# Plot
# -----------------------------------------------------------------------------
response_colours <- c(PD = "#D62728", nonPD = "#1F77B4")
response_labels  <- c(PD = "PD (progressive)", nonPD = "non-PD (SD+PR+CR)")

plot_obj <- ggplot() +
  geom_ribbon(data = binned_means,
              aes(x, ymin = mean_spp1 - se, ymax = mean_spp1 + se, fill = Response),
              alpha = 0.18) +
  geom_line(data = binned_means, aes(x, mean_spp1, colour = Response), linewidth = 0.9) +
  geom_point(data = binned_means, aes(x, mean_spp1, colour = Response), size = 1.8) +
  geom_line(data = pace_lines, aes(x, y, colour = Response),
            linetype = "dashed", linewidth = 1.1) +
  scale_colour_manual(values = response_colours, labels = response_labels,
                      name = "Response") +
  scale_fill_manual(values = response_colours, guide = "none") +
  coord_cartesian(xlim = c(0, 65), ylim = c(0, 4.8)) +
  labs(
    title = "SPP1 in macrophages vs tumour-neighbour density, by ICI response",
    subtitle = paste0("Solid = binned means +/- SE (the data);  dashed = PACE ",
                      "shrunken slope (PD -0.060, non-PD +0.007; interaction -0.067, lfsr 0)"),
    x = "Tumour-neighbour density (Gaussian kernel)",
    y = "SPP1 expression (CP10K, log1p)"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top",
        panel.grid.minor = element_blank(),
        plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(size = 9))

ggsave(output_pdf, plot_obj, width = 10, height = 7, device = cairo_pdf)
cat("Wrote", output_pdf, "\n")
cat(sprintf("PACE slopes drawn:  non-PD %+.4f   PD %+.4f\n",
            pace_slope["nonPD"], pace_slope["PD"]))
