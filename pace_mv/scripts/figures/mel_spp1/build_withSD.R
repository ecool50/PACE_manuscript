## Mel SPP1-in-macrophages vs tumour-neighbour-density, all 25 patients, binary response (PD+SD vs PR+CR)
## Computed directly from h5ad. Validated against mac_spp1_joined.rds (18-patient PD/RESP).
suppressMessages({
  library(zellkonverter); library(SingleCellExperiment)
  library(dbscan); library(Matrix); library(dplyr); library(ggplot2)
})
set.seed(1)
setwd("/Users/elijahwillie/Documents/Academic/PhD/lateral_spillover/pace_mv"); fig_dir <- "plots/mel_spp1"; cache_dir <- "data/figure_cache"
px_size <- 0.12028
h_bio   <- 30
radius  <- 90

## ---- 1. LOAD ----
sce <- readH5AD("data/simvi_melanoma/Melanoma_5612.h5ad", reader="R")
assayNames(sce)[assayNames(sce) == "X"] <- "counts"
rownames(sce) <- make.names(rownames(sce))
stopifnot("SPP1" %in% rownames(sce))

cd <- as.data.frame(colData(sce))
keep <- cd$BEST_RESPONSE_BY_SCAN %in% c("PD","SD","PR","CR")
sce <- sce[, keep]; cd <- cd[keep, ]

## collapse tumour subtypes
ct <- as.character(cd$celltype)
ct[grepl("^Tumor_", ct)] <- "Tumour"
cd$celltype_c <- ct

## imageID, coords (microns), binary response
cd$imageID  <- paste(cd$SPID, cd$fov, sep="_")
sp <- reducedDim(sce, "spatial")
cd$x <- sp[,1] * px_size
cd$y <- sp[,2] * px_size
resp_chr <- as.character(cd$BEST_RESPONSE_BY_SCAN)
cd$Response <- ifelse(resp_chr %in% c("PD","SD"), "NonResp", "Resp")
cd$Response <- factor(cd$Response, levels = c("NonResp","Resp"))

## ---- 4. SPP1 CP10K (over all genes present) ----
cnt <- assay(sce, "counts")
nCount <- Matrix::colSums(cnt)
spp1_count <- as.numeric(cnt["SPP1", ])
cd$nCount     <- nCount
cd$SPP1_count <- spp1_count
cd$SPP1_cp10k <- log1p(1e4 * spp1_count / nCount)

## ---- 3. TUMOUR-NEIGHBOUR DENSITY per imageID (Gaussian, h_bio=30, radius 90) ----
cd$tumour_density <- NA_real_
images <- unique(cd$imageID)
for (im in images) {
  idx <- which(cd$imageID == im)
  coords_im <- as.matrix(cd[idx, c("x","y")])
  is_tum <- cd$celltype_c[idx] == "Tumour"
  nn <- frNN(coords_im, eps = radius)
  dens <- numeric(length(idx))
  for (j in seq_along(idx)) {
    nb <- nn$id[[j]]                 # neighbours within radius (self excluded by frNN)
    if (length(nb) == 0) { dens[j] <- 0; next }
    nb_tum <- nb[is_tum[nb]]
    if (length(nb_tum) == 0) { dens[j] <- 0; next }
    d <- nn$dist[[j]][is_tum[nb]]
    dens[j] <- sum(exp(-(d^2) / (h_bio^2)))
  }
  cd$tumour_density[idx] <- dens
}
stopifnot(!any(is.na(cd$tumour_density)))

## ---- 5. subset to macrophages ----
cd$cell_id   <- cd$cell_ID
cd$patientID <- cd$SPID
mac <- cd[cd$celltype_c == "macrophage",
          c("cell_id","patientID","fov","imageID","Response",
            "BEST_RESPONSE_BY_SCAN","tumour_density","SPP1_cp10k","SPP1_count","nCount")]
mac <- mac %>% group_by(patientID) %>% mutate(n_mac_patient = n()) %>% ungroup()
mac <- as.data.frame(mac)

## ---- VALIDATION GATE ----
cat("\n================ VALIDATION ================\n")
pat_resp <- unique(cd[,c("SPID","BEST_RESPONSE_BY_SCAN")])
cat("Total patients:", nrow(pat_resp), "\n")
cat("By category:\n"); print(table(as.character(pat_resp$BEST_RESPONSE_BY_SCAN)))
pat_bin <- unique(cd[,c("SPID","Response")])
cat("Binary patients:\n"); print(table(as.character(pat_bin$Response)))
cat("Macrophages per arm:\n"); print(table(as.character(mac$Response)))
cat("Total macrophages:", nrow(mac), "\n")

## compare to reference (18-patient)
ref <- readRDS(file.path(cache_dir, "mac_spp1_joined.rds"))
ref$key <- paste(ref$patientID, ref$cell_id, sep="::")
mac$key <- paste(mac$patientID, mac$cell_id, sep="::")
m <- merge(ref, mac, by="key", suffixes=c("_ref","_new"))
cat("\nOverlap cells with reference:", nrow(m), "of", nrow(ref), "ref rows\n")
cat("cor(SPP1_cp10k new vs ref):",
    round(cor(m$SPP1_cp10k_new, m$SPP1_cp10k_ref),4), "\n")
cat("cor(tumour_density new vs ref):",
    round(cor(m$tumour_density_new, m$tumour_density_ref),4), "\n")
mac$key <- NULL

saveRDS(mac, file.path(cache_dir, "mac_spp1_withSD.rds"))
cat("\nSaved mac_spp1_withSD.rds\n")

## ============ PLOTS ============
tab10 <- c(NonResp="#d62728", Resp="#1f77b4")
theme_set(theme_bw(base_size=12))

## ---- E: population binned trend (8 quantile bins per arm) ----
mk_bins <- function(df) {
  df %>% group_by(Response) %>%
    mutate(bin = cut(tumour_density,
                     breaks = quantile(tumour_density, probs=seq(0,1,length.out=9), na.rm=TRUE),
                     include.lowest=TRUE, labels=FALSE)) %>%
    group_by(Response, bin) %>%
    summarise(x = mean(tumour_density),
              spp1 = mean(SPP1_cp10k),
              se = sd(SPP1_cp10k)/sqrt(n()),
              n = n(), .groups="drop")
}
binned <- mk_bins(mac)
pE <- ggplot(binned, aes(x, spp1, colour=Response, group=Response)) +
  geom_line(linewidth=0.9) +
  geom_point(size=2) +
  geom_errorbar(aes(ymin=spp1-se, ymax=spp1+se), width=0, linewidth=0.6) +
  geom_text(aes(label=n), vjust=-1.1, size=2.6, show.legend=FALSE) +
  scale_colour_manual(values=tab10, name="Response",
                      labels=c(NonResp="Non-responder (PD+SD)", Resp="Responder (PR+CR)")) +
  labs(x = "Tumour-neighbour density (Gaussian, h=30 µm)",
       y = "Macrophage SPP1 (log CP10K)",
       title = "Macrophage SPP1 vs tumour-neighbour density",
       subtitle = "Population binned trend (8 quantile bins per arm), all 25 patients") +
  theme(legend.position="top")
ggsave(file.path(fig_dir, "E_population_binary_withSD.pdf"), pE,
       width=7, height=5.2, device=cairo_pdf)

## ---- F: donor gate (per-patient slopes) ----
slope_tab <- function(min_mac) {
  mac %>% group_by(patientID, Response) %>%
    filter(n() >= min_mac) %>%
    summarise(slope = coef(lm(SPP1_cp10k ~ tumour_density))[2],
              n_mac = n(), .groups="drop")
}
sl3  <- slope_tab(3)
sl30 <- slope_tab(30)

wilc <- function(s) {
  if (length(unique(s$Response)) < 2) return(NA_real_)
  wilcox.test(slope ~ Response, data=s)$p.value
}
p3  <- wilc(sl3)
p30 <- wilc(sl30)

med_by_arm <- function(s) s %>% group_by(Response) %>%
  summarise(med = median(slope), n_pat = n(), .groups="drop")
cat("\n================ DONOR GATE ================\n")
cat(">=3 macrophages/patient: Wilcoxon p =", round(p3,4), "\n")
print(as.data.frame(med_by_arm(sl3)))
cat("\n>=30 macrophages/patient: Wilcoxon p =", round(p30,4), "\n")
print(as.data.frame(med_by_arm(sl30)))

## SuperPlot (ungated >=3)
grpsum <- sl3 %>% group_by(Response) %>%
  summarise(m = mean(slope),
            lo = m - qt(.975, n()-1)*sd(slope)/sqrt(n()),
            hi = m + qt(.975, n()-1)*sd(slope)/sqrt(n()),
            n_pat = n(), .groups="drop")
set.seed(2)
pF <- ggplot(sl3, aes(Response, slope, colour=Response)) +
  geom_hline(yintercept=0, linetype="dotted", colour="grey50") +
  geom_jitter(aes(size=n_mac), width=0.13, alpha=0.7) +
  geom_point(data=grpsum, aes(y=m), shape=95, size=14, show.legend=FALSE) +
  geom_errorbar(data=grpsum, aes(y=m, ymin=lo, ymax=hi), width=0.12,
                linewidth=0.8, show.legend=FALSE) +
  geom_text(data=grpsum, aes(y=Inf, label=paste0("n=",n_pat)),
            vjust=1.5, size=3.2, show.legend=FALSE) +
  scale_colour_manual(values=tab10, guide="none") +
  scale_x_discrete(labels=c(NonResp="Non-responder\n(PD+SD)", Resp="Responder\n(PR+CR)")) +
  scale_size_continuous(name="n macrophages", range=c(1.5,7)) +
  annotate("text", x=1.5, y=max(sl3$slope), vjust=1,
           label=paste0("Wilcoxon p = ", signif(p3,3)), size=3.4) +
  labs(x=NULL, y="Per-patient slope  (SPP1 ~ tumour density)",
       title="Donor gate: per-patient SPP1–tumour-density slopes",
       subtitle="≥3 macrophages/patient; point size = n macrophages") +
  theme(legend.position="right")
ggsave(file.path(fig_dir, "F_donor_gate_binary_withSD.pdf"), pF,
       width=7, height=5.5, device=cairo_pdf)

cat("\nSaved E_population_binary_withSD.pdf and F_donor_gate_binary_withSD.pdf\n")
cat("DONE\n")
