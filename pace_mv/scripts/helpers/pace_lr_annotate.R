## pace_lr_annotate.R — ligand-receptor annotation of PACE pair drivers.
##
## Interpretation layer on top of the canonical drivers (NOT a competing LR
## method). For each focal<-neighbour pair, each significant driver gene (a FOCAL
## gene that responds to neighbour density) is checked against CellChatDB: if the
## driver is a ligand or a receptor subunit, we name the cognate partner and
## report the partner's expression in the NEIGHBOUR cell type, so the spatial
## coupling reads as a candidate signaling axis. Because PACE drivers are
## contamination-corrected, an axis built on them is cleaner than raw spatial
## LR co-expression (which is confounded by segmentation spillover).
##
## Direction:
##   driver = receptor -> neighbour(ligand)  ->  focal(driver receptor)
##   driver = ligand   -> focal(driver ligand) ->  neighbour(receptor)
## Either way the cognate partner is evaluated in the NEIGHBOUR (paracrine read).
## Homotypic pairs (ligand gene == receptor gene, e.g. CDH1, PECAM1) are flagged
## as adhesion, not paracrine signaling.

pace_lr_annotate <- function(mv, pairs, lfsr_thresh = 0.05, db = NULL) {
  if (is.null(db)) db <- CellChat::CellChatDB.human
  panel <- mv$gene_set
  fit   <- mv$fit
  shr   <- mv$shrunken_long
  cx    <- db$complex
  inter <- db$interaction

  ## Resolve an interaction side (a gene symbol or a complex name) to its genes.
  resolve_side <- function(x) {
    if (x %in% rownames(cx)) {
      subunits <- as.character(cx[x, ])
      subunits[subunits != "" & !is.na(subunits)]
    } else x
  }
  ligand_genes   <- lapply(inter$ligand,   resolve_side)
  receptor_genes <- lapply(inter$receptor, resolve_side)

  ## Per-cell-type mean fitted expression (contamination-included mu).
  blk   <- which(vapply(fit$re_meta$blocks, `[[`, character(1), "group_col") == "celltype")
  types <- fit$re_meta$blocks[[blk]]$group_levels
  cells_by_type <- setNames(fit$re_meta$cells_by_grp_list[[blk]], types)
  ct_mean <- sapply(types, function(t) Matrix::colMeans(fit$mu[cells_by_type[[t]], , drop = FALSE]))
  rownames(ct_mean) <- panel

  contains_gene <- function(gene_lists, gene) which(vapply(gene_lists, function(v) gene %in% v, logical(1)))

  rows <- list()
  for (p in pairs) {
    focal <- p[[1]]; neighbour <- p[[2]]
    ## Rank drivers within the pair by |shrunken slope| so each LR hit can report
    ## how prominent a PACE driver it is (top-ranked = "in keeping"; tail = marginal).
    drivers <- shr |>
      dplyr::filter(focal == !!focal, neighbour == !!neighbour, lfsr < lfsr_thresh) |>
      dplyr::distinct(gene, .keep_all = TRUE) |>
      dplyr::arrange(dplyr::desc(abs(estimate_shrunk)))
    n_drivers <- nrow(drivers)

    for (k in seq_len(n_drivers)) {
      gene <- drivers$gene[k]
      ligand_i   <- contains_gene(ligand_genes,   gene)
      receptor_i <- contains_gene(receptor_genes, gene)
      hits <- rbind(
        if (length(ligand_i))   data.frame(i = ligand_i,   role = "ligand",   stringsAsFactors = FALSE),
        if (length(receptor_i)) data.frame(i = receptor_i, role = "receptor", stringsAsFactors = FALSE))
      if (is.null(hits)) next
      for (h in seq_len(nrow(hits))) {
        i <- hits$i[h]; role <- hits$role[h]
        partner_genes <- if (role == "ligand") receptor_genes[[i]] else ligand_genes[[i]]
        homotypic     <- setequal(ligand_genes[[i]], receptor_genes[[i]])
        partner_on    <- length(partner_genes) > 0 && all(partner_genes %in% panel)

        ## Direction-aware partner expression: the cognate partner can sit on the
        ## neighbour (paracrine) OR the focal (autocrine), so evaluate BOTH and
        ## take the compartment where it is most expressed. partner_spec = that
        ## expression relative to the partner's max across all cell types, so a
        ## gene barely expressed anywhere in the pair (e.g. MRC1 in B/T cells) is
        ## caught even when it is on-panel.
        if (partner_on) {
          partner_profile <- colMeans(ct_mean[partner_genes, , drop = FALSE])  # per cell type
          pe_nbr    <- partner_profile[[neighbour]]
          pe_foc    <- partner_profile[[focal]]
          pe_target <- max(pe_nbr, pe_foc)
          where     <- if (pe_foc >= pe_nbr) "focal (autocrine)" else "neighbour (paracrine)"
          pspec     <- pe_target / max(partner_profile)
        } else {
          pe_nbr <- pe_foc <- pe_target <- pspec <- NA_real_; where <- NA_character_
        }

        rows[[length(rows) + 1]] <- data.frame(
          focal = focal, neighbour = neighbour,
          driver = gene, role = role,
          partner = paste(partner_genes, collapse = "+"),
          pathway = inter$pathway_name[i],
          category = inter$annotation[i],                 # Secreted Signaling / Cell-Cell Contact / ECM-Receptor
          class   = if (homotypic) "adhesion (homotypic)" else "paracrine",
          partner_on_panel = partner_on,
          partner_expr_neighbour = round(pe_nbr, 3),
          partner_expr_focal     = round(pe_foc, 3),
          partner_expr_target    = round(pe_target, 3),
          partner_where          = where,                 # where the cognate partner actually sits
          partner_spec           = round(pspec, 3),        # target expr / partner max across cell types
          driver_b      = round(drivers$estimate_shrunk[k], 3),
          direction     = if (drivers$estimate_shrunk[k] > 0) "up near nbr" else "down near nbr",
          driver_rank   = sprintf("%d/%d", k, n_drivers),
          driver_pctile = round(100 * (1 - k / n_drivers)),   # 99 = top driver of the pair
          driver_lfsr   = signif(drivers$lfsr[k], 2),
          stringsAsFactors = FALSE)
      }
    }
  }
  if (!length(rows)) return(NULL)
  out <- dplyr::bind_rows(rows)
  out <- out[!duplicated(out[, c("focal", "neighbour", "driver", "role", "partner")]), ]
  out[order(out$focal, out$neighbour, -out$partner_on_panel,
            -tidyr::replace_na(out$partner_expr_neighbour, -1)), ]
}

## Standard ligand-receptor dot/bubble plot (CellPhoneDB / LIANA style) of the
## complete axes: focal<-neighbour pairs on x, ligand->receptor axes on y, dot
## size = PACE driver prominence (percentile), colour = driver slope (blue = down
## / red = up near the neighbour). Bidirectional axes (both ligand and receptor
## are drivers of a pair) are collapsed to the more prominent side.
## Two credibility filters (both adjustable):
##   categories       : keep only these CellChatDB interaction categories.
##                      Default "Secreted Signaling" drops the soft Cell-Cell
##                      Contact set (CD45, adhesion, checkpoint) — e.g. the
##                      spurious PTPRC->MRC1.
##   min_partner_spec : keep only axes where the cognate partner is genuinely
##                      expressed in the focal OR neighbour (specificity =
##                      target expr / partner max across cell types >= this).
##                      Direction-aware, so real autocrine axes (EDN1->EDNRB,
##                      receptor on the focal) survive while receptor-not-in-
##                      target artifacts drop. Default 0.2 is a tunable knob.
pace_lr_dotplot <- function(lr, title = "Ligand-receptor axes behind spatial couplings",
                            categories = "Secreted Signaling", min_partner_spec = 0.2) {
  d <- lr |>
    dplyr::filter(partner_on_panel, class == "paracrine",
                  category %in% categories, partner_spec >= min_partner_spec) |>
    dplyr::mutate(
      ligand   = ifelse(role == "ligand",   driver, partner),
      receptor = ifelse(role == "receptor", driver, partner),
      axis = paste0(ligand, " → ", receptor),
      pair = paste0(focal, " ← ", neighbour)) |>
    dplyr::group_by(pair, axis) |>
    dplyr::slice_max(driver_pctile, n = 1, with_ties = FALSE) |>
    dplyr::ungroup()

  ggplot2::ggplot(d, ggplot2::aes(pair, axis, size = driver_pctile, colour = driver_b)) +
    ggplot2::geom_point() +
    ggplot2::scale_colour_gradient2(low = "#2166AC", mid = "grey85", high = "#B2182B",
      midpoint = 0, name = "driver slope\n(blue down / red up)") +
    ggplot2::scale_size_continuous(range = c(2, 8), name = "driver prominence\n(percentile)") +
    ggplot2::labs(x = "Focal ← neighbour pair", y = "Ligand → receptor axis", title = title) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1),
                   plot.title  = ggplot2::element_text(hjust = 0.5),
                   panel.grid.minor = ggplot2::element_blank())
}
