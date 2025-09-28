suppressPackageStartupMessages({
  library(Seurat)
  library(GSEABase)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(purrr)
  library(tibble)
})

setwd("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/fig3-stroma/scoreEpiCellLeadingEdge/geneSets")

seu <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_All_sub.rds")
gene_set_name <- NULL

## Optional: batch list example (uncomment to use later)
# gmt_list <- c(
#   "GOBP_RESPONSE_TO_INTERLEUKIN_4.v2025.1.Mm.gmt",
#   "GOMF_CHEMOKINE_RECEPTOR_BINDING.v2025.1.Mm.gmt",
#   "HALLMARK_HEDGEHOG_SIGNALING.v2025.1.Mm.gmt",
#   "HALLMARK_INFLAMMATORY_RESPONSE.v2025.1.Mm.gmt",
#   "HALLMARK_NOTCH_SIGNALING.v2025.1.Mm.gmt",
#   "HALLMARK_TNFA_SIGNALING_VIA_NFKB.v2025.1.Mm.gmt",
#   "HALLMARK_WNT_BETA_CATENIN_SIGNALING.v2025.1.Mm.gmt",
#   "REACTOME_SIGNALING_BY_ERBB4.v2025.1.Mm.gmt",
#   "REACTOME_SIGNALING_BY_FGFR.v2025.1.Mm.gmt"
# )

condition_col <- "orig.ident"
celltype_col  <- "subcluster"
assay_to_use  <- "RNA"
slot_to_use   <- "data"
condition_levels <- c("ND", "HFD")

condition_colors <- c("ND" = "#74c5be", "HFD" = "#e95503")

load_gmt_one_set <- function(path, gene_set_name = NULL) {
  gsc <- GSEABase::getGmt(path)
  gs_names <- vapply(gsc, GSEABase::setName, character(1))

  if (length(gsc) == 1L && is.null(gene_set_name)) {
    return(list(name = gs_names[1], genes = unique(GSEABase::geneIds(gsc[[1]]))))
  }
  if (is.null(gene_set_name)) {
    stop(paste0(
      "The GMT contains multiple gene sets. Please set `gene_set_name`.\nOptions: ",
      paste(gs_names, collapse = ", ")
    ))
  }
  idx <- which(gs_names == gene_set_name)
  if (length(idx) == 0) {
    stop(paste0(
      "Gene set not found in GMT: ", gene_set_name,
      "\nOptions: ", paste(gs_names, collapse = ", ")
    ))
  }
  list(name = gs_names[idx[1]], genes = unique(GSEABase::geneIds(gsc[[idx[1]]])))
}

score_and_compare <- function(seu,
                              cell_type,
                              gmt_file,
                              gene_set_name = NULL,
                              condition_col = "orig.ident",
                              celltype_col  = "subcluster",
                              assay = "RNA",
                              slot  = "data",
                              min_genes = 5,
                              condition_levels = c("ND", "HFD"),
                              condition_colors = c("ND"="#74c5be","HFD"="#e95503")) {
  DefaultAssay(seu) <- assay

  # subset by cell type
  cells_keep <- rownames(seu@meta.data)[seu@meta.data[[celltype_col]] == cell_type]
  if (length(cells_keep) == 0) stop("No cells found for cell type in column '", celltype_col, "': ", cell_type)
  obj <- subset(seu, cells = cells_keep)

  # load gene set and intersect with available genes
  gs <- load_gmt_one_set(gmt_file, gene_set_name)
  features_use <- intersect(rownames(obj), gs$genes)
  if (length(features_use) < min_genes) {
    warning("Only ", length(features_use), " signature genes found in the object (min_genes=", min_genes, "). ",
            "Check gene symbol case or species mapping.")
  }

  # scoring (Seurat AddModuleScore)
  obj <- Seurat::AddModuleScore(
    object   = obj,
    features = list(features_use),
    name     = "sig",
    assay    = assay,
    slot     = slot
  )
  score_col <- grep("^sig", colnames(obj@meta.data), value = TRUE)[1]

  # assemble data frame for stats and plotting
  df <- obj@meta.data %>%
    transmute(
      Condition = factor(.data[[condition_col]], levels = condition_levels),
      Score = .data[[score_col]]
    ) %>%
    filter(!is.na(Condition))

  if (nlevels(df$Condition) < 2 || any(table(df$Condition) == 0)) {
    stop("Within this cell type, at least one of ND/HFD has zero cells.")
  }

  # two-sided Wilcoxon rank-sum test
  wt <- wilcox.test(Score ~ Condition, data = df, exact = FALSE)

  # plot (your style)
  y_lab <- paste0(" module score")
  p <- ggplot(df, aes(x = Condition, y = Score, fill = Condition)) +
    geom_violin(scale = "width", trim = F, adjust = 1, width = 0.8, color = NA) +
    geom_boxplot(width = 0.12, fill = "white", alpha = 1,
                 outlier.size = 0.8, outlier.color = "black",
                 outlier.shape = 19, linewidth = 0.3) +
    scale_fill_manual(values = condition_colors) +
    theme_classic(base_family = "Arial") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 12),
          axis.text.x  = element_text(size = 10),
          axis.text.y  = element_text(size = 10),
          legend.position = "none") +
    labs(y = y_lab, title = "") +
    ggpubr::stat_compare_means(comparisons = list(c("ND","HFD")),
                               method = "wilcox.test", label = "p.format")

  list(
    data = df,
    plot = p,
    p.value = wt$p.value,
    signature = gs$name,
    n_genes_used = length(features_use),
    cell_type = cell_type
  )
}

# res <- score_and_compare(
#   seu,
#   cell_type        = celltype,
#   gmt_file         = gmt_file,
#   gene_set_name    = gene_set_name,
#   condition_col    = condition_col,
#   celltype_col     = celltype_col,
#   assay            = assay_to_use,
#   slot             = slot_to_use,
#   condition_levels = condition_levels,
#   condition_colors = condition_colors
# )
# 
# print(res$plot)
# res$p.value
# res$n_genes_used
# res$signature

# Basal, epi migration
target_celltype <- "Basal"
gmt_file <- "GOBP_POSITIVE_REGULATION_OF_EPITHELIAL_CELL_MIGRATION.v2025.1.Mm.gmt"
basal_epi_migration <- score_and_compare(
  seu,
  cell_type        = target_celltype,
  gmt_file         = gmt_file,
  gene_set_name    = gene_set_name,
  condition_col    = condition_col,
  celltype_col     = celltype_col,
  assay            = assay_to_use,
  slot             = slot_to_use,
  condition_levels = condition_levels,
  condition_colors = condition_colors
)

print(basal_epi_migration$plot)

# HormSens, epi migration
target_celltype <- "HormSens"
gmt_file <- "GOBP_POSITIVE_REGULATION_OF_EPITHELIAL_CELL_MIGRATION.v2025.1.Mm.gmt"
hormsens_epi_migration <- score_and_compare(
  seu,
  cell_type        = target_celltype,
  gmt_file         = gmt_file,
  gene_set_name    = gene_set_name,
  condition_col    = condition_col,
  celltype_col     = celltype_col,
  assay            = assay_to_use,
  slot             = slot_to_use,
  condition_levels = condition_levels,
  condition_colors = condition_colors
)

print(hormsens_epi_migration$plot)

# LumProg, epi migration
target_celltype <- "LumProg"
gmt_file <- "GOBP_POSITIVE_REGULATION_OF_EPITHELIAL_CELL_MIGRATION.v2025.1.Mm.gmt"
lumprog_epi_migration <- score_and_compare(
  seu,
  cell_type        = target_celltype,
  gmt_file         = gmt_file,
  gene_set_name    = gene_set_name,
  condition_col    = condition_col,
  celltype_col     = celltype_col,
  assay            = assay_to_use,
  slot             = slot_to_use,
  condition_levels = condition_levels,
  condition_colors = condition_colors
)
print(lumprog_epi_migration$plot)

# Basal, migration
tararget_celltype <- "Basal"
gmt_file <- "GOBP_CELL_MIGRATION.v2025.1.Mm.gmt"
basal_migration <- score_and_compare(
  seu,
  cell_type        = target_celltype,
  gmt_file         = gmt_file,
  gene_set_name    = gene_set_name,
  condition_col    = condition_col,
  celltype_col     = celltype_col,
  assay            = assay_to_use,
  slot             = slot_to_use,
  condition_levels = condition_levels,
  condition_colors = condition_colors
)

print(basal_migration$plot)

# HormSens, migration
target_celltype <- "HormSens"
gmt_file <- "GOBP_CELL_MIGRATION.v2025.1.Mm.gmt"
hormsens_migration <- score_and_compare(
  seu,
  cell_type        = target_celltype,
  gmt_file         = gmt_file,
  gene_set_name    = gene_set_name,
  condition_col    = condition_col,
  celltype_col     = celltype_col,
  assay            = assay_to_use,
  slot             = slot_to_use,
  condition_levels = condition_levels,
  condition_colors = condition_colors
)
print(hormsens_migration$plot)

# LumProg, migration
target_celltype <- "LumProg"
gmt_file <- "GOBP_CELL_MIGRATION.v2025.1.Mm.gmt"
lumprog_migration <- score_and_compare(
  seu,
  cell_type        = target_celltype,
  gmt_file         = gmt_file,
  gene_set_name    = gene_set_name,
  condition_col    = condition_col,
  celltype_col     = celltype_col,
  assay            = assay_to_use,
  slot             = slot_to_use,
  condition_levels = condition_levels,
  condition_colors = condition_colors
)
print(lumprog_migration$plot)

# Basal, cell leading edge
target_celltype <- "Basal"
gmt_file <- "GOCC_CELL_LEADING_EDGE.v2025.1.Mm.gmt"
basal_cell_leading_edge <- score_and_compare(
  seu,
  cell_type        = target_celltype,
  gmt_file         = gmt_file,
  gene_set_name    = gene_set_name,
  condition_col    = condition_col,
  celltype_col     = celltype_col,
  assay            = assay_to_use,
  slot             = slot_to_use,
  condition_levels = condition_levels,
  condition_colors = condition_colors
)
print(basal_cell_leading_edge$plot)

# HormSens, cell leading edge
target_celltype <- "HormSens"
gmt_file <- "GOCC_CELL_LEADING_EDGE.v2025.1.Mm.gmt"
hormsens_cell_leading_edge <- score_and_compare(
  seu,
  cell_type        = target_celltype,
  gmt_file         = gmt_file,
  gene_set_name    = gene_set_name,
  condition_col    = condition_col,
  celltype_col     = celltype_col,
  assay            = assay_to_use,
  slot             = slot_to_use,
  condition_levels = condition_levels,
  condition_colors = condition_colors
)

print(hormsens_cell_leading_edge$plot)

# LumProg, cell leading edge
target_celltype <- "LumProg"
gmt_file <- "GOCC_CELL_LEADING_EDGE.v2025.1.Mm.gmt"
lumprog_cell_leading_edge <- score_and_compare(
  seu,
  cell_type        = target_celltype,
  gmt_file         = gmt_file,
  gene_set_name    = gene_set_name,
  condition_col    = condition_col,
  celltype_col     = celltype_col,
  assay            = assay_to_use,
  slot             = slot_to_use,
  condition_levels = condition_levels,
  condition_colors = condition_colors
)
print(lumprog_cell_leading_edge$plot)

# patch together
library(patchwork)
final_figure <- (basal_epi_migration$plot + labs(title = "Basal")) +
  (hormsens_epi_migration$plot + labs(title = "HormSens")) +
  (lumprog_epi_migration$plot + labs(title = "LumProg")) +
  (basal_migration$plot + labs(title = "Basal")) +
  (hormsens_migration$plot + labs(title = "HormSens")) +
  (lumprog_migration$plot + labs(title = "LumProg")) +
  (basal_cell_leading_edge$plot + labs(title = "Basal")) +
  (hormsens_cell_leading_edge$plot + labs(title = "HormSens")) +
  (lumprog_cell_leading_edge$plot + labs(title = "LumProg")) +
  plot_layout(ncol = 3, guides = "collect") &
  theme(legend.position = 'bottom')
final_figure

##################################
######### dont run this ##########
##################################

suppressPackageStartupMessages({
  library(tidyr)     # for pivot_wider
  library(stringr)
})

# Small helper to format FDR labels
.format_fdr <- function(padj) {
  dplyr::case_when(
    padj < 1e-3 ~ "FDR<0.001",
    TRUE        ~ paste0("FDR=", signif(padj, 2))
  )
}

# One-row-per-(cell_type × gene set) result; dumbbell plot of ND vs HFD summaries
integrated_dumbbell <- function(seu,
                                plan_tbl,            # tibble with columns: cell_type, gmt_file, gs_name (or NA)
                                summary_stat = c("median","mean"),
                                condition_col = "orig.ident",
                                celltype_col  = "subcluster",
                                assay         = "RNA",
                                slot          = "data",
                                condition_levels = c("ND","HFD"),
                                condition_colors = c("ND"="#74c5be","HFD"="#e95503")) {

  summary_stat <- match.arg(summary_stat)

  safe_score <- purrr::possibly(score_and_compare, otherwise = NULL)

  results <- purrr::pmap(
    plan_tbl,
    function(cell_type, gmt_file, gs_name) {
      safe_score(
        seu,
        cell_type        = cell_type,
        gmt_file         = gmt_file,
        gene_set_name    = if (is.na(gs_name) || is.null(gs_name)) NULL else gs_name,
        condition_col    = condition_col,
        celltype_col     = celltype_col,
        assay            = assay,
        slot             = slot,
        condition_levels = condition_levels,
        condition_colors = condition_colors
      )
    }
  ) |> purrr::compact()

  if (length(results) == 0) stop("No successful results to plot.")

  # Build a tidy summary table (medians/means per condition + p, FDR, etc.)
  sum_tbl <- purrr::map_dfr(results, function(res) {
    df <- res$data
    gsum <- df |>
      dplyr::group_by(Condition) |>
      dplyr::summarise(
        median = median(Score),
        mean   = mean(Score),
        n      = dplyr::n(),
        .groups = "drop"
      ) |>
      tidyr::pivot_wider(names_from = Condition,
                         values_from = c(median, mean, n),
                         names_sep = ".")

    dplyr::tibble(
      cell_type = res$cell_type,
      gene_set  = res$signature,
      p         = res$p.value,
      n_genes   = res$n_genes_used
    ) |>
      dplyr::bind_cols(gsum)
  })

  # Choose your summary (median or mean) for the dumbbell
  a_col <- paste0(summary_stat, ".", condition_levels[1])  # e.g., "median.ND"
  b_col <- paste0(summary_stat, ".", condition_levels[2])  # e.g., "median.HFD"

  # BH adjust, effect size (difference), pretty labels, ordering
  sum_tbl <- sum_tbl |>
    dplyr::mutate(
      padj      = p.adjust(p, method = "BH"),
      diff_val  = .data[[b_col]] - .data[[a_col]],
      row_lab   = paste0(cell_type, " \u2022 ", gene_set),     # “Basal • GO TERM”
      star      = dplyr::case_when(padj < 0.001 ~ "***",
                                   padj < 0.01  ~ "**",
                                   padj < 0.05  ~ "*",
                                   TRUE         ~ ""),
      p_text    = paste0(star, " ", .format_fdr(padj))
    ) |>
    dplyr::arrange(diff_val) |>
    dplyr::mutate(row_lab = factor(row_lab, levels = row_lab))

  # Space on the right for p/FDR labels
  x_min <- min(c(sum_tbl[[a_col]], sum_tbl[[b_col]]), na.rm = TRUE)
  x_max <- max(c(sum_tbl[[a_col]], sum_tbl[[b_col]]), na.rm = TRUE)
  pad   <- 0.08 * (x_max - x_min + 1e-8)
  sum_tbl$p_x <- x_max + pad

  # Build the dumbbell/forest hybrid
  p <- ggplot2::ggplot(sum_tbl) +
    ggplot2::geom_segment(
      ggplot2::aes(x = .data[[a_col]], xend = .data[[b_col]], y = row_lab, yend = row_lab),
      linewidth = 0.3, color = "#b8b8b8"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = .data[[a_col]], y = row_lab),
      size = 2.3, color = condition_colors[[condition_levels[1]]]
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = .data[[b_col]], y = row_lab),
      size = 2.3, color = condition_colors[[condition_levels[2]]]
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = p_x, y = row_lab, label = p_text),
      hjust = 0, size = 3.1
    ) +
    ggplot2::theme_classic(base_family = "Arial") +
    ggplot2::theme(
      axis.text.y  = ggplot2::element_text(size = 9),
      plot.title   = ggplot2::element_text(size = 12, face = "bold"),
      plot.subtitle= ggplot2::element_text(size = 10)
    ) +
    ggplot2::coord_cartesian(xlim = c(x_min, x_max + pad * 3), clip = "off")

  list(summary_table = sum_tbl, plot = p)
}

plan_tbl <- tibble::tribble(
  ~cell_type, ~gmt_file,                                                     ~gs_name,
  "Basal",   "GOBP_POSITIVE_REGULATION_OF_EPITHELIAL_CELL_MIGRATION.v2025.1.Mm.gmt", NA,
  "HormSens","GOBP_POSITIVE_REGULATION_OF_EPITHELIAL_CELL_MIGRATION.v2025.1.Mm.gmt", NA,
  "LumProg", "GOBP_POSITIVE_REGULATION_OF_EPITHELIAL_CELL_MIGRATION.v2025.1.Mm.gmt", NA,
  "Basal",   "GOBP_CELL_MIGRATION.v2025.1.Mm.gmt",                                 NA,
  "HormSens","GOBP_CELL_MIGRATION.v2025.1.Mm.gmt",                                 NA,
  "LumProg", "GOBP_CELL_MIGRATION.v2025.1.Mm.gmt",                                 NA,
  "Basal",   "GOCC_CELL_LEADING_EDGE.v2025.1.Mm.gmt",                              NA,
  "HormSens","GOCC_CELL_LEADING_EDGE.v2025.1.Mm.gmt",                              NA,
  "LumProg", "GOCC_CELL_LEADING_EDGE.v2025.1.Mm.gmt",                              NA
)

suppressPackageStartupMessages({
  library(tidyr)
  library(stringr)
})

.format_fdr <- function(padj) {
  dplyr::case_when(
    padj < 1e-3 ~ "FDR<0.001",
    TRUE        ~ paste0("FDR=", signif(padj, 2))
  )
}

##################################
########### run this #############
##################################

integrated_dumbbell <- function(seu,
                                plan_tbl,
                                summary_stat = c("median", "mean"),
                                condition_col = "orig.ident",
                                celltype_col  = "subcluster",
                                assay         = "RNA",
                                slot          = "data",
                                condition_levels = c("ND","HFD"),
                                condition_colors = c("ND"="#74c5be",
                                                     "HFD"="#e95503"),
                                min_genes = 5) {

  summary_stat <- match.arg(summary_stat)
  DefaultAssay(seu) <- assay

  uniq_sets <- plan_tbl |>
    dplyr::distinct(gmt_file, gs_name)

  gs_meta <- vector("list", nrow(uniq_sets))
  obj <- seu
  for (i in seq_len(nrow(uniq_sets))) {
    gmt_file <- uniq_sets$gmt_file[i]
    gs_name  <- uniq_sets$gs_name[i]

    gs <- load_gmt_one_set(gmt_file, if (is.na(gs_name)) NULL else gs_name)
    feats <- intersect(rownames(obj), gs$genes)
    if (length(feats) < min_genes) {
      warning("Only ", length(feats), " signature genes found (min_genes=", min_genes, ") for ", gs$name)
    }

    before <- colnames(obj@meta.data)
    obj <- Seurat::AddModuleScore(
      object   = obj,
      features = list(feats),
      name     = paste0("igs", i),
      assay    = assay,
      slot     = slot
    )
    after <- colnames(obj@meta.data)
    newcol <- setdiff(after, before)
    
    score_col <- newcol[grep(paste0("^igs", i), newcol)][1]

    gs_meta[[i]] <- list(
      gmt_file  = gmt_file,
      gs_name   = gs$name,
      score_col = score_col,
      n_genes   = length(feats)
    )
  }
  gs_meta <- dplyr::bind_rows(gs_meta)

  res_rows <- purrr::pmap(
    plan_tbl,
    function(cell_type, gmt_file, gs_name) {
      meta_row <- gs_meta |>
        dplyr::filter(.data$gmt_file == !!gmt_file) |>
        dplyr::slice(1)

      df <- obj@meta.data |>
        dplyr::transmute(
          cell_type  = .data[[celltype_col]],
          Condition  = factor(.data[[condition_col]], levels = condition_levels),
          Score      = .data[[meta_row$score_col]]
        ) |>
        dplyr::filter(.data$cell_type == !!cell_type, !is.na(Condition))

      if (nrow(df) == 0 || nlevels(df$Condition) < 2 || any(table(df$Condition) == 0)) {
        return(NULL)
      }

      wt <- wilcox.test(Score ~ Condition, data = df, exact = FALSE)

      gsum <- df |>
        dplyr::group_by(Condition) |>
        dplyr::summarise(
          median = median(Score),
          mean   = mean(Score),
          n      = dplyr::n(),
          .groups = "drop"
        ) |>
        tidyr::pivot_wider(names_from = Condition,
                           values_from = c(median, mean, n),
                           names_sep = ".")

      dplyr::tibble(
        cell_type = cell_type,
        gene_set  = meta_row$gs_name,
        p         = wt$p.value,
        n_genes   = meta_row$n_genes
      ) |>
        dplyr::bind_cols(gsum)
    }
  ) |> purrr::compact() |> dplyr::bind_rows()

  if (nrow(res_rows) == 0) stop("No successful results to plot.")

  a_col <- paste0(summary_stat, ".", condition_levels[1])  # ND
  b_col <- paste0(summary_stat, ".", condition_levels[2])  # HFD

  sum_tbl <- res_rows |>
    dplyr::mutate(
      padj     = p.adjust(p, method = "BH"),
      diff_val = .data[[b_col]] - .data[[a_col]],
      star     = dplyr::case_when(padj < 0.001 ~ "***",
                                  padj < 0.01  ~ "**",
                                  padj < 0.05  ~ "*",
                                  TRUE         ~ ""),
      p_text   = paste0(star, " ", .format_fdr(padj))
    )

  gene_order <- plan_tbl |>
    dplyr::mutate(gs = NULL) |>
    dplyr::distinct(gmt_file) |>
    dplyr::left_join(gs_meta, by = "gmt_file") |>
    dplyr::pull(gs_name) |>
    unique()

  cell_order <- c("Basal", "HormSens", "LumProg")

  sum_tbl <- sum_tbl |>
    dplyr::mutate(
      row_lab = paste0(gene_set, " \u2022 ", cell_type)
    )

  y_levels <- as.character(unlist(lapply(gene_order, function(gs)
    paste0(gs, " \u2022 ", cell_order))))
  y_levels <- y_levels[y_levels %in% sum_tbl$row_lab]

  sum_tbl <- sum_tbl |>
    dplyr::mutate(row_lab = factor(row_lab, levels = rev(y_levels)))

  x_min <- min(c(sum_tbl[[a_col]], sum_tbl[[b_col]]), na.rm = TRUE)
  x_max <- max(c(sum_tbl[[a_col]], sum_tbl[[b_col]]), na.rm = TRUE)
  pad   <- 0.10 * (x_max - x_min + 1e-8)
  sum_tbl$p_x <- x_max + pad

  p <- ggplot2::ggplot(sum_tbl) +
    ggplot2::geom_segment(
      ggplot2::aes(x = .data[[a_col]], xend = .data[[b_col]], y = row_lab, yend = row_lab),
      linewidth = 0.3, color = "#b8b8b8"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = .data[[a_col]], y = row_lab),
      size = 2, color = condition_colors[[condition_levels[1]]]
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = .data[[b_col]], y = row_lab),
      size = 2, color = condition_colors[[condition_levels[2]]]
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = p_x, y = row_lab, label = p_text),
      hjust = 0, size = 3.1
    ) +
    ggplot2::labs(
      x = paste0("Module score (", summary_stat, ")"),
      y = NULL
    ) +
    ggplot2::theme_classic(base_family = "Arial") +
    ggplot2::theme(
      axis.text.y   = ggplot2::element_text(size = 8, color = "black"),
      axis.text.x   = ggplot2::element_text(size = 8, color = "black"),
      axis.line.y = ggplot2::element_line(size = 0.3),
      axis.line.x = ggplot2::element_line(size = 0.3),
      axis.ticks.y  = ggplot2::element_line(size = 0.3),
      axis.ticks.x  = ggplot2::element_line(size = 0.3),
      plot.margin   = grid::unit(c(5, 28, 5, 5), "pt")
    ) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.40))) +
    ggplot2::coord_cartesian(clip = "off")

  list(summary_table = sum_tbl, plot = p)
}

res <- integrated_dumbbell(
  seu,
  plan_tbl,
  summary_stat     = "median",     # or "mean"
  condition_col    = condition_col,
  celltype_col     = celltype_col,
  assay            = assay_to_use,
  slot             = slot_to_use,
  condition_levels = condition_levels,
  condition_colors = condition_colors
)

print(res$plot)

res$summary_table