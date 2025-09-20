suppressPackageStartupMessages({
  library(Seurat)
  library(GSEABase)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(purrr)
  library(tibble)
})

setwd("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/3geneSetScoring")

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
  y_lab <- paste0(gs$name, " module score")
  p <- ggplot(df, aes(x = Condition, y = Score, fill = Condition)) +
    geom_violin(scale = "width", trim = TRUE, adjust = 1, width = 0.8, color = NA) +
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
    labs(y = y_lab, title = paste0(cell_type, " — ", gs$name)) +
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

res <- score_and_compare(
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

print(res$plot)
res$p.value
res$n_genes_used
res$signature

# Basal, wnt4
target_celltype <- "Basal"
gmt_file <- "HALLMARK_WNT_BETA_CATENIN_SIGNALING.v2025.1.Mm.gmt"
basal_wnt4 <- score_and_compare(
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

print(basal_wnt4$plot)

# stroma, inf
target_celltype <- "Stroma"
gmt_file <- "HALLMARK_INFLAMMATORY_RESPONSE.v2025.1.Mm.gmt"
stroma_inf <- score_and_compare(
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

print(stroma_inf$plot)

# HS RANKL
target_celltype <- "HormSens"
gmt_file <- "HALLMARK_TNFA_SIGNALING_VIA_NFKB.v2025.1.Mm.gmt"
hormsens_rankl <- score_and_compare(
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
print(hormsens_rankl$plot)

# Notch
target_celltype <- "Basal"
gmt_file <- "HALLMARK_NOTCH_SIGNALING.v2025.1.Mm.gmt"
basal_notch <- score_and_compare(
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
print(basal_notch$plot)

target_celltype <- "Immune"
gmt_file <- "HALLMARK_NOTCH_SIGNALING.v2025.1.Mm.gmt"
hormsens_rankl <- score_and_compare(
  seu,
  cell_type        = target_celltype,
  gmt_file         = gmt_file,
  gene_set_name    = gene_set_name,
  condition_col    = condition_col,
  celltype_col     = "cell_type",
  assay            = assay_to_use,
  slot             = slot_to_use,
  condition_levels = condition_levels,
  condition_colors = condition_colors
)
print(hormsens_rankl$plot)
