# =========== Prep ===========
obj_path <- "/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_All_sub.rds"
gs_dir   <- "/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/fig3-stroma/epiDNAAUCellScoring/geneSets"
out_dir  <- "/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/fig3-stroma/epiDNAAUCellScoring"

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(GSVA)
  library(GSEABase)
  library(ggsci)
})

# =========== Load Seurat object ===========
obj <- readRDS(obj_path)
subcol <- "cell_type"
obj$combo <- paste0(obj$orig.ident, "_", obj@meta.data[[subcol]])

# =========== Expression matrix for GSVA ===========
DefaultAssay(obj) <- "RNA"

Idents(obj) <- "cell_type"
obj <- subset(obj, idents = c("Basal","LumProg","HormSens"))

expr <- as.matrix(GetAssayData(obj, slot = "data"))

dup <- duplicated(rownames(expr))
if (any(dup)) expr <- expr[!dup, , drop = FALSE]
rownames(expr) <- toupper(rownames(expr))

# =========== Read gene sets (GMT) ===========
gmt_files <- list.files(gs_dir, pattern = "\\.gmt$", full.names = TRUE)

read_gmt_as_list <- function(gmt_path) {
  gsc <- getGmt(gmt_path)
  lst <- lapply(gsc, function(gs) toupper(geneIds(gs)))
  names(lst) <- paste0(basename(gmt_path), "::", sapply(gsc, setName))
  lst
}
gs_list <- unlist(lapply(gmt_files, read_gmt_as_list), recursive = FALSE)

# min_sz <- 10
# max_sz <- 500
# gs_list <- lapply(gs_list, function(v) intersect(v, rownames(expr)))
# gs_list <- gs_list[sapply(gs_list, function(v) length(v) >= min_sz & length(v) <= max_sz)]
# if (length(gs_list) == 0) stop("Stop")

# # =========== Run GSVA ===========
# set.seed(1234)
# gp <- GSVA::gsvaParam(
#   exprData = expr,
#   geneSets = gs_list,
#   kcdf     = "Gaussian",
#   minSize  = min_sz,
#   maxSize  = max_sz,
#   maxDiff  = TRUE
# )
# 
# gsva_scores <- GSVA::gsva(gp)
# 
# saveRDS(gsva_scores, file = file.path(out_dir, "GSVA_scores_by_cell.rds"))
# 
# # =========== Long format + stats ===========
# scores_df <- as.data.frame(t(gsva_scores))
# scores_df$cell  <- rownames(scores_df)
# scores_df$cell_type <- obj$cell_type[rownames(scores_df)]
# 
# long_df <- scores_df |>
#   tidyr::pivot_longer(cols = -c(cell, cell_type),
#                       names_to = "gene_set",
#                       values_to = "score") |>
#   dplyr::filter(!is.na(cell_type))
# 
# long_df_norm <- long_df |>
#   group_by(gene_set) |>
#   mutate(
#     mu    = mean(score, na.rm = TRUE),
#     sigma = sd(score,  na.rm = TRUE),
#     score_z = ifelse(is.finite(sigma) & sigma > 0, (score - mu) / sigma, 0),
#     smin  = min(score, na.rm = TRUE),
#     smax  = max(score, na.rm = TRUE),
#     score_minmax = ifelse(smax > smin, (score - smin) / (smax - smin), 0)
#   ) |>
#   ungroup() |>
#   dplyr::select(-mu, -sigma, -smin, -smax)
# 
# summary_df <- long_df_norm |>
#   group_by(gene_set, cell_type) |>
#   summarise(
#     n = dplyr::n(),
#     mean_raw   = mean(score,        na.rm = TRUE),
#     median_raw = median(score,      na.rm = TRUE),
#     mean_z     = mean(score_z,      na.rm = TRUE),
#     median_z   = median(score_z,    na.rm = TRUE),
#     mean_mm    = mean(score_minmax, na.rm = TRUE),
#     median_mm  = median(score_minmax, na.rm = TRUE),
#     .groups = "drop"
#   )
# 
# p_z <- ggplot(long_df_norm, aes(x = cell_type, y = score_z, fill = cell_type)) +
#   geom_boxplot(outlier.size = 0.4, width = 0.7) +
#   facet_wrap(~ gene_set, scales = "free_y", ncol = 2) +
#   scale_fill_npg(guide = "none") +
#   labs(title = "GSVA (z-score) by cell_type", x = "cell_type", y = "GSVA z-score") +
#   theme_bw(base_size = 11) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1),
#         strip.text = element_text(size = 9))

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr); library(tidyr); library(ggplot2)
  library(GSEABase)
  library(ggsci)
  library(AUCell)
})

# =========== Run AUCell ===========
set.seed(1234)

rankings <- AUCell_buildRankings(expr, 
                                 plotStats = FALSE, 
                                 verbose = TRUE)

min_sz <- 5
max_sz <- 500
gs_list <- lapply(gs_list, function(v) intersect(v, rownames(rankings)))
gs_list <- gs_list[sapply(gs_list, function(v) length(v) >= min_sz & length(v) <= max_sz)]

nGenes <- nrow(rankings)
aucMaxRank <- min(nGenes, max(50, ceiling(0.05 * nGenes)))

cellsAUC <- AUCell_calcAUC(gs_list, rankings, aucMaxRank = aucMaxRank, nCores = 1)
auc_mat  <- as.matrix(getAUC(cellsAUC))

scores_df <- as.data.frame(t(auc_mat))
scores_df$cell  <- rownames(scores_df)
scores_df$combo <- obj$combo[rownames(scores_df)]

long_df <- scores_df |>
  tidyr::pivot_longer(cols = -c(cell, combo),
                      names_to = "gene_set",
                      values_to = "score") |>
  dplyr::filter(!is.na(combo))

long_df_norm <- long_df |>
  group_by(gene_set) |>
  mutate(
    mu    = mean(score, na.rm = TRUE),
    sigma = sd(score,  na.rm = TRUE),
    score_z = ifelse(is.finite(sigma) & sigma > 0, (score - mu) / sigma, 0),
    smin  = min(score, na.rm = TRUE),
    smax  = max(score, na.rm = TRUE),
    score_minmax = ifelse(smax > smin, (score - smin) / (smax - smin), 0)
  ) |>
  ungroup() |>
  dplyr::select(-mu, -sigma, -smin, -smax)

summary_df <- long_df_norm |>
  group_by(gene_set, combo) |>
  summarise(
    n = dplyr::n(),
    mean_raw   = mean(score,        na.rm = TRUE),
    median_raw = median(score,      na.rm = TRUE),
    mean_z     = mean(score_z,      na.rm = TRUE),
    median_z   = median(score_z,    na.rm = TRUE),
    mean_mm    = mean(score_minmax, na.rm = TRUE),
    median_mm  = median(score_minmax, na.rm = TRUE),
    .groups = "drop"
  )

p_z <- 
  ggplot(
    long_df_norm, 
         aes(x = combo, y = score_z, fill = combo)
         ) +
  geom_boxplot(outlier.size = 0.35, 
               width = 0.65,
               staplewidth = 1,
               color = "black",
               size = 0.3) +
  # stat_summary(
  #   fun = median,
  #   geom = "crossbar",
  #   width = 0.65,
  #   fatten = 0,
  #   color = "white",
  #   size = 0.1
  # ) +
  facet_wrap(~ gene_set, scales = "free_y", ncol = 2) +
  scale_fill_npg(guide = "none") +
  labs(title = "AUCell AUC (z-score) by combo", x = "Combo", y = "AUC z-score") +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
        axis.text.y = element_text(colour = "black"),
        strip.text = element_text(size = 9),
        axis.line.y = element_line(color = "black", size = 0.3),
        axis.line.x = element_line(color = "black", size = 0.3),
        axis.ticks = element_line(color = "black", size = 0.3))
print(p_z)

make_gs_labels <- function(x) {
  x1 <- sub("^.*::", "", x)

  x2 <- gsub("_", " ", x1)

  x3 <- sub("^\\S+\\s*", "", x2)
  x3 <- trimws(x3)

  x3[x3 == ""] <- x2[x3 == ""]
  x3
}

gs_levels <- unique(long_df_norm$gene_set)
gs_labmap <- setNames(make_gs_labels(gs_levels), gs_levels)

p_z <- 
  ggplot(long_df_norm, aes(x = combo, y = score_z, fill = combo)) +
  geom_boxplot(
    outlier.size = 0.2, width = 0.65, staplewidth = 1,
    color = "black", size = 0.25
  ) +
  facet_wrap(
    ~ gene_set, scales = "free_y", ncol = 2,
    labeller = labeller(gene_set = gs_labmap)
  ) +
  scale_fill_npg(guide = "none") +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
    axis.text.y = element_text(colour = "black"),
    strip.text  = element_text(size = 5),
    axis.line.y = element_line(color = "black", size = 0.3),
    axis.line.x = element_line(color = "black", size = 0.3),
    axis.ticks  = element_line(color = "black", size = 0.3)
  )

print(p_z) # pdf: 2.74 x 7.33
