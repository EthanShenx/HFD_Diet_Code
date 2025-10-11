setwd("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/fig6-adipo")

library(Seurat)
library(SingleCellExperiment)
library(slingshot)
library(tidyverse)
library(scater)
library(RColorBrewer)
library(grDevices)

set.seed(12345)
options(stringsAsFactors = FALSE)

Adipo_All <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_All_Adipo_sub.rds")
Adipo_ND <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_ND_Adipo_sub.rds")
Adipo_HFD <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_HFD_Adipo_sub.rds")

outdir <- "slingshot_results"
dir.create(outdir, showWarnings = FALSE)

use_cluster_col <- "cell_type"
start_cluster   <- NA
approx_points   <- 150
npcs_for_sling  <- 20

## ==================================================================================
## 1) ALL
## ==================================================================================
seurat_obj <- Adipo_All
obj_tag    <- "All"

## cluster labels (prefer 'seurat_clusters', otherwise Idents())
if (!is.null(use_cluster_col) && use_cluster_col %in% colnames(seurat_obj@meta.data)) {
  cl <- seurat_obj@meta.data[[use_cluster_col]]
} else {
  cl <- as.character(Idents(seurat_obj))
}
cl <- as.factor(cl)
seurat_obj$slingshot_cluster <- cl

## choose reduction with PCA priority; DR used by Slingshot (can be >2D), DR2D for 2D storage
rd_names <- names(seurat_obj@reductions)
if ("pca" %in% rd_names) {
  pca_mat <- seurat_obj@reductions$pca@cell.embeddings
  coords_for_sling <- pca_mat[, seq_len(min(ncol(pca_mat), npcs_for_sling)), drop = FALSE]
  coords_2d        <- pca_mat[, 1:2, drop = FALSE]
  red_source       <- "pca"
} else if ("umap" %in% rd_names) {
  umap_mat <- seurat_obj@reductions$umap@cell.embeddings
  coords_for_sling <- umap_mat[, 1:2, drop = FALSE]
  coords_2d        <- umap_mat[, 1:2, drop = FALSE]
  red_source       <- "umap"
} else if ("tsne" %in% rd_names) {
  tsne_mat <- seurat_obj@reductions$tsne@cell.embeddings
  coords_for_sling <- tsne_mat[, 1:2, drop = FALSE]
  coords_2d        <- tsne_mat[, 1:2, drop = FALSE]
  red_source       <- "tsne"
} else if ("harmony" %in% rd_names) {
  harm_mat <- seurat_obj@reductions$harmony@cell.embeddings
  coords_for_sling <- harm_mat[, 1:2, drop = FALSE]
  coords_2d        <- harm_mat[, 1:2, drop = FALSE]
  red_source       <- "harmony"
} else {
  stop("No PCA/UMAP/tSNE/Harmony reduction found in Adipo_All.")
}
npcs_used <- ncol(coords_for_sling)

## build SCE and run slingshot
sce <- as.SingleCellExperiment(seurat_obj)
reducedDim(sce, "DR")   <- as.matrix(coords_for_sling)  # space for Slingshot
reducedDim(sce, "DR2D") <- as.matrix(coords_2d)         # 2D coordinates for storage
colData(sce)$cluster    <- cl
start_clus_arg <- if (is.na(start_cluster)) NULL else as.character(start_cluster)

sce <- slingshot(
  sce,
  clusterLabels = "cluster",
  reducedDim    = "DR",
  start.clus    = start_clus_arg,
  approx_points = approx_points
)

## write back pseudotime (per-lineage) + mean to Seurat meta.data
pt_mat <- slingPseudotime(sce, na = TRUE)            # cells x lineages
if (!is.null(pt_mat) && ncol(pt_mat) > 0) {
  common_cells <- rownames(pt_mat)
  for (i in seq_len(ncol(pt_mat))) {
    v <- pt_mat[common_cells, i]
    seurat_obj[[paste0("slingPT_L", i)]] <- v[colnames(seurat_obj)]
  }
  pt_mean <- apply(pt_mat, 1, function(x) mean(x, na.rm = TRUE))
  seurat_obj$slingPT_mean <- pt_mean[colnames(seurat_obj)]
  ## optional: lineage assignment by min pseudotime
  seurat_obj$slingPT_lineage <- vapply(
    1:nrow(pt_mat),
    function(r) {
      vals <- pt_mat[r, ]
      if (all(is.na(vals))) return(NA_character_)
      paste0("L", which.min(vals))
    }, character(1)
  )[colnames(seurat_obj)]
} else {
  seurat_obj$slingPT_mean    <- NA_real_
  seurat_obj$slingPT_lineage <- NA_character_
}

## store extra Slingshot results into @misc (lightweight list)
seurat_obj@misc$slingshot <- list(
  object_tag       = obj_tag,
  reduction_source = red_source,
  npcs_used        = npcs_used,
  DR               = reducedDim(sce, "DR"),
  DR2D             = reducedDim(sce, "DR2D"),
  cluster_labels   = cl,
  start_cluster    = start_clus_arg,
  pseudotime       = pt_mat,
  curve_weights    = slingCurveWeights(sce),
  lineages         = slingLineages(sce),
  mst              = slingMST(sce),
  curves           = lapply(slingCurves(sce), function(crv) {
    list(s = crv$s, ord = crv$ord, w = crv$w, lambda = crv$lambda)
  }),
  sds              = SlingshotDataSet(sce)  # full SlingshotDataSet for later plotting if needed
)

## write back to original variable
Adipo_All <- seurat_obj

## ==================================================================================
## 2) ND
## ==================================================================================
seurat_obj <- Adipo_ND
obj_tag    <- "ND"

if (!is.null(use_cluster_col) && use_cluster_col %in% colnames(seurat_obj@meta.data)) {
  cl <- seurat_obj@meta.data[[use_cluster_col]]
} else {
  cl <- as.character(Idents(seurat_obj))
}
cl <- as.factor(cl)
seurat_obj$slingshot_cluster <- cl

rd_names <- names(seurat_obj@reductions)
if ("pca" %in% rd_names) {
  pca_mat <- seurat_obj@reductions$pca@cell.embeddings
  coords_for_sling <- pca_mat[, seq_len(min(ncol(pca_mat), npcs_for_sling)), drop = FALSE]
  coords_2d        <- pca_mat[, 1:2, drop = FALSE]
  red_source       <- "pca"
} else if ("umap" %in% rd_names) {
  umap_mat <- seurat_obj@reductions$umap@cell.embeddings
  coords_for_sling <- umap_mat[, 1:2, drop = FALSE]
  coords_2d        <- umap_mat[, 1:2, drop = FALSE]
  red_source       <- "umap"
} else if ("tsne" %in% rd_names) {
  tsne_mat <- seurat_obj@reductions$tsne@cell.embeddings
  coords_for_sling <- tsne_mat[, 1:2, drop = FALSE]
  coords_2d        <- tsne_mat[, 1:2, drop = FALSE]
  red_source       <- "tsne"
} else if ("harmony" %in% rd_names) {
  harm_mat <- seurat_obj@reductions$harmony@cell.embeddings
  coords_for_sling <- harm_mat[, 1:2, drop = FALSE]
  coords_2d        <- harm_mat[, 1:2, drop = FALSE]
  red_source       <- "harmony"
} else {
  stop("No PCA/UMAP/tSNE/Harmony reduction found in Adipo_ND.")
}
npcs_used <- ncol(coords_for_sling)

sce <- as.SingleCellExperiment(seurat_obj)
reducedDim(sce, "DR")   <- as.matrix(coords_for_sling)
reducedDim(sce, "DR2D") <- as.matrix(coords_2d)
colData(sce)$cluster    <- cl
start_clus_arg <- if (is.na(start_cluster)) NULL else as.character(start_cluster)

sce <- slingshot(sce, clusterLabels = "cluster", reducedDim = "DR",
                 start.clus = start_clus_arg, approx_points = approx_points)

pt_mat <- slingPseudotime(sce, na = TRUE)
if (!is.null(pt_mat) && ncol(pt_mat) > 0) {
  common_cells <- rownames(pt_mat)
  for (i in seq_len(ncol(pt_mat))) {
    v <- pt_mat[common_cells, i]
    seurat_obj[[paste0("slingPT_L", i)]] <- v[colnames(seurat_obj)]
  }
  pt_mean <- apply(pt_mat, 1, function(x) mean(x, na.rm = TRUE))
  seurat_obj$slingPT_mean <- pt_mean[colnames(seurat_obj)]
  seurat_obj$slingPT_lineage <- vapply(
    1:nrow(pt_mat),
    function(r) {
      vals <- pt_mat[r, ]
      if (all(is.na(vals))) return(NA_character_)
      paste0("L", which.min(vals))
    }, character(1)
  )[colnames(seurat_obj)]
} else {
  seurat_obj$slingPT_mean    <- NA_real_
  seurat_obj$slingPT_lineage <- NA_character_
}

seurat_obj@misc$slingshot <- list(
  object_tag       = obj_tag,
  reduction_source = red_source,
  npcs_used        = npcs_used,
  DR               = reducedDim(sce, "DR"),
  DR2D             = reducedDim(sce, "DR2D"),
  cluster_labels   = cl,
  start_cluster    = start_clus_arg,
  pseudotime       = pt_mat,
  curve_weights    = slingCurveWeights(sce),
  lineages         = slingLineages(sce),
  mst              = slingMST(sce),
  curves           = lapply(slingCurves(sce), function(crv) {
    list(s = crv$s, ord = crv$ord, w = crv$w, lambda = crv$lambda)
  }),
  sds              = SlingshotDataSet(sce)
)

Adipo_ND <- seurat_obj

## ==================================================================================
## 3) HFD
## ==================================================================================
seurat_obj <- Adipo_HFD
obj_tag    <- "HFD"

if (!is.null(use_cluster_col) && use_cluster_col %in% colnames(seurat_obj@meta.data)) {
  cl <- seurat_obj@meta.data[[use_cluster_col]]
} else {
  cl <- as.character(Idents(seurat_obj))
}
cl <- as.factor(cl)
seurat_obj$slingshot_cluster <- cl

rd_names <- names(seurat_obj@reductions)
if ("pca" %in% rd_names) {
  pca_mat <- seurat_obj@reductions$pca@cell.embeddings
  coords_for_sling <- pca_mat[, seq_len(min(ncol(pca_mat), npcs_for_sling)), drop = FALSE]
  coords_2d        <- pca_mat[, 1:2, drop = FALSE]
  red_source       <- "pca"
} else if ("umap" %in% rd_names) {
  umap_mat <- seurat_obj@reductions$umap@cell.embeddings
  coords_for_sling <- umap_mat[, 1:2, drop = FALSE]
  coords_2d        <- umap_mat[, 1:2, drop = FALSE]
  red_source       <- "umap"
} else if ("tsne" %in% rd_names) {
  tsne_mat <- seurat_obj@reductions$tsne@cell.embeddings
  coords_for_sling <- tsne_mat[, 1:2, drop = FALSE]
  coords_2d        <- tsne_mat[, 1:2, drop = FALSE]
  red_source       <- "tsne"
} else if ("harmony" %in% rd_names) {
  harm_mat <- seurat_obj@reductions$harmony@cell.embeddings
  coords_for_sling <- harm_mat[, 1:2, drop = FALSE]
  coords_2d        <- harm_mat[, 1:2, drop = FALSE]
  red_source       <- "harmony"
} else {
  stop("No PCA/UMAP/tSNE/Harmony reduction found in Adipo_HFD.")
}
npcs_used <- ncol(coords_for_sling)

sce <- as.SingleCellExperiment(seurat_obj)
reducedDim(sce, "DR")   <- as.matrix(coords_for_sling)
reducedDim(sce, "DR2D") <- as.matrix(coords_2d)
colData(sce)$cluster    <- cl
start_clus_arg <- if (is.na(start_cluster)) NULL else as.character(start_cluster)

sce <- slingshot(sce, clusterLabels = "cluster", reducedDim = "DR",
                 start.clus = start_clus_arg, approx_points = approx_points)

pt_mat <- slingPseudotime(sce, na = TRUE)
if (!is.null(pt_mat) && ncol(pt_mat) > 0) {
  common_cells <- rownames(pt_mat)
  for (i in seq_len(ncol(pt_mat))) {
    v <- pt_mat[common_cells, i]
    seurat_obj[[paste0("slingPT_L", i)]] <- v[colnames(seurat_obj)]
  }
  pt_mean <- apply(pt_mat, 1, function(x) mean(x, na.rm = TRUE))
  seurat_obj$slingPT_mean <- pt_mean[colnames(seurat_obj)]
  seurat_obj$slingPT_lineage <- vapply(
    1:nrow(pt_mat),
    function(r) {
      vals <- pt_mat[r, ]
      if (all(is.na(vals))) return(NA_character_)
      paste0("L", which.min(vals))
    }, character(1)
  )[colnames(seurat_obj)]
} else {
  seurat_obj$slingPT_mean    <- NA_REAL_
  seurat_obj$slingPT_lineage <- NA_character_
}

seurat_obj@misc$slingshot <- list(
  object_tag       = obj_tag,
  reduction_source = red_source,
  npcs_used        = npcs_used,
  DR               = reducedDim(sce, "DR"),
  DR2D             = reducedDim(sce, "DR2D"),
  cluster_labels   = cl,
  start_cluster    = start_clus_arg,
  pseudotime       = pt_mat,
  curve_weights    = slingCurveWeights(sce),
  lineages         = slingLineages(sce),
  mst              = slingMST(sce),
  curves           = lapply(slingCurves(sce), function(crv) {
    list(s = crv$s, ord = crv$ord, w = crv$w, lambda = crv$lambda)
  }),
  sds              = SlingshotDataSet(sce)
)

Adipo_HFD <- seurat_obj

library(Seurat)
library(tidyverse)

## ---------- ND data frames ----------
## 2D coordinates (DR2D) for ND
nd_coords <- as.data.frame(Adipo_ND@misc$slingshot$DR2D)
colnames(nd_coords)[1:2] <- c("DR1","DR2")
nd_coords$cell <- rownames(nd_coords)

## pseudotime (mean) for ND, aligned by cell
nd_pt <- Adipo_ND$slingPT_mean
nd_df <- nd_coords %>%
  mutate(pseudotime = nd_pt[cell],
         cluster = as.character(Adipo_ND$slingshot_cluster[cell]),
         condition = "ND")

## Slingshot curves for ND (take first two dims)
## curves is a list; each element has $s (curve coordinates in DR space)
nd_curves_list <- Adipo_ND@misc$slingshot$curves
nd_curves_df <- bind_rows(lapply(seq_along(nd_curves_list), function(i){
  s <- as.data.frame(nd_curves_list[[i]]$s)
  if (ncol(s) >= 2) s <- s[,1:2] else stop("ND curve has <2 dims.")
  colnames(s) <- c("DR1","DR2")
  s$lineage <- paste0("Lineage", i)
  s$condition <- "ND"
  s
}))

## ---------- HFD data frames ----------
hfd_coords <- as.data.frame(Adipo_HFD@misc$slingshot$DR2D)
colnames(hfd_coords)[1:2] <- c("DR1","DR2")
hfd_coords$cell <- rownames(hfd_coords)

hfd_pt <- Adipo_HFD$slingPT_mean
hfd_df <- hfd_coords %>%
  mutate(pseudotime = hfd_pt[cell],
         cluster = as.character(Adipo_HFD$slingshot_cluster[cell]),
         condition = "HFD")

hfd_curves_list <- Adipo_HFD@misc$slingshot$curves
hfd_curves_df <- bind_rows(lapply(seq_along(hfd_curves_list), function(i){
  s <- as.data.frame(hfd_curves_list[[i]]$s)
  if (ncol(s) >= 2) s <- s[,1:2] else stop("HFD curve has <2 dims.")
  colnames(s) <- c("DR1","DR2")
  s$lineage <- paste0("Lineage", i)
  s$condition <- "HFD"
  s
}))

both_cells <- bind_rows(nd_df, hfd_df)
both_curves <- bind_rows(nd_curves_df, hfd_curves_df)

pt_range <- range(both_cells$pseudotime, na.rm = TRUE)

x_range <- range(both_cells$DR1, na.rm = TRUE)
y_range <- range(both_cells$DR2, na.rm = TRUE)

p_compare <- ggplot(both_cells, aes(DR1, DR2, color = pseudotime)) +
  geom_point(size = 0.35, alpha = 0.9) +
  geom_path(data = both_curves,
            aes(DR1, DR2, group = interaction(condition, lineage)),
            inherit.aes = FALSE, linewidth = 0.7) +
  scale_color_viridis_c(na.value = "grey80", limits = pt_range) +
  scale_x_continuous(limits = x_range) +
  scale_y_continuous(limits = y_range) +
  coord_equal() +
  facet_wrap(~ condition, ncol = 2, scales = "fixed") +  # <- fixed
  theme_classic(base_size = 12) +
  labs(title = "Adipocyte differentiation trajectories: ND vs HFD",
       x = "DR1", y = "DR2", color = "Pseudotime")

print(p_compare)
