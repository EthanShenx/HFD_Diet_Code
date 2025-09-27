# =========== Prep ===========
setwd("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/0Subdivision")

library(Seurat)
library(dplyr)
library(plyr)
library(glmGamPoi)
library(patchwork)
library(Nebulosa)
library(SingleR)
library(celldex)
library(ggsci)
library(scico)

All <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_All_Stroma_sub.rds")

Stroma_All$combo <- paste0(Stroma_All$orig.ident, "_", Stroma_All$subcluster)

pairwise_wilcox_for_feature <- function(obj, feature, group.by = "combo") {
  dat <- FetchData(obj, vars = c(feature, group.by))
  colnames(dat) <- c("expr", "group")
  dat$group <- factor(dat$group)
  groups <- levels(dat$group)

  cmb <- t(combn(groups, 2))
  if (is.null(dim(cmb))) {

    cmb <- matrix(cmb, ncol = 2, byrow = TRUE)
  }

  rows <- vector("list", nrow(cmb))
  p_raw <- numeric(nrow(cmb))

  for (i in seq_len(nrow(cmb))) {
    g1 <- cmb[i, 1]; g2 <- cmb[i, 2]
    x  <- dat$expr[dat$group == g1]
    y  <- dat$expr[dat$group == g2]

    wt <- suppressWarnings(wilcox.test(x, y, exact = FALSE))
    p_raw[i] <- wt$p.value

    W <- unname(wt$statistic)

    cd <- NA_real_
    if (requireNamespace("effsize", quietly = TRUE)) {
      cd <- tryCatch(effsize::cliff.delta(x, y)$estimate, error = function(e) NA_real_)
    }

    rows[[i]] <- data.frame(
      feature = feature,
      group1  = g1,
      group2  = g2,
      n1      = length(x),
      n2      = length(y),
      W       = W,
      p_value = wt$p.value,
      cliffs_delta = cd,
      stringsAsFactors = FALSE
    )
  }

  res <- bind_rows(rows)
  res$p_adj_BH <- p.adjust(res$p_value, method = "BH")
  res <- res %>%
    arrange(p_adj_BH, p_value, feature, group1, group2)
  return(res)
}

pairwise_wilcox_for_features <- function(obj, features, group.by = "combo") {
  features <- unique(features)
  results <- purrr::map_dfr(features, ~pairwise_wilcox_for_feature(obj, .x, group.by = group.by))
  return(results)
}

palette_10 <- c(
  "#af9789", # muted brown
  "#e29898", # soft pink
  "#52a57f", # green
  "#7b6ca3", # purple
  "#f0a641", # orange
  "#c8735e", # terracotta
  "#af5878", # rose
  "#a0c277", # sage
  "#6671a7", # blue
  "#8abdbf"  # teal
)
grp_levels <- levels(factor(Stroma_All$combo))
n_groups   <- length(grp_levels)
cols_use <- palette_10[seq_len(min(n_groups, length(palette_10)))]

##############################
############ EGFR ###########
##############################

egfr <- VlnPlot(Stroma_All, 
        features = "Egfr", 
        pt.size = 0, 
        group.by = "combo",
        cols = cols_use) +
  geom_boxplot(
    width = 0.125,
    fill = "white",
    alpha = 1,
    outlier.size = 0.8,
    outlier.color = "black",
    outlier.shape = 19,
    linewidth = 0.3
  )

vio_layers <- which(sapply(egfr$layers, function(l) inherits(l$geom, "GeomViolin")))
for (i in vio_layers) {
  egfr$layers[[i]]$aes_params$colour <- NA
  egfr$layers[[i]]$aes_params$linewidth <- 0
  egfr$layers[[i]]$aes_params$size <- 0
}

print(egfr) # pdf: 4.77 x 2.69

res_egfr <- pairwise_wilcox_for_features(Stroma_All, 
                                         features = "Egfr", 
                                         group.by = "combo")
res_egfr
  
##############################
############ IGF1 ###########
##############################
igf1 <- VlnPlot(Stroma_All, 
        features = "Igf1", 
        pt.size = 0, 
        group.by = "combo",
        cols = cols_use) +
  geom_boxplot(
    width = 0.125,
    fill = "white",
    alpha = 1,
    outlier.size = 0.8,
    outlier.color = "black",
    outlier.shape = 19,
    linewidth = 0.3
  )

vio_layers <- which(sapply(igf1$layers, function(l) inherits(l$geom, "GeomViolin")))
for (i in vio_layers) {
  igf1$layers[[i]]$aes_params$colour <- NA
  igf1$layers[[i]]$aes_params$linewidth <- 0
  igf1$layers[[i]]$aes_params$size <- 0
}

print(igf1) # pdf: 4.77 x 2.69

res_igf1 <- pairwise_wilcox_for_features(Stroma_All, 
                                         features = "Igf1", 
                                         group.by = "combo")
res_igf1

##############################
############ TIMP3 ###########
##############################
timp3 <- VlnPlot(Stroma_All, 
        features = "Timp3", 
        pt.size = 0, 
        group.by = "combo",
        cols = cols_use) +
  geom_boxplot(
    width = 0.125,
    fill = "white",
    alpha = 1,
    outlier.size = 0.8,
    outlier.color = "black",
    outlier.shape = 19,
    linewidth = 0.3
  )

vio_layers <- which(sapply(timp3$layers, function(l) inherits(l$geom, "GeomViolin")))
for (i in vio_layers) {
  timp3$layers[[i]]$aes_params$colour <- NA
  timp3$layers[[i]]$aes_params$linewidth <- 0
  timp3$layers[[i]]$aes_params$size <- 0
}

print(timp3) # pdf: 4.77 x 2.69

res_timp3 <- pairwise_wilcox_for_features(Stroma_All, 
                                         features = "Igf1", 
                                         group.by = "combo")
res_timp3

##############################
############ FGF2 ###########
##############################
fgf2 <- VlnPlot(Stroma_All, 
        features = "Fgf2", 
        pt.size = 0, 
        group.by = "combo",
        cols = cols_use) +
  geom_boxplot(
    width = 0.125,
    fill = "white",
    alpha = 1,
    outlier.size = 0.8,
    outlier.color = "black",
    outlier.shape = 19,
    linewidth = 0.3
  )

vio_layers <- which(sapply(fgf2$layers, function(l) inherits(l$geom, "GeomViolin")))
for (i in vio_layers) {
  fgf2$layers[[i]]$aes_params$colour <- NA
  fgf2$layers[[i]]$aes_params$linewidth <- 0
  fgf2$layers[[i]]$aes_params$size <- 0
}

print(fgf2) # pdf: 4.77 x 2.69

res_fgf2 <- pairwise_wilcox_for_features(Stroma_All, 
                                         features = "Igf1", 
                                         group.by = "combo")
res_fgf2

##############################
############ FGF10 ###########
##############################
fgf10 <- VlnPlot(Stroma_All, 
        features = "Fgf10", 
        pt.size = 0, 
        group.by = "combo",
        cols = cols_use) +
  geom_boxplot(
    width = 0.125,
    fill = "white",
    alpha = 1,
    outlier.size = 0.8,
    outlier.color = "black",
    outlier.shape = 19,
    linewidth = 0.3
  )

vio_layers <- which(sapply(fgf10$layers, function(l) inherits(l$geom, "GeomViolin")))
for (i in vio_layers) {
  fgf10$layers[[i]]$aes_params$colour <- NA
  fgf10$layers[[i]]$aes_params$linewidth <- 0
  fgf10$layers[[i]]$aes_params$size <- 0
}

print(fgf10) # pdf: 4.77 x 2.69

res_fgf10 <- pairwise_wilcox_for_features(Stroma_All, 
                                         features = "Igf1", 
                                         group.by = "combo")
res_fgf10

##############################
############ FGF7 ############
##############################
fgf7 <- VlnPlot(Stroma_All, 
        features = "Fgf7", 
        pt.size = 0, 
        group.by = "combo",
        cols = cols_use) +
  geom_boxplot(
    width = 0.125,
    fill = "white",
    alpha = 1,
    outlier.size = 0.8,
    outlier.color = "black",
    outlier.shape = 19,
    linewidth = 0.3
  )

vio_layers <- which(sapply(fgf7$layers, function(l) inherits(l$geom, "GeomViolin")))
for (i in vio_layers) {
  fgf7$layers[[i]]$aes_params$colour <- NA
  fgf7$layers[[i]]$aes_params$linewidth <- 0
  fgf7$layers[[i]]$aes_params$size <- 0
}

print(fgf7) # pdf: 4.77 x 2.69

res_fgf7 <- pairwise_wilcox_for_features(Stroma_All, 
                                         features = "Igf1", 
                                         group.by = "combo")
res_fgf7

##############################
############ GH ############
##############################
gh <- VlnPlot(Stroma_All, 
        features = "Ghr", 
        pt.size = 0, 
        group.by = "combo",
        cols = cols_use) +
  geom_boxplot(
    width = 0.125,
    fill = "white",
    alpha = 1,
    outlier.size = 0.8,
    outlier.color = "black",
    outlier.shape = 19,
    linewidth = 0.3
  )

vio_layers <- which(sapply(gh$layers, function(l) inherits(l$geom, "GeomViolin")))
for (i in vio_layers) {
  gh$layers[[i]]$aes_params$colour <- NA
  gh$layers[[i]]$aes_params$linewidth <- 0
  gh$layers[[i]]$aes_params$size <- 0
}

print(gh) # pdf: 4.77 x 2.69

res_gh <- pairwise_wilcox_for_features(Stroma_All, 
                                         features = "Igf1", 
                                         group.by = "combo")
res_gh

##########################################################################################
########################################### SI ###########################################
##########################################################################################

##############################
############ FGF1 ###########
##############################
fgf1 <- VlnPlot(Stroma_All, 
        features = "Fgf1", 
        pt.size = 0, 
        group.by = "combo",
        cols = cols_use) +
  geom_boxplot(
    width = 0.125,
    fill = "white",
    alpha = 1,
    outlier.size = 0.8,
    outlier.color = "black",
    outlier.shape = 19,
    linewidth = 0.3
  )

vio_layers <- which(sapply(fgf1$layers, function(l) inherits(l$geom, "GeomViolin")))
for (i in vio_layers) {
  fgf1$layers[[i]]$aes_params$colour <- NA
  fgf1$layers[[i]]$aes_params$linewidth <- 0
  fgf1$layers[[i]]$aes_params$size <- 0
}

print(fgf1) # pdf: 4.77 x 2.69

res_fgf1 <- pairwise_wilcox_for_features(Stroma_All, 
                                         features = "Igf1", 
                                         group.by = "combo")
res_fgf1

##############################
############ FGF9 ###########
##############################
fgf9 <- VlnPlot(Stroma_All, 
        features = "Fgf9", 
        pt.size = 0, 
        group.by = "combo",
        cols = cols_use) +
  geom_boxplot(
    width = 0.125,
    fill = "white",
    alpha = 1,
    outlier.size = 0.8,
    outlier.color = "black",
    outlier.shape = 19,
    linewidth = 0.3
  )

vio_layers <- which(sapply(fgf9$layers, function(l) inherits(l$geom, "GeomViolin")))
for (i in vio_layers) {
  fgf9$layers[[i]]$aes_params$colour <- NA
  fgf9$layers[[i]]$aes_params$linewidth <- 0
  fgf9$layers[[i]]$aes_params$size <- 0
}

print(fgf9) # pdf: 4.77 x 2.69

res_fgf9 <- pairwise_wilcox_for_features(Stroma_All, 
                                         features = "Igf1", 
                                         group.by = "combo")
res_fgf9

##############################
############ FGF11 ###########
##############################
fgf11 <- VlnPlot(Stroma_All, 
        features = "Fgf11", 
        pt.size = 0, 
        group.by = "combo",
        cols = cols_use) +
  geom_boxplot(
    width = 0.125,
    fill = "white",
    alpha = 1,
    outlier.size = 0.8,
    outlier.color = "black",
    outlier.shape = 19,
    linewidth = 0.3
  )

vio_layers <- which(sapply(fgf11$layers, function(l) inherits(l$geom, "GeomViolin")))
for (i in vio_layers) {
  fgf11$layers[[i]]$aes_params$colour <- NA
  fgf11$layers[[i]]$aes_params$linewidth <- 0
  fgf11$layers[[i]]$aes_params$size <- 0
}

print(fgf11) # pdf: 4.77 x 2.69

res_fgf11 <- pairwise_wilcox_for_features(Stroma_All, 
                                         features = "Igf1", 
                                         group.by = "combo")
res_fgf11

##############################
############ FGF13 ###########
##############################
fgf13 <- VlnPlot(Stroma_All, 
        features = "Fgf13", 
        pt.size = 0, 
        group.by = "combo",
        cols = cols_use) +
  geom_boxplot(
    width = 0.125,
    fill = "white",
    alpha = 1,
    outlier.size = 0.8,
    outlier.color = "black",
    outlier.shape = 19,
    linewidth = 0.3
  )

vio_layers <- which(sapply(fgf13$layers, function(l) inherits(l$geom, "GeomViolin")))
for (i in vio_layers) {
  fgf13$layers[[i]]$aes_params$colour <- NA
  fgf13$layers[[i]]$aes_params$linewidth <- 0
  fgf13$layers[[i]]$aes_params$size <- 0
}

print(fgf13) # pdf: 4.77 x 2.69

res_fgf13 <- pairwise_wilcox_for_features(Stroma_All, 
                                         features = "Igf1", 
                                         group.by = "combo")
res_fgf13

##############################
############ FGF18 ###########
##############################
fgf18 <- VlnPlot(Stroma_All, 
        features = "Fgf18", 
        pt.size = 0, 
        group.by = "combo",
        cols = cols_use) +
  geom_boxplot(
    width = 0.125,
    fill = "white",
    alpha = 1,
    outlier.size = 0.8,
    outlier.color = "black",
    outlier.shape = 19,
    linewidth = 0.3
  )

vio_layers <- which(sapply(fgf18$layers, function(l) inherits(l$geom, "GeomViolin")))
for (i in vio_layers) {
  fgf18$layers[[i]]$aes_params$colour <- NA
  fgf18$layers[[i]]$aes_params$linewidth <- 0
  fgf18$layers[[i]]$aes_params$size <- 0
}

print(fgf18) # pdf: 4.77 x 2.69

res_fgf18 <- pairwise_wilcox_for_features(Stroma_All, 
                                         features = "Igf1", 
                                         group.by = "combo")
res_fgf18
