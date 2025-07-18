setwd("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/6downstreamEffects")
library(Seurat)
library(biomaRt)
library(msigdbr)
library(AUCell)
library(dplyr)
library(patchwork)
library(tibble)
library(pheatmap)

ND <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_ND.rds")
HFD <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_HFD.rds")
Idents(ND) <- "cell_type"
Idents(HFD) <- "cell_type"

# 获取 Reactome 人类基因集
msig_mouse <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:REACTOME")
hippo_mouse <- msig_mouse %>%
  filter(gs_name == "REACTOME_SIGNALING_BY_HIPPO") %>%
  pull(gene_symbol) %>%
  unique()
fgfr_mouse <- msig_mouse %>%
  filter(gs_name == "REACTOME_SIGNALING_BY_FGFR") %>%
  pull(gene_symbol) %>%
  unique()

hallmarks <- c(
  PI3K_AKT  = "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
  MAPK_ERK  = "HALLMARK_KRAS_SIGNALING_UP",
  JAK_STAT  = "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  mTORC1    = "HALLMARK_MTORC1_SIGNALING",
  MYC       = "HALLMARK_MYC_TARGETS_V1",
  G2M_E2F   = "HALLMARK_G2M_CHECKPOINT"
)

msig_db <- msigdbr(species = "Mus musculus", collection = "H")
gene_sets <- lapply(hallmarks, \(gs)
msig_db |>
  filter(gs_name == gs) |>
  pull(gene_symbol) |>
  unique())
gene_sets <- c(
  gene_sets,
  list(HIPPO = hippo_mouse, FGFR2 = fgfr_mouse)
)

ND <- AddModuleScore(
  object   = ND,
  features = gene_sets,
  name     = names(gene_sets)
)

exprMat <- GetAssayData(ND, layer = "data") # v5 用 layer 参数
rankings <- AUCell_buildRankings(exprMat, plotStats = FALSE)
auc <- AUCell_calcAUC(gene_sets, rankings)
ND[[]] <- cbind(ND[[]], t(getAUC(auc))) # 写回 meta.data

##############################################################################
## 3. 可视化：UMAP & Violin
##############################################################################
FeaturePlot(
  ND,
  features = c(
    "PI3K_AKT1", "MAPK_ERK2", "JAK_STAT3",
    "mTORC14", "MYC5", "G2M_E2F6", "HIPPO", "FGFR2"
  ),
  order = TRUE, cols = c("grey90", "red")
)

VlnPlot(
  ND,
  features = c(
    "PI3K_AKT1", "MAPK_ERK2", "JAK_STAT3",
    "mTORC14", "MYC5", "G2M_E2F6", "HIPPO", "FGFR2"
  ),
  group.by = "cell_type", pt.size = 0
) + ggplot2::theme(legend.position = "none")


### For HFD ###

HFD <- AddModuleScore(
  object   = HFD,
  features = gene_sets,
  name     = names(gene_sets)
)

exprMat <- GetAssayData(HFD, layer = "data") # v5 用 layer 参数
rankings <- AUCell_buildRankings(exprMat, plotStats = FALSE)
auc <- AUCell_calcAUC(gene_sets, rankings)
HFD[[]] <- cbind(HFD[[]], t(getAUC(auc))) # 写回 meta.data

score_cols <- c(
  "PI3K_AKT1", "MAPK_ERK2", "JAK_STAT3",
  "mTORC14", "MYC5", "G2M_E2F6", "HIPPO", "FGFR2"
)

###### Comparison
library(dplyr)
library(tibble)
library(pheatmap)

## 获取 ND 均值
nd_avg <- ND@meta.data %>%
  group_by(cell_type) %>%
  summarise_at(vars(all_of(score_cols)), ~ mean(.x, na.rm = TRUE)) %>%
  column_to_rownames("cell_type") %>%
  as.matrix()

## 获取 HFD 均值
hfd_avg <- HFD@meta.data %>%
  group_by(cell_type) %>%
  summarise_at(vars(all_of(score_cols)), ~ mean(.x, na.rm = TRUE)) %>%
  column_to_rownames("cell_type") %>%
  as.matrix()

## 两组有无完全一致的 cell_type？
intersect_types <- intersect(rownames(nd_avg), rownames(hfd_avg))
nd_avg <- nd_avg[intersect_types, , drop = FALSE]
hfd_avg <- hfd_avg[intersect_types, , drop = FALSE]

## 计算 HFD - ND 差值
delta_mat <- hfd_avg - nd_avg

## 画热图（红高、蓝低表示 HFD 高于 ND）
pheatmap::pheatmap(
  delta_mat,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  scale = "row", # 不建议再行归一化，否则看不清差异
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  border_color = NA,
  main = "HFD - ND Pathway Scores by Cell Type"
)
