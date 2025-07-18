setwd("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/4CellCellCommunication/03CellCall")

###### library ######
library(dplyr)
library(magrittr)
library(writexl)
library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(NMF)
library(presto)
library(ggalluvial)
library(reticulate)
library(gridExtra)
library(uwot)
library(cellcall)

###### ND ######

###### object convertion ######
ND_seurat_object <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_ND.rds")

ND_cellcall <- CreateObject_fromSeurat(
  Seurat.object = ND_seurat_object,
  slot = "counts",
  cell_type = "cell_type",
  data_source = "UMI",
  scale.factor = 10^6,
  Org = "Mus musculus"
)

###### cell-cell communication score ######
ND_cellcall <- TransCommuProfile(
  object = ND_cellcall,
  pValueCor = 0.05,
  CorValue = 0.1,
  topTargetCor = 1,
  p.adjust = 0.05,
  use.type = "median",
  probs = 0.9,
  method = "weighted",
  IS_core = TRUE,
  Org = "Mus musculus"
)

n <- ND_cellcall@data$expr_l_r_log2_scale

###### cell-cell hyper-pathway activity ######
pathway.hyper.list <- lapply(colnames(n), function(i) {
  print(i)
  tmp <- getHyperPathway(data = n, object = ND_cellcall, cella_cellb = i, Org = "Mus musculus")
  return(tmp)
})
###### hyper-pathway visualization ######
myPub.df <- getForBubble(pathway.hyper.list, cella_cellb = colnames(n))
myPub.df <- subset(myPub.df, grepl(".+-LumProg$", cc))
plot1 <- plotBubble(myPub.df)
plot1

###### heatmap: L-R pair visualization between cell types ######
viewPheatmap(
  object = ND_cellcall, slot = "expr_l_r_log2_scale", show_rownames = T,
  show_colnames = T, treeheight_row = 0, treeheight_col = 10,
  cluster_rows = T, cluster_cols = F, fontsize = 12, angle_col = "45",
  main = "score"
)
###### TF enrichment plot ######
lumprog.tf <- names(ND_cellcall@data$gsea.list$LumProg@geneSets)
lumprog.tf
getGSEAplot(
  gsea.list = ND_cellcall@data$gsea.list,
  geneSetID = c("E2f3", "Elk1", "Irf1", "Irf4", "Nr1h3", "Ppargc1a", "Rxra", "Sp1", "Srebf1", "Stat1", "Stat5b"),
  myCelltype = "Adipo",
  fc.list = ND_cellcall@data$fc.list,
  selectedGeneID = ND_cellcall@data$gsea.list$Adipo@geneSets$Rb1[1:10],
  mycol = NULL
)

###### TGs activatibility to TFs ######
egmt <- ND_cellcall@data$gsea.list$LumProg

egmt.df <- data.frame(egmt)
head(egmt.df[, 1:6])
flag.index <- which(egmt.df$p.adjust < 0.05)

ridgeplot.DIY(
  x = egmt, fill = "p.adjust", showCategory = flag.index, core_enrichment = T,
  orderBy = "NES", decreasing = FALSE
)

###### HFD ######

###### object convertion ######
HFD_seurat_object <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_HFD.rds")

HFD_cellcall <- CreateObject_fromSeurat(
  Seurat.object = HFD_seurat_object,
  slot = "counts",
  cell_type = "cell_type",
  data_source = "UMI",
  scale.factor = 10^6,
  Org = "Mus musculus"
)

###### cell-cell communication score ######
HFD_cellcall <- TransCommuProfile(
  object = HFD_cellcall,
  pValueCor = 0.05,
  CorValue = 0.1,
  topTargetCor = 1,
  p.adjust = 0.05,
  use.type = "median",
  probs = 0.9,
  method = "weighted",
  IS_core = TRUE,
  Org = "Mus musculus"
)

n <- HFD_cellcall@data$expr_l_r_log2_scale

###### cell-cell hyper-pathway activity ######
pathway.hyper.list <- lapply(colnames(n), function(i) {
  print(i)
  tmp <- getHyperPathway(data = n, object = HFD_cellcall, cella_cellb = i, Org = "Mus musculus")
  return(tmp)
})

###### hyper-pathway visualization ######
myPub.df <- getForBubble(pathway.hyper.list, cella_cellb = colnames(n))
myPub.df <- subset(myPub.df, grepl(".+-LumProg$", cc))
plot1 <- plotBubble(myPub.df)
plot1

###### heatmap: L-R pair visualization between cell types ######
viewPheatmap(
  object = HFD_cellcall, slot = "expr_l_r_log2_scale", show_rownames = T,
  show_colnames = T, treeheight_row = 0, treeheight_col = 10,
  cluster_rows = T, cluster_cols = F, fontsize = 12, angle_col = "45",
  main = "score"
)
###### TF enrichment plot ######
lumprog.tf <- names(HFD_cellcall@data$gsea.list$LumProg@geneSets)
lumprog.tf
getGSEAplot(
  gsea.list = HFD_cellcall@data$gsea.list,
  geneSetID = c("E2f3", "Elk1", "Irf1", "Irf4", "Nr1h3", "Ppargc1a", "Rxra", "Sp1", "Srebf1", "Stat1", "Stat5b"),
  myCelltype = "Adipo",
  fc.list = HFD_cellcall@data$fc.list,
  selectedGeneID = HFD_cellcall@data$gsea.list$Adipo@geneSets$Rb1[1:10],
  mycol = NULL
)

###### TGs activatibility to TFs ######
egmt <- HFD_cellcall@data$gsea.list$LumProg

egmt.df <- data.frame(egmt)
head(egmt.df[, 1:6])
flag.index <- which(egmt.df$p.adjust < 0.05)

ridgeplot.DIY(
  x = egmt, fill = "p.adjust", showCategory = flag.index, core_enrichment = T,
  orderBy = "NES", decreasing = FALSE
)
