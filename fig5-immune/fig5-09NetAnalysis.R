# =========================
# CellChat centrality: role-specific heatmaps for ND / HFD
# =========================
library(CellChat)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)   # 用于更好看的热图；若没装可换 pheatmap
library(circlize)

setwd("/Users/renyixiang/Desktop/data/23BMI/ND_HFD_MG_snRNAseq/Figure7_immune/CellChat")

# ---- 读取对象 ----
HFD <- readRDS("cellchat_HFD_mac_state.rds")
ND  <- readRDS("cellchat_ND_mac_state.rds")

# ---- 你要看的 pathway（示例，自己替换）----
# 例：pathways.show <- c("NRG", "EGF", "TNF")
pathways.show <- c("ADGRL","Cholesterol","APP","DHEA","FN1","NRG","KIT","FGF")  # 如果设为 NULL，就展示所有已计算的通路

# ---- 角色设定 ----
senders_keep   <- c("MacS01", "MacS03")
receivers_keep <- c("Basal", "LumProg", "HormSens")

# 根据对象的分群名自动推断 "其他细胞"（作为 mediator / influencer）
others_of <- function(cellchat){
  all_groups <- levels(cellchat@idents$group)
  setdiff(all_groups, union(senders_keep, receivers_keep))
}

# ---- 计算中心性（如果已经计算过也可重复运行，CellChat 会覆盖存储）----
compute_cent <- function(cc){
  cc <- netAnalysis_computeCentrality(cc, slot.name = "netP")
  return(cc)
}
HFD <- compute_cent(HFD)
ND  <- compute_cent(ND)

netAnalysis_signalingRole_network(ND, signaling = pathways.show)
netAnalysis_signalingRole_network(HFD, signaling = pathways.show)

save.image("NetAnalysis.RData")
