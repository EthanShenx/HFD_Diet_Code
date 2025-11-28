library(Seurat)
library(CellChat)
library(dplyr)
library(future)
library(biomaRt)

cellchat <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/10SpatialTranscriptomics/CellChat_object/ST-patient-79_cellchat.rds")

old.cluster.name <- levels(cellchat@idents)
new.cluster.name <- c(
  "Fibroblast",
  "B",
  "Pericyte",
  "Lymph Endo",
  "Vas Endo",
  "Lum Epi",
  "Basal Epi",
  "Adipo",
  "unknown"
)
cellchat <- updateClusterLabels(
  object            = cellchat,
  old.cluster.name  = old.cluster.name,     # 也可以省略，用默认 levels(object@idents)
  new.cluster.name  = new.cluster.name,
  new.cluster.metaname = "celltype_new"    # 新的meta列名，可自定义
)

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = rowSums(cellchat@net$count), weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = rowSums(cellchat@net$weight), weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

netVisual_heatmap(cellchat, measure = "count", color.heatmap = "Blues")

pathways.show <- c("PDGF")

par(mfrow=c(1,1), xpd = TRUE)

netVisual_aggregate(cellchat, 
                    signaling = pathways.show, 
                    layout = "circle")

netVisual_aggregate(cellchat, 
                    signaling = pathways.show, 
                    layout = "spatial", 
                    edge.width.max = 2,
                    vertex.size.max = 0.8, 
                    alpha.image = 0.2, 
                    vertex.label.cex = 3.5,
                    remove.isolate = T,
                    top = 0.7)

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

par(mfrow=c(1,1))
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

levels(cellchat@idents)
table(cellchat@idents)

