library(Seurat)
library(CellChat)
library(dplyr)
library(future)
library(biomaRt)

plan(multicore, workers = 8)

cellchat_sc <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/10SpatialTranscriptomics/CellChat_object/cellchat_human_atlas.rds")

groupSize <- as.numeric(table(cellchat_sc@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat_sc@net$count, 
                 vertex.weight = groupSize, 
                 weight.scale = T, 
                 label.edge= F, 
                 title.name = "Number of interactions")
netVisual_circle(cellchat_sc@net$weight, 
                 vertex.weight = groupSize, 
                 weight.scale = T, 
                 label.edge= F, 
                 title.name = "Interaction weights/strength")

par(mfrow=c(1,1))
netVisual_heatmap(cellchat_sc, signaling = pathways.show, color.heatmap = "Reds")
#> Do heatmap based on a single object
#> 

netVisual_bubble(cellchat_sc, sources.use = "fibro-matrix", signaling = c("COLLAGEN"), remove.isolate = FALSE)

netVisual_aggregate(cellchat_sc, 
                    signaling = "COLLAGEN", 
                    layout = "circle")

netVisual_heatmap(cellchat_sc, signaling = "COLLAGEN", color.heatmap = "Reds")

mat <- cellchat_sc@net$weight
par(mfrow = c(4,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
