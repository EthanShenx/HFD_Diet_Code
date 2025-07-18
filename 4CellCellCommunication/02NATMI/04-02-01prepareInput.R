setwd("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/4CellCellCommunication/02NATMI")
ND <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_ND.rds")
HFD <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_HFD.rds")
library(Seurat)
library(dplyr)
library(tidyverse)

write.csv(100 * (exp(as.matrix(GetAssayData(object = ND, assay = "RNA", slot = "data"))) - 1), "ND_expressionMatrix.csv", row.names = T) 

write.csv(100 * (exp(as.matrix(GetAssayData(object = HFD, assay = "RNA", slot = "data"))) - 1), "HFD_expressionMatrix.csv", row.names = T) 

Idents(ND) <- "cell_type"
Idents(HFD) <- "cell_type"

write.csv(Idents(object = ND),"ND_metadata.csv", row.names = T)

write.csv(Idents(object = HFD),"HFD_metadata.csv", row.names = T)
