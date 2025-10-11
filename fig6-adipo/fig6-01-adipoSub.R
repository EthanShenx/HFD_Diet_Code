harmony_all <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_All_sub_sub.rds")

library(Seurat)
library(dplyr)
library(plyr)
library(glmGamPoi)
library(patchwork)
library(Nebulosa)
library(SingleR)
library(celldex)

Adipo_All <- subset(harmony_all, subset = subcluster == "Adipo")

Idents(Adipo_All) <- "orig.ident"
Adipo_All_ND <- subset(Adipo_All, idents = "ND")
Adipo_All_HFD <- subset(Adipo_All, idents = "HFD")
saveRDS(Adipo_All_ND, 
        "/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_ND_Adipo_sub.rds")
saveRDS(Adipo_All_HFD,
        "/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_HFD_Adipo_sub.rds")
saveRDS(Adipo_All,
        "/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_All_Adipo_sub.rds")
