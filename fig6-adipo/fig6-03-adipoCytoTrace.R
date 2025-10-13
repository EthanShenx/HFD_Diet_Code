setwd("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/fig6-adipo")

library(Seurat)
library(scater)
library(RColorBrewer)
library(CytoTRACE2) 
library(dplyr)
library(plyr)
library(glmGamPoi)
library(patchwork)
library(Nebulosa)
library(SingleR)
library(celldex)

harmony_all <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_All_sub.rds")

All <- subset(harmony_all, subset = subcluster == "Adipo" | subcluster == "Stroma")

Idents(All) <- "orig.ident"
ND <- subset(All, idents = "ND")
HFD <- subset(All, idents = "HFD")

ND_expression_data <- Seurat::GetAssayData(ND, 
                                   slot = "counts")
HFD_expression_data <- Seurat::GetAssayData(HFD, 
                                   slot = "counts")
All_expression_data <- Seurat::GetAssayData(All, 
                                   slot = "counts")

####################################################
####################################################
####################################################

# running CytoTRACE 2 main function - cytotrace2 - with default parameters
ND_cytotrace2_result <- cytotrace2(ND_expression_data)

# extract annotation data
ND_annotation <- ND$subcluster

# generate prediction and phenotype association plots with plotData function
ND_plots <- plotData(cytotrace2_result = ND_cytotrace2_result, 
                     is_seurat = T,
                  expression_data = ND_expression_data
                  )

ND_plots$CytoTRACE2_Potency_UMAP
ND_plots$CytoTRACE2_Relative_UMAP

ND_expression_data <- Seurat::GetAssayData(Adipo_ND, 
                                   slot = "counts")

####################################################
####################################################
####################################################

# running CytoTRACE 2 main function - cytotrace2 - with default parameters
HFD_cytotrace2_result <- cytotrace2(HFD_expression_data)

# extract annotation data
HFD_annotation <- HFD$subcluster

# generate prediction and phenotype association plots with plotData function
HFD_plots <- plotData(cytotrace2_result = HFD_cytotrace2_result,
                  expression_data = HFD_expression_data
                  )

HFD_plots$CytoTRACE2_Potency_UMAP
HFD_plots$CytoTRACE2_Relative_UMAP

####################################################
####################################################
####################################################

# running CytoTRACE 2 main function - cytotrace2 - with default parameters
HFD_cytotrace2_result <- cytotrace2(HFD_expression_data)

# extract annotation data
HFD_annotation <- HFD$subcluster

# generate prediction and phenotype association plots with plotData function
HFD_plots <- plotData(cytotrace2_result = HFD_cytotrace2_result, 
                     is_seurat = F,
                  expression_data = HFD_expression_data
                  )

HFD_plots$CytoTRACE2_Potency_UMAP
HFD_plots$CytoTRACE2_Relative_UMAP

ND_plots$CytoTRACE2_Potency_UMAP / HFD_plots$CytoTRACE2_Potency_UMAP
