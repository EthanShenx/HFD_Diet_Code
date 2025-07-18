library(Seurat)
library(harmony)
library(patchwork)
library(dplyr)
library(ggplot2)
library(sctransform)
library(stringr)
library(ggrepel)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

ND <- readRDS('D:/data/23BMI/ND_HFD_MG_snRNAseq/7.1data/ND_corrected.rds')
HFD <- readRDS('D:/data/23BMI/ND_HFD_MG_snRNAseq/7.1data/HFD_corrected.rds')
HFD@meta.data$cell_type <- Idents(HFD)
ND@meta.data$cell_type <- Idents(ND)

cell_types_to_keep <- c("Adipo", "Stroma","Immune", "LumProg", "Endo", "Basal", "HormSens")
ND <- subset(ND, idents = cell_types_to_keep)

ND <- RenameCells(ND, add.cell.id = "ND")
HFD <- RenameCells(HFD, add.cell.id = "HFD")
data <- list(ND, HFD)

scobj <- merge(x=data[[1]], y=data[-1], project="diet")
scobj <- NormalizeData(scobj)
scobj <- FindVariableFeatures(scobj, selection.method = "vst", nfeatures = 2000)
scobj <- ScaleData(scobj, features = VariableFeatures(scobj))
scobj <- RunPCA(scobj, features = VariableFeatures(object = scobj))
ElbowPlot(scobj)
scobj <- RunHarmony(scobj, reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
scobj <- RunUMAP(scobj, reduction = "harmony", dims = 1:15,reduction.name = "umap")

scobj$named_cluster <- paste0(scobj$orig.ident, "_", Idents(scobj))
Idents(scobj) <- "named_cluster"

DimPlot(scobj, reduction = "umap",group.by = "orig.ident")
DimPlot(scobj, reduction = "umap",group.by = "named_cluster")

saveRDS(scobj, file = "D:/data/23BMI/ND_HFD_MG_snRNAseq/harmony_all.rds")

devtools::install_github('immunogenomics/presto')
Adipo.DEGs <- FindMarkers(scobj, 
                          ident.1 = "HFD_Adipo",
                          ident.2 = "ND_Adipo") 

Stroma.DEGs <- FindMarkers(scobj, 
                           ident.1 = "HFD_Stroma",
                           ident.2 = "ND_Stroma") 
Immune.DEGs <- FindMarkers(scobj, 
                           ident.1 = "HFD_Immune",
                           ident.2 = "ND_Immune") 
Endo.DEGs <- FindMarkers(scobj, 
                         ident.1 = "HFD_Endo",
                         ident.2 = "ND_Endo") 
LumProg.DEGs <- FindMarkers(scobj, 
                            ident.1 = "HFD_LumProg",
                            ident.2 = "ND_LumProg") 
HormSens.DEGs <- FindMarkers(scobj, 
                             ident.1 = "HFD_HormSens",
                             ident.2 = "ND_HormSens") 
Basal.DEGs <- FindMarkers(scobj, 
                          ident.1 = "HFD_Basal",
                          ident.2 = "ND_Basal") 

Adipo.DEGs$cell_type <- "Adipo"
Stroma.DEGs$cell_type <- "Stroma"
Basal.DEGs$cell_type <- "Basal"
HormSens.DEGs$cell_type <- "HormSens"
LumProg.DEGs$cell_type <- "LumProg"
Endo.DEGs$cell_type <- "Endo"
Immune.DEGs$cell_type <- "Immune"

write.csv(Immune.DEGs, row.names = T, file = "Immune_DEGs.csv")
write.csv(Endo.DEGs, row.names = T, file = "Endo_DEGs.csv")
write.csv(Adipo.DEGs, row.names = T, file = "Adipo_DEGs.csv")
write.csv(Basal.DEGs, row.names = T, file = "Basal_DEGs.csv")
write.csv(HormSens.DEGs, row.names = T, file = "HormSens_DEGs.csv")
write.csv(LumProg.DEGs, row.names = T, file = "LumProg_DEGs.csv")
write.csv(Stroma.DEGs, row.names = T, file = "Stroma_DEGs.csv")

