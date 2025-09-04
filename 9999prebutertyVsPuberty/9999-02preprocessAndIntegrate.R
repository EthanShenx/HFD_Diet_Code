setwd("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/9999prebutertyVsPuberty")

library(Seurat)
library(harmony)
library(dplyr)
library(patchwork) 
library(ggplot2)
library(harmony)
library(glue)
library(future)
library(Matrix)
set.seed(1234)

rds_files <- list.files(
  path = "./elseSeuratObj",
  pattern = "\\.rds$",
  recursive = TRUE,
  full.names = TRUE
)

for (f in rds_files) {
  objname <- tools::file_path_sans_ext(basename(f))
  assign(objname, readRDS(f))
}

objs <- list(
  GSM2510820_2w_C1_prepub = GSM2510820_GSM2510963_pre,
  GSM2510964_4w9_C1_pub = GSM2510964_GSM2511099_pub,
  GSM2512829_4w9_C1_pub  = GSM2512829_GSM2513049_pub,
  GSM4994963_2w5_10x_prepub = GSM4994963_pre,
  GSM4994964_5w_TEB_10x_pub = GSM4994964_pub,
  GSM4994965_5w_Duct_10x_pub = GSM4994965_pub,
  GSM2759554_5w_10x_pub = GSM2759554_pub,
  GSM2759555_5w_10x_pub = GSM2759555_pub
)

objs <- lapply(objs, function(o){
  DefaultAssay(o) <- "RNA"
  o[["percent.mt"]] <- PercentageFeatureSet(o, pattern="^mt-|^MT-")
  o
})

meta_from_name <- function(nm){
  parts <- strsplit(nm, "_")[[1]]
  list(
    sample_id = nm,
    week = parts[2],
    platform = ifelse(grepl("C1", nm), "C1", "10x"),
    stage = ifelse(grepl("2w|2w5|Pre", nm, ignore.case=TRUE), "Prepuberty","Puberty")
  )
}

objs <- lapply(names(objs), function(nm){
  o <- objs[[nm]]
  m <- meta_from_name(nm)
  o$sample_id <- m$sample_id
  o$platform  <- m$platform
  o$stage     <- m$stage
  o$batch     <- paste(m$platform, m$stage, sep="|")
  o
})

names(objs) <- sapply(objs, \(o) unique(o$sample_id))

objs <- lapply(objs, function(o){
  upper <- ifelse(unique(o$platform)[1] == "C1", 25000, 8000)
  subset(o, subset = nFeature_RNA >= 200 & nFeature_RNA <= upper & percent.mt < 20)
})


merged <- Reduce(function(x,y) merge(x, y), objs)
DefaultAssay(merged) <- "RNA"
merged <- NormalizeData(merged)
merged <- FindVariableFeatures(merged, nfeatures = 3000)
merged <- ScaleData(merged, vars.to.regress = c("nCount_RNA", "percent.mt"))
merged <- RunPCA(merged, npcs = 60)
# ElbowPlot(merged, ndims = 60)

merged <- RunHarmony(merged, group.by.vars = "sample_id", dims = 1:60)
merged <- RunUMAP(merged, reduction = "harmony", dims = 1:50)
merged <- FindNeighbors(merged, reduction = "harmony", dims = 1:50)
merged <- FindClusters(merged, dims = 1:30, resolution = 0.001)

DimPlot(merged, group.by = "sample_id", repel = TRUE, alpha = 0.5)
DimPlot(merged, group.by = "platform", repel = TRUE, alpha = 0.5)
DimPlot(merged, group.by = "stage")
