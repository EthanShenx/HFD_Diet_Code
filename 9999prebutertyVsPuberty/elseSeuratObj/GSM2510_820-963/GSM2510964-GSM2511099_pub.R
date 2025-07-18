library(GEOquery)
library(Seurat)
library(Matrix)
library(tidyverse)
library(org.Mm.eg.db)
library(rlang)
library(DoubletFinder)

setwd("/Users/coellearth/Desktop/HFD_Paper/elseSeuratObj/GSM2510_820-963")
data_dir <- "/Users/coellearth/Desktop/HFD_Paper/elseSeuratObj/GSM2510_820-963"
dir.create(data_dir, showWarnings = FALSE)

gsm_ids <- c(
  paste0("GSM2510", 964:999),
  paste0("GSM2511", sprintf("%03d", 0:99))
)

options(timeout = 300)

for (gsm in gsm_ids) {
  message("Downloading: ", gsm)
  
  gsm_data <- tryCatch(getGEO(gsm, GSEMatrix = FALSE), error = function(e) NULL)
  if (is.null(gsm_data)) next
  
  supp <- gsm_data@header$supplementary_file
  if (length(supp) == 0 || !grepl("\\.txt\\.gz$", supp)) next
  
  supp <- sub("^ftp://", "https://", supp)  
  
  dest_dir <- file.path(data_dir, gsm)
  dir.create(dest_dir, showWarnings = FALSE)
  dest_file <- file.path(dest_dir, basename(supp))
  
  if (!file.exists(dest_file)) {
    download.file(supp, destfile = dest_file, mode = "wb", method = "auto")
  }
}

expr_list <- list()

gsm_dirs <- list.dirs(data_dir, recursive = FALSE)

for (dir in gsm_dirs) {
  txt_file <- list.files(dir, pattern = "\\.txt\\.gz$", full.names = TRUE)
  if (length(txt_file) == 0) next
  
  df <- read.table(gzfile(txt_file[1]), header = TRUE, sep = "\t")
  gene_ids <- as.character(df[, 1])
  counts <- df[, 3]
  sample_name <- basename(dir)
  
  df_expr <- data.frame(gene_id = gene_ids,
                        counts,
                        stringsAsFactors = FALSE)
  colnames(df_expr)[2] <- sample_name  
  expr_list[[sample_name]] <- df_expr
}

expr_mat <- purrr::reduce(expr_list, 
                          full_join, 
                          by = "gene_id")

rownames(expr_mat) <- expr_mat$gene_id

expr_mat <- expr_mat[, -1]

expr_mat[is.na(expr_mat)] <- 0

expr_mat <- as.matrix(expr_mat)

gene_symbols <- mapIds(org.Mm.eg.db,
                       keys = rownames(expr_mat),
                       column = "SYMBOL",
                       keytype = "ENTREZID",
                       multiVals = "first")

expr_mat <- expr_mat[!is.na(gene_symbols), ]
rownames(expr_mat) <- gene_symbols[rownames(expr_mat)]

seurat_obj <- CreateSeuratObject(counts = expr_mat,
                                 project = "GSM2510964-GSM2511099",
                                 min.cells = 3,
                                 min.features = 200)

seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")

seurat_obj <- subset(seurat_obj, 
                     subset = nFeature_RNA > 200 & 
                       nFeature_RNA < 6000 & 
                       percent.mt < 6)
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, 
                                   selection.method = "vst", 
                                   nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, 
                     features = VariableFeatures(object = seurat_obj))

seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.4)

GSM2510964_GSM2511099_pub <- seurat_obj

GSM2510964_GSM2511099_pub@meta.data$stage <- "puberty"

GSM2510964_GSM2511099_pub@meta.data$stage <- as.factor(GSM2510964_GSM2511099_pub@meta.data$stage)

saveRDS(GSM2510964_GSM2511099_pub, file = "GSM2510964-GSM2511099_pub.rds")

message("✅ Seurat object successfully created and saved as 'Fluidigm_Seurat.rds'")
