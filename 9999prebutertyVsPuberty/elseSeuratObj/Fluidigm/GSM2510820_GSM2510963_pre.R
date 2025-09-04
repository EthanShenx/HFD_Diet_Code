library(GEOquery)
library(Seurat)
library(Matrix)
library(tidyverse)
library(org.Mm.eg.db)
library(rlang)
library(DoubletFinder)

setwd("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/9999prebutertyVsPuberty/elseSeuratObj/Fluidigm")
data_dir <- "/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/9999prebutertyVsPuberty/elseSeuratObj/Fluidigm/data"
dir.create(data_dir, showWarnings = FALSE)

gsm_ids <- paste0("GSM2510", 820:963)

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
                                 project = "GSM2510820-GSM2510963",
                                 min.cells = 3,
                                 min.features = 200)

seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")

GSM2510820_GSM2510963_pre <- seurat_obj

GSM2510820_GSM2510963_pre@meta.data$stage <- "pre-puberty"
GSM2510820_GSM2510963_pre@meta.data$platform <- "fluidigm"

GSM2510820_GSM2510963_pre@meta.data$stage <- as.factor(GSM2510820_GSM2510963_pre@meta.data$stage)
GSM2510820_GSM2510963_pre@meta.data$platform <- as.factor(GSM2510820_GSM2510963_pre@meta.data$platform)

saveRDS(GSM2510820_GSM2510963_pre, file = "GSM2510820_GSM2510963_pre.rds")
