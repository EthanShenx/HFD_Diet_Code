library(Seurat)
library(Matrix)
library(dplyr)

ND <- readRDS("D:/data/23BMI/ND_HFD_MG_snRNAseq/Figure7_immune/ND_immune.rds")
HFD <- readRDS("D:/data/23BMI/ND_HFD_MG_snRNAseq/Figure7_immune/HFD_immune.rds")

# Turn counts to CPM
toCPM <- function(mat) {
  lib <- Matrix::colSums(mat)
  cpm <- t( t(mat) / pmax(lib, 1) * 1e6 )
  return(cpm)
}

export_for_ecotyper_discovery <- function(obj, out_dir, assay="RNA") {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  DefaultAssay(obj) <- assay
  expr <- GetAssayData(obj, slot="counts")
  expr <- toCPM(expr)
  
  keep_genes <- Matrix::rowSums(expr > 0) >= 5
  expr <- expr[keep_genes, , drop=FALSE]
  
  colnames(expr) <- gsub("-", "_", colnames(expr))
  colnames(expr) <- gsub(" ", "_", colnames(expr))
  colnames(expr) <- gsub("/", "_", colnames(expr))
  
  expr_mat <- as.matrix(expr)
  
  expr_df <- as.data.frame(expr_mat)
  expr_df <- cbind(Gene = rownames(expr_mat), expr_df)
  rownames(expr_df) <- NULL
  
  expr_file_txt <- file.path(out_dir, "expr_matrix.txt")
  write.table(expr_df, expr_file_txt, sep="\t", quote=FALSE, 
              row.names=FALSE, col.names=TRUE)
  
  system(paste0("gzip -f ", expr_file_txt))
  
  message("Count matrix has finished")
  
  cell_types <- if (!is.null(obj$subcluster)) {
    as.character(obj$subcluster)
  } else {
    as.character(Idents(obj))
  }

  cell_types <- gsub(" ", "_", cell_types)
  cell_types <- gsub("/", "_", cell_types)
  
  cell_ids <- gsub("-", "_", colnames(obj))
  cell_ids <- gsub(" ", "_", cell_ids)
  cell_ids <- gsub("/", "_", cell_ids)
  
  annotation <- data.frame(
    ID       = cell_ids,
    CellType = cell_types,
    Sample   = obj$orig.ident,
    stringsAsFactors = FALSE
  )
  
  anno_file <- file.path(out_dir, "annotation.txt")
  write.table(annotation, anno_file, sep="\t", quote=FALSE, row.names=FALSE)
  
}

export_for_ecotyper_discovery(
  ND,
  out_dir = "D:/data/23BMI/ND_HFD_MG_snRNAseq/Figure7_immune/ND_discovery"
)

export_for_ecotyper_discovery(
  HFD,
  out_dir = "D:/data/23BMI/ND_HFD_MG_snRNAseq/Figure7_immune/HFD_discovery"
)