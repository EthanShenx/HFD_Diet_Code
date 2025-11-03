library(Seurat)
library(Matrix)
library(dplyr)

# 读取合并的免疫细胞数据
All_immune <- readRDS("D:/data/23BMI/ND_HFD_MG_snRNAseq/Figure7_immune/All_immune.rds")

# Turn counts to CPM
toCPM <- function(mat) {
  lib <- Matrix::colSums(mat)
  cpm <- t( t(mat) / pmax(lib, 1) * 1e6 )
  return(cpm)
}

# 导出函数（修改版）
export_merged_for_ecotyper_discovery <- function(obj, out_dir, assay="RNA") {
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
  
  message("Count matrix export finished: ", file.path(out_dir, "expr_matrix.txt.gz"))
  
  cell_type <- if (!is.null(obj$subcluster)) {
    as.character(obj$subcluster)
  } else {
    as.character(Idents(obj))
  }
  cell_type <- gsub(" ", "_", cell_type)
  cell_type <- gsub("/", "_", cell_type)
  
  sample_combined <- paste0(obj$orig.ident, "_", cell_type)
  sample_combined <- gsub(" ", "_", sample_combined)
  sample_combined <- gsub("/", "_", sample_combined)
  
  cell_ids <- gsub("-", "_", colnames(obj))
  cell_ids <- gsub(" ", "_", cell_ids)
  cell_ids <- gsub("/", "_", cell_ids)
  
  annotation <- data.frame(
    ID       = cell_ids,
    CellType = cell_type,
    Sample   = sample_combined,
    stringsAsFactors = FALSE
  )
  
  anno_file <- file.path(out_dir, "annotation.txt")
  write.table(annotation, anno_file, sep="\t", quote=FALSE, row.names=FALSE)
  
  message("Annotation export finished: ", anno_file)
  message("\nSummary:")
  message("  - Genes: ", nrow(expr_df) - 1)
  message("  - Cells: ", ncol(expr_df) - 1)
  message("  - Unique CellTypes: ", length(unique(cell_type)))
  message("  - Unique Samples: ", length(unique(sample_combined)))
  
  message("\nCellType composition:")
  print(table(cell_type))
  
  message("\nSample composition:")
  print(table(sample_combined))
}

export_merged_for_ecotyper_discovery(
  All_immune,
  out_dir = "D:/data/23BMI/ND_HFD_MG_snRNAseq/Figure7_immune/ND_HFD_merged_discovery"
)