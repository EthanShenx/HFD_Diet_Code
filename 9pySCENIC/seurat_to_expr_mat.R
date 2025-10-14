library(Seurat)
library(dplyr)

seurat_obj <- readRDS("harmony_All_sub_sub.rds")

# 设置 assay
DefaultAssay(seurat_obj) <- "RNA"

# 获取原始 counts 或 normalized counts
# 通常使用原始 counts 或 TPM/CPM
expr_mat <- as.data.frame(GetAssayData(seurat_obj, slot = "counts"))

# 行是基因，列是细胞，但 pySCENIC 要求行是细胞，列是基因
expr_mat_t <- t(expr_mat)

# 行名为细胞条形码，列名为基因
# 保存为 tsv
write.table(expr_mat_t,
            file = "expr_mat_for_pyscenic.tsv",
            sep = "\t", quote = FALSE, col.names = NA)
