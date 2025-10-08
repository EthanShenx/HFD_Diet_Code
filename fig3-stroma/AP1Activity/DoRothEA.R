SEURAT_RDS      <- "/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_All_sub_sub.rds"
SUBCLUSTER_COL  <- "subcluster"
MATURE_LABEL    <- "HormSens"
DIET_COL        <- "orig.ident"
SAMPLE_COL      <- NA
SPECIES         <- "mouse"
AP1_TFS         <- c("Jun","Junb","Jund","Fos","Fosb","Fosl1","Fosl2")
F3_GENE         <- "F3"

pkgs <- c("Seurat", "dorothea", "decoupleR", "viper", "tidyverse", "Matrix",
          "lme4", "lmerTest", "ggpubr", "patchwork")
for(p in pkgs){
  if(!requireNamespace(p, quietly = TRUE)){
    message("Installing ", p, " ...")
    if(p %in% c("dorothea","decoupleR")){
      if(!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(p, update = FALSE, ask = FALSE)
    } else {
      install.packages(p, dependencies = TRUE)
    }
  }
}
invisible(lapply(pkgs, library, character.only = TRUE))

seu <- readRDS(SEURAT_RDS)
DefaultAssay(seu) <- "RNA"

regulon_all <- dorothea_hs
reg_ap1 <- regulon_all %>%
  dplyr::filter(confidence %in% c("A"),
                tf %in% AP1_TFS)

emat <- Seurat::GetAssayData(seu, assay = "RNA", slot = "data")

genes_keep <- intersect(rownames(emat), unique(reg_ap1$target))
emat_sub   <- emat[genes_keep, , drop = FALSE]      

emat_dense <- as.matrix(emat_sub)

tf_acts <- dorothea::run_viper(
  emat_dense,
  reg_ap1,
  options = list(minsize = 5, eset.filter = FALSE)
)

tfs_present <- intersect(AP1_TFS, rownames(tf_acts))
seu$AP1_activity <- colMeans(tf_acts[tfs_present, , drop = FALSE])
seu$F3_expr      <- as.numeric(seu[["RNA"]]@data[F3_GENE, ])

################################################
################################################
################################################

AP1_TFS_UP <- toupper(AP1_TFS)
rownames(tf_acts) <- toupper(rownames(tf_acts))
tfs_present <- intersect(AP1_TFS_UP, rownames(tf_acts))
seu$AP1_activity <- colMeans(tf_acts[tfs_present, , drop = FALSE])

stopifnot(F3_GENE %in% rownames(seu[["RNA"]]@data))
seu$F3_expr <- as.numeric(seu[["RNA"]]@data[F3_GENE, ])

seu$sample_id <- as.character(seu@meta.data[[DIET_COL]])

seu$diet_group <- ifelse(grepl("HFD", seu$sample_id, ignore.case = TRUE), "HFD", "ND")
seu$diet_group <- factor(seu$diet_group, levels = c("ND","HFD"))

Idents(seu) <- SUBCLUSTER_COL
mature <- subset(seu, idents = MATURE_LABEL)

meta_m <- mature@meta.data

wilcox_ap1 <- wilcox.test(AP1_activity ~ diet_group, data = meta_m)
wilcox_f3  <- wilcox.test(F3_expr      ~ diet_group, data = meta_m)

corr_cell <- suppressWarnings(cor.test(meta_m$AP1_activity, meta_m$F3_expr, method = "spearman"))

library(dplyr); library(ggplot2)
pb <- meta_m %>%
  group_by(sample_id, diet_group) %>%
  summarise(n_cells = n(),
            AP1_activity = mean(AP1_activity, na.rm = TRUE),
            F3_expr      = mean(F3_expr,      na.rm = TRUE),
            .groups = "drop")

corr_pb <- if(nrow(pb) >= 3) suppressWarnings(cor.test(pb$AP1_activity, pb$F3_expr, method = "spearman")) else NULL

p1 <- ggplot(meta_m, aes(x = diet_group, y = AP1_activity)) +
  geom_violin(trim = TRUE) + geom_boxplot(width = 0.15, outlier.shape = NA) +
  labs(title = "AP-1 activity", x = "Diet", y = "AP1_activity (VIPER)") +
  theme_classic()

p3 <- ggplot(meta_m, aes(x = AP1_activity, y = F3_expr, color = diet_group)) +
  geom_point(alpha = 0.35, size = 0.6) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(title = "AP-1 activity vs F3 expression",
       x = "AP1_activity (VIPER)", y = paste0(F3_GENE, " (log-normalized)")) +
  theme_classic()


cat("\n===== 统计结果（HormSens）=====\n")
cat("AP-1活性：HFD vs CTRL (Wilcoxon)  p =", signif(wilcox_ap1$p.value, 3), "\n")
cat("F3表达 ：HFD vs CTRL (Wilcoxon)  p =", signif(wilcox_f3$p.value, 3),  "\n")
cat("细胞层相关：AP-1活性 vs F3表达 (Spearman) rho =", signif(corr_cell$estimate,3),
    ", p =", signif(corr_cell$p.value,3), "\n")
if(!is.null(corr_pb)){
  cat("样本层相关：AP-1活性 vs F3表达 (Spearman) rho =", signif(corr_pb$estimate,3),
      ", p =", signif(corr_pb$p.value,3), "\n")
} else {
  cat("样本层相关：样本数不足（<3），跳过。\n")
}

nfkb_tfs <- c("Rela","Nfkb1","Nfkb2","Relb","Rel")
reg_nfkb <- (if (tolower(SPECIES)=="mouse") dorothea_mm else dorothea_hs) %>%
  filter(confidence %in% c("A","B","C"), tf %in% nfkb_tfs)
genes_keep2 <- intersect(rownames(seu[["RNA"]]@data), unique(reg_nfkb$target))
emat2 <- as.matrix(seu[["RNA"]]@data[genes_keep2, , drop=FALSE])
nfkb_acts <- dorothea::run_viper(emat2, reg_nfkb, options = list(minsize=5, eset.filter=FALSE))
nfkb_present <- intersect(nfkb_tfs, rownames(nfkb_acts))
seu$NFKB_activity <- colMeans(nfkb_acts[nfkb_present, , drop=FALSE])
mature$NFKB_activity <- seu$NFKB_activity[colnames(mature)]

cor_nfkb <- suppressWarnings(cor.test(mature$NFKB_activity, mature$F3_expr, method="spearman"))
cor_nfkb
