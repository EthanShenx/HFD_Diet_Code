library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(AUCell)
library(patchwork)
library(viridis)
library(ggpubr)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(ggridges)
library(rstatix)

# 设置工作目录
setwd("/Users/renyixiang/Desktop/data/23BMI/ND_HFD_MG_snRNAseq/Figure7_immune")

# 读取数据
obj <- readRDS("All_for_CCC.rds")
DefaultAssay(obj) <- "RNA"
Idents(obj) <- "celltype_state"

# ============================================================================
# 1. 定义基因集
# ============================================================================

# Rora Downstream(Rora low --> high ):
Rora_downstream_genes <- c(
  "E2f1", "Cdc6", "E2f8", "E2f7",
  "Il6", "Cxcl1", "Cxcl2", "Cxcl15","Ccl2",
  "Ccl5","Ccl7","Tnf","Il1a","Cxcl10",
  "Nfkb1","Socs3","Ccnd1","Trem1"
)

Rora_downstream_genes <- intersect(Rora_downstream_genes, rownames(obj))
# ============================================================================
# 2. AUCell Scoring
# ============================================================================

gene_sets <- list(Rora_downstream = Rora_downstream_genes)
expr_matrix <- GetAssayData(obj, slot = "data", assay = "RNA")

# 构建基因排名矩阵
cells_rankings <- AUCell_buildRankings(expr_matrix, 
                                       nCores = 4,
                                       plotStats = FALSE)

# 计算 AUC 值
cells_AUC <- AUCell_calcAUC(gene_sets, cells_rankings)
auc_scores <- getAUC(cells_AUC)

# 添加到 Seurat 对象
obj$Rora_downstream_AUC <- auc_scores["Rora_downstream", ]

# 提取上皮细胞亚群
epi_rora <- subset(obj, celltype_state %in% c("Basal", "LumProg", "HormSens"))

diet_colors <- c(
  "ND" = "#3C5488",   
  "HFD" = "#DC0000"   
)

celltype_colors <- c(
  "Basal" = "#E64B35",
  "LumProg" = "#4DBBD5",
  "HormSens" = "#00A087"
)

plot_data <- data.frame(
  AUC = epi_rora$Rora_downstream_AUC,
  Diet = factor(epi_rora$orig.ident, levels = c("ND", "HFD")),
  CellType = epi_rora$celltype_state
)

plot_data_facet <- data.frame(
  AUC = epi_rora$Rora_downstream_AUC,
  CellType = epi_rora$celltype_state,
  Diet = factor(epi_rora$orig.ident, levels = c("ND", "HFD"))
)

stat_test_overall <- compare_means(
  AUC ~ Diet, 
  data = plot_data,
  method = "wilcox.test"
)

stat_test_by_celltype <- plot_data_facet %>%
  group_by(CellType) %>%
  rstatix::wilcox_test(AUC ~ Diet) %>%
  rstatix::add_significance() %>%
  rstatix::add_xy_position(x = "Diet")

summary_data <- plot_data_facet %>%
  group_by(CellType, Diet) %>%
  summarise(
    mean_AUC = mean(AUC),
    se_AUC = sd(AUC) / sqrt(n()),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Diet,
    values_from = c(mean_AUC, se_AUC)
  ) %>%
  mutate(
    fold_change = mean_AUC_HFD / mean_AUC_ND,
    log2FC = log2(fold_change)
  )

significance <- stat_test_by_celltype %>%
  dplyr::select(CellType, p, p.signif)

summary_data <- left_join(summary_data, significance, by = "CellType")

summary_paired <- plot_data_facet %>%
  group_by(CellType, Diet) %>%
  summarise(
    mean_AUC = mean(AUC),
    se_AUC = sd(AUC) / sqrt(n()),
    .groups = "drop"
  )
summary_paired$Diet <- factor(summary_paired$Diet, levels = c("ND", "HFD"))


p4 <- ggplot(plot_data_facet, aes(x = Diet, y = AUC, fill = Diet)) +
  geom_violin(alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA) +
  scale_fill_manual(values = diet_colors) +
  facet_wrap(~CellType, scales = "free_y", nrow = 1) +
  stat_pvalue_manual(
    stat_test_by_celltype,
    label = "p = {p}",
    size = 3.5,
    bracket.size = 0.5,
    tip.length = 0.01
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
    axis.title.x = element_blank(),
    legend.position = "top",
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(size = 11, face = "bold")
  ) +
  labs(
    title = "Pathway Activity by Cell Type: ND vs HFD",
    y = "AUC Score"
  )
p4

p6 <- ggplot(summary_paired, aes(x = Diet, y = mean_AUC, group = CellType, color = CellType)) +
  geom_line(size = 1.2, alpha = 0.8) +
  geom_point(size = 4, alpha = 0.9) +
  geom_errorbar(aes(ymin = mean_AUC - se_AUC, 
                    ymax = mean_AUC + se_AUC),
                width = 0.1, alpha = 0.6) +
  scale_color_manual(values = celltype_colors) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 11, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10)
  ) +
  labs(
    title = "Pathway Activity Change: ND → HFD",
    x = "Diet Condition",
    y = "Mean AUC Score",
    color = "Cell Type"
  )
p6

key_genes <- Rora_downstream_genes[1:min(10, length(Rora_downstream_genes))]
gene_expr <- FetchData(epi_rora, vars = c(key_genes, "celltype_state", "orig.ident", "Rora_downstream_AUC"))


avg_expr <- gene_expr %>%
  group_by(celltype_state, orig.ident) %>%
  summarise(across(all_of(key_genes), mean), .groups = "drop") %>%
  unite("Group", celltype_state, orig.ident, sep = "_")

mat <- as.matrix(avg_expr[, key_genes])
rownames(mat) <- avg_expr$Group
mat <- t(scale(t(mat)))  

group_info <- data.frame(
  Group = avg_expr$Group,
  CellType = sapply(strsplit(avg_expr$Group, "_"), `[`, 1),
  Diet = sapply(strsplit(avg_expr$Group, "_"), `[`, 2)
)

row_anno <- rowAnnotation(
  CellType = group_info$CellType,
  Diet = group_info$Diet,
  col = list(
    CellType = celltype_colors,
    Diet = diet_colors
  ),
  annotation_name_side = "top"
)

ht <- Heatmap(
  mat,
  name = "Z-score",
  col = colorRamp2(c(-2, 0, 2), c("#3B4992", "white", "#EE0000")),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10),
  left_annotation = row_anno,
  column_title = "Gene Expression: ND vs HFD",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  row_title = "Cell Type + Diet",
  row_title_gp = gpar(fontsize = 12, fontface = "bold"),
  heatmap_legend_param = list(
    title = "Expression\n(Z-score)",
    legend_direction = "vertical"
  )
)
ht
