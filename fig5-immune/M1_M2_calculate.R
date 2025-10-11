library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)


M1_genes <- c("Tnf","Cd86","Ccl2","Ccl4")
M2_genes <- c("Retnla","Clec10a","Fn1")

Immune <- readRDS("D:/data/23BMI/ND_HFD_MG_snRNAseq/Figure7_immune/All_immune.rds")
Idents(Immune) <- "subcluster"
cells_to_keep <- subset(Immune, idents = c("Monocyte/Inf mac","Tissue-resident mac"))

cells_to_keep <- NormalizeData(
  cells_to_keep,
  normalization.method = "LogNormalize",
  scale.factor = 1e4
)

cells_to_keep <- FindVariableFeatures(
  cells_to_keep,
  selection.method = "vst",
  nfeatures = 3000
)

cells_to_keep <- ScaleData(
  cells_to_keep,
  vars.to.regress = c("percent.mt","nCount_RNA"),
  features = rownames(cells_to_keep),
  do.scale = TRUE,
  do.center = TRUE
)

cells_to_keep <- RunPCA(
  cells_to_keep,
  features = VariableFeatures(cells_to_keep),
  npcs = 50,
  ndims.print = 1:5, nfeatures.print = 30,
  reduction.name = "pca", reduction.key = "PC_"
)

cells_to_keep <- FindNeighbors(
  cells_to_keep,
  dims = 1:30,
  k.param = 20,
  nn.method = "annoy",
  annoy.n.trees = 50
)

cells_to_keep <- FindClusters(
  cells_to_keep,
  resolution = 0.2,
  algorithm = 1,
  n.start = 10,
  n.iter = 10
)

cells_to_keep <- RunUMAP(
  cells_to_keep,
  dims = 1:30,
  n.neighbors = 50,
  min.dist = 0.6,
  spread = 1.0,
  metric = "cosine",
  n.components = 2,
  seed.use = 42,
  learning.rate = 1.0,
  repulsion.strength = 1.0,
  local.connectivity = 1
)

# M1/2 score
cells_to_keep <- AddModuleScore(cells_to_keep, features = list(M1_genes), name = "M1_Score")
cells_to_keep <- AddModuleScore(cells_to_keep, features = list(M2_genes), name = "M2_Score")
ND <- subset(cells_to_keep, subset = orig.ident == "ND")
HFD <- subset(cells_to_keep, subset = orig.ident == "HFD")

mean(ND$M1_Score1)
mean(HFD$M1_Score1)
mean(ND$M1_Score1)
mean(HFD$M2_Score1)


Idents(cells_to_keep) <- "orig.ident"
expr <- GetAssayData(cells_to_keep, slot = "data")  # normalized data
gene_list <- c(M1_genes, M2_genes)
gene_list <- gene_list[gene_list %in% rownames(expr)]  

df <- data.frame(t(as.matrix(expr[gene_list, ]))) %>%
  mutate(orig.ident = cells_to_keep$orig.ident)

df_long <- df %>%
  pivot_longer(
    cols = all_of(gene_list),
    names_to = "Gene",
    values_to = "Expression"
  )

df_long$GeneType <- ifelse(df_long$Gene %in% M1_genes, "M1", "M2")

df_long$Gene <- factor(df_long$Gene, levels = gene_list)
df_long$orig.ident <- factor(df_long$orig.ident, levels = c("ND", "HFD"))

df_long$Gene_Condition <- paste0(df_long$Gene, "_", df_long$orig.ident)

y_levels <- as.vector(sapply(gene_list, function(g) paste0(g, "_", c("ND", "HFD"))))
df_long$Gene_Condition <- factor(df_long$Gene_Condition, levels = y_levels)

p <- ggplot(df_long, aes(x = Expression, y = Gene_Condition, fill = GeneType)) +
  geom_violin(scale = "width", trim = TRUE) +
  scale_fill_manual(
    values = c("M1" = "#3b85db", "M2" = "#9db648"),
    labels = c("M1 genes", "M2 genes")
  ) +
  scale_y_discrete(
    labels = function(x) {
      gsub(".*_", "", x)
    }
  ) +
  labs(
    x = "Normalized expression", 
    y = NULL, 
    title = "M1/M2 signature genes in ND vs HFD macrophages",
    fill = "Gene Type"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 10),
    panel.grid.minor = element_blank(),
    legend.position = "top"
  ) +
  facet_grid(Gene ~ ., scales = "free_y", space = "free_y", switch = "y")

p

ggsave("M1_M2_expression.pdf", plot = p, width = 10, height = 18)
