library(Seurat)
library(rlang)
library(dplyr)
library(tidyr)
library(ggdist)
library(ggpubr)
library(ggsignif)
library(extrafont)
loadfonts(quiet = TRUE)
base_family <- "Arial"
library(ggplot2)
library(scales)

All_diet <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_All_sub.rds")

combine_diet_df <- function(cell_type, gene, slot = "data") {
  cells_diet <- Cells(All_diet)[All_diet$cell_type == cell_type]
  sub <- subset(All_diet, cells = cells_diet)
  DefaultAssay(sub) <- DefaultAssay(All_diet)

  sub$orig.ident <- dplyr::recode(toupper(sub$orig.ident), "ND"="ND", "HFD"="HFD")
  sub$orig.ident <- factor(sub$orig.ident, levels = c("ND", "HFD"))

  mat <- GetAssayData(sub, slot = slot)
  if (!gene %in% rownames(mat)) stop(sprintf("Gene '%s' not found in slot '%s'", gene, slot))
  expr <- as.numeric(mat[gene, ])

  df <- data.frame(
    Expression = expr,
    Condition  = sub$orig.ident
  )
  df$Condition <- factor(df$Condition, levels = c("ND", "HFD"))
  df
}

slot_used <- "data"

#===== Areg HormSens =====

df <- combine_diet_df(cell_type = "HormSens", gene = "Areg", slot = slot_used)

df <- na.omit(df)
wt <- wilcox.test(Expression ~ Condition, data = df, exact = FALSE)
wt$p.value

condition_colors <- c("ND" = "#74c5be", "HFD" = "#e95503")

y_lab <- "Expression"

p <- ggplot(df, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(scale = "width", trim = TRUE, adjust = 1, width = 0.8, color = NA) +
  geom_boxplot(width = 0.12, fill = "white", alpha = 1,
               outlier.size = 0.8, outlier.color = "black",
               outlier.shape = 19, linewidth = 0.3) +
  scale_fill_manual(values = condition_colors) +
  theme_classic(base_family = "Arial") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.x  = element_text(size = 10),
        axis.text.y  = element_text(size = 10),
        legend.position = "none") +
  labs(y = y_lab) +
  ggpubr::stat_compare_means(comparisons = list(c("ND","HFD")),
                             method = "wilcox.test", label = "p.format")

p

#===== Esr-alpha HormSens =====

df <- combine_diet_df(cell_type = "HormSens", gene = "Esr1", slot = slot_used)

df <- na.omit(df)
wt <- wilcox.test(Expression ~ Condition, data = df, exact = FALSE)
wt$p.value

condition_colors <- c("ND" = "#74c5be", "HFD" = "#e95503")

y_lab <- "Expression"

p <- ggplot(df, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(scale = "width", trim = TRUE, adjust = 1, width = 0.8, color = NA) +
  geom_boxplot(width = 0.12, fill = "white", alpha = 1,
               outlier.size = 0.8, outlier.color = "black",
               outlier.shape = 19, linewidth = 0.3) +
  scale_fill_manual(values = condition_colors) +
  theme_classic(base_family = "Arial") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.x  = element_text(size = 10),
        axis.text.y  = element_text(size = 10),
        legend.position = "none") +
  labs(y = y_lab) +
  ggpubr::stat_compare_means(comparisons = list(c("ND","HFD")),
                             method = "wilcox.test", label = "p.format")

p

#===== Egfr Stroma =====

df <- combine_diet_df(cell_type = "Stroma", gene = "Egfr", slot = slot_used)

df <- na.omit(df)
wt <- wilcox.test(Expression ~ Condition, data = df, exact = FALSE)
wt$p.value

condition_colors <- c("ND" = "#74c5be", "HFD" = "#e95503")

y_lab <- "Expression"

p <- ggplot(df, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(scale = "width", trim = TRUE, adjust = 1, width = 0.8, color = NA) +
  geom_boxplot(width = 0.12, fill = "white", alpha = 1,
               outlier.size = 0.8, outlier.color = "black",
               outlier.shape = 19, linewidth = 0.3) +
  scale_fill_manual(values = condition_colors) +
  theme_classic(base_family = "Arial") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.x  = element_text(size = 10),
        axis.text.y  = element_text(size = 10),
        legend.position = "none") +
  labs(y = y_lab) +
  ggpubr::stat_compare_means(comparisons = list(c("ND","HFD")),
                             method = "wilcox.test", label = "p.format")

p

# ===== Adam17 HormSens =====
df <- combine_diet_df(cell_type = "HormSens", gene = "Adam17", slot = slot_used)

df <- na.omit(df)
wt <- wilcox.test(Expression ~ Condition, data = df, exact = FALSE)
wt$p.value

condition_colors <- c("ND" = "#74c5be", "HFD" = "#e95503")

y_lab <- "Expression"

p <- ggplot(df, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(scale = "width", trim = TRUE, adjust = 1, width = 0.8, color = NA) +
  geom_boxplot(width = 0.12, fill = "white", alpha = 1,
               outlier.size = 0.8, outlier.color = "black",
               outlier.shape = 19, linewidth = 0.3) +
  scale_fill_manual(values = condition_colors) +
  theme_classic(base_family = "Arial") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.x  = element_text(size = 10),
        axis.text.y  = element_text(size = 10),
        legend.position = "none") +
  labs(y = y_lab) +
  ggpubr::stat_compare_means(comparisons = list(c("ND","HFD")),
                             method = "wilcox.test", label = "p.format")

p

