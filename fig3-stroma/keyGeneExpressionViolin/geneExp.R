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

Areg <- ggplot(df, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(scale = "width", trim = F, adjust = 1, width = 0.8, color = NA) +
  geom_boxplot(width = 0.12, fill = "white", alpha = 1,
               outlier.size = 0.8, outlier.color = "black",
               outlier.shape = 19, linewidth = 0.3) +
  scale_fill_manual(values = condition_colors) +
  theme_classic(base_family = "Arial") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        xis.text.x  = element_text(size = 10, colour = "black"),
        xis.text.y  = element_text(size = 10, colour = "black"),
        legend.position = "none") +
  labs(y = y_lab) +
  ggpubr::stat_compare_means(comparisons = list(c("ND","HFD")),
                             method = "wilcox.test", label = "p.format", size = 3, vjust = -0.5)

Areg

#===== Sp1 HormSens =====

df <- combine_diet_df(cell_type = "HormSens", gene = "Sp1", slot = slot_used)

df <- na.omit(df)
wt <- wilcox.test(Expression ~ Condition, data = df, exact = FALSE)
wt$p.value

condition_colors <- c("ND" = "#74c5be", "HFD" = "#e95503")

y_lab <- "Expression"

Areg <- ggplot(df, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(scale = "width", trim = F, adjust = 1, width = 0.8, color = NA) +
  geom_boxplot(width = 0.12, fill = "white", alpha = 1,
               outlier.size = 0.8, outlier.color = "black",
               outlier.shape = 19, linewidth = 0.3) +
  scale_fill_manual(values = condition_colors) +
  theme_classic(base_family = "Arial") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        xis.text.x  = element_text(size = 10, colour = "black"),
        xis.text.y  = element_text(size = 10, colour = "black"),
        legend.position = "none") +
  labs(y = y_lab) +
  ggpubr::stat_compare_means(comparisons = list(c("ND","HFD")),
                             method = "wilcox.test", label = "p.format", size = 3, vjust = -0.5)

Areg

#===== F3 HormSens =====

# 如未安装：
# install.packages("gghalves")

library(ggplot2)
library(ggpubr)
library(gghalves)

df <- combine_diet_df(cell_type = "HormSens", gene = "F3", slot = slot_used)
df <- na.omit(df)

wt <- wilcox.test(Expression ~ Condition, data = df, exact = FALSE)
wt$p.value

condition_colors <- c("ND" = "#74c5be", "HFD" = "#e95503")
y_lab <- "Expression"

F3 <- ggplot(df, aes(x = Condition, y = Expression)) +

  gghalves::geom_half_violin(aes(fill = Condition),
                             side = "r",
                             trim = FALSE,
                             adjust = 1,
                             width = 0.8,
                             color = NA) +

  gghalves::geom_half_point(side = "l",
                            range_scale = 0.5,
                            position = position_jitter(width = 0.065, height = 0),
                            color = "black",
                            alpha = 0.9,
                            size = 0.4) +

  geom_boxplot(width = 0.12,
               fill = "white",
               alpha = 1,
               outlier.shape = NA,
               linewidth = 0.3) +
  scale_fill_manual(values = condition_colors) +
  theme_classic(base_family = "Arial") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.x  = element_text(size = 10, colour = "black"),
        axis.text.y  = element_text(size = 10, colour = "black"),
        legend.position = "none") +
  labs(y = y_lab) +
  ggpubr::stat_compare_means(comparisons = list(c("ND","HFD")),
                             method = "wilcox.test",
                             label = "p.format",
                             size = 3,
                             vjust = -0.5)

F3

#===== test HormSens =====

# 如未安装：
# install.packages("gghalves")

library(ggplot2)
library(ggpubr)
library(gghalves)

df <- combine_diet_df(cell_type = "HormSens", gene = "Areg", slot = slot_used)
df <- na.omit(df)

wt <- wilcox.test(Expression ~ Condition, data = df, exact = FALSE)
wt$p.value

condition_colors <- c("ND" = "#74c5be", "HFD" = "#e95503")
y_lab <- "Expression"

Test <- ggplot(df, aes(x = Condition, y = Expression)) +

  gghalves::geom_half_violin(aes(fill = Condition),
                             side = "r",
                             trim = FALSE,
                             adjust = 1,
                             width = 0.8,
                             color = NA) +

  gghalves::geom_half_point(side = "l",
                            range_scale = 0.5,
                            position = position_jitter(width = 0.065, height = 0),
                            color = "black",
                            alpha = 0.9,
                            size = 0.4) +

  geom_boxplot(width = 0.12,
               fill = "white",
               alpha = 1,
               outlier.shape = NA,
               linewidth = 0.3) +
  scale_fill_manual(values = condition_colors) +
  theme_classic(base_family = "Arial") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.x  = element_text(size = 10, colour = "black"),
        axis.text.y  = element_text(size = 10, colour = "black"),
        legend.position = "none") +
  labs(y = y_lab) +
  ggpubr::stat_compare_means(comparisons = list(c("ND","HFD")),
                             method = "wilcox.test",
                             label = "p.format",
                             size = 3,
                             vjust = -0.5)

Test

#===== Esr1 HormSens =====

df <- combine_diet_df(cell_type = "HormSens", gene = "Esr1", slot = slot_used)

df <- na.omit(df)
wt <- wilcox.test(Expression ~ Condition, data = df, exact = FALSE)
wt$p.value

condition_colors <- c("ND" = "#74c5be", "HFD" = "#e95503")

y_lab <- "Expression"

Esr1 <- ggplot(df, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(scale = "width", trim = F, adjust = 1, width = 0.8, color = NA) +
  geom_boxplot(width = 0.12, fill = "white", alpha = 1,
               outlier.size = 0.8, outlier.color = "black",
               outlier.shape = 19, linewidth = 0.3) +
  scale_fill_manual(values = condition_colors) +
  theme_classic(base_family = "Arial") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        xis.text.x  = element_text(size = 10, colour = "black"),
        xis.text.y  = element_text(size = 10, colour = "black"),
        legend.position = "none") +
  labs(y = y_lab) +
  ggpubr::stat_compare_means(comparisons = list(c("ND","HFD")),
                             method = "wilcox.test", label = "p.format", size = 3, vjust = -0.5) 

Esr1

#===== Egfr Stroma =====

df <- combine_diet_df(cell_type = "Stroma", gene = "Egfr", slot = slot_used)

df <- na.omit(df)
wt <- wilcox.test(Expression ~ Condition, data = df, exact = FALSE)
wt$p.value

condition_colors <- c("ND" = "#74c5be", "HFD" = "#e95503")

y_lab <- "Expression"

Egfr <- ggplot(df, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(scale = "width", trim = F, adjust = 1, width = 0.8, color = NA) +
  geom_boxplot(width = 0.12, fill = "white", alpha = 1,
               outlier.size = 0.8, outlier.color = "black",
               outlier.shape = 19, linewidth = 0.3) +
  scale_fill_manual(values = condition_colors) +
  theme_classic(base_family = "Arial") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        xis.text.x  = element_text(size = 10, colour = "black"),
        xis.text.y  = element_text(size = 10, colour = "black"),
        legend.position = "none") +
  labs(y = y_lab) +
  ggpubr::stat_compare_means(comparisons = list(c("ND","HFD")),
                             method = "wilcox.test", label = "p.format", size = 3, vjust = -0.5) 

Egfr

# ===== Adam17 HormSens =====
df <- combine_diet_df(cell_type = "HormSens", gene = "Adam17", slot = slot_used)

df <- na.omit(df)
wt <- wilcox.test(Expression ~ Condition, data = df, exact = FALSE)
wt$p.value

condition_colors <- c("ND" = "#74c5be", "HFD" = "#e95503")

y_lab <- "Expression"

Adam17 <- ggplot(df, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(scale = "width", trim = F, adjust = 1, width = 0.8, color = NA) +
  geom_boxplot(width = 0.12, fill = "white", alpha = 1,
               outlier.size = 0.8, outlier.color = "black",
               outlier.shape = 19, linewidth = 0.3) +
  scale_fill_manual(values = condition_colors) +
  theme_classic(base_family = "Arial") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        xis.text.x  = element_text(size = 10, colour = "black"),
        xis.text.y  = element_text(size = 10, colour = "black"),
        legend.position = "none") +
  labs(y = y_lab) +
  ggpubr::stat_compare_means(comparisons = list(c("ND","HFD")),
                             method = "wilcox.test", label = "p.format", size = 3, vjust = -0.5) 

Adam17

# ===== Timp3 HormSens =====
df <- combine_diet_df(cell_type = "HormSens", gene = "Timp3", slot = slot_used)
df <- na.omit(df)
wt <- wilcox.test(Expression ~ Condition, data = df, exact = FALSE)
wt$p.value
condition_colors <- c("ND" = "#74c5be", "HFD" = "#e95503")
y_lab <- "Expression"
Timp3 <- ggplot(df, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(scale = "width", trim = F, adjust = 1, width = 0.8, color = NA) +
  geom_boxplot(width = 0.12, fill = "white", alpha = 1,
               outlier.size = 0.8, outlier.color = "black",
               outlier.shape = 19, linewidth = 0.3) +
  scale_fill_manual(values = condition_colors) +
  theme_classic(base_family = "Arial") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        xis.text.x  = element_text(size = 10, colour = "black"),
        xis.text.y  = element_text(size = 10, colour = "black"),
        legend.position = "none") +
  labs(y = y_lab) +
  ggpubr::stat_compare_means(comparisons = list(c("ND","HFD")),
                             method = "wilcox.test", label = "p.format", size = 3, vjust = -0.5) 
Timp3

# ===== Cited1 HormSens =====
df <- combine_diet_df(cell_type = "HormSens", gene = "Cited1", slot = slot_used)
df <- na.omit(df)
wt <- wilcox.test(Expression ~ Condition, data = df, exact = FALSE)
wt$p.value
condition_colors <- c("ND" = "#74c5be", "HFD" = "#e95503")
y_lab <- "Expression"
Cited1 <- ggplot(df, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(scale = "width", trim = F, adjust = 1, width = 0.8, color = NA) +
  geom_boxplot(width = 0.12, fill = "white", alpha = 1,
               outlier.size = 0.8, outlier.color = "black",
               outlier.shape = 19, linewidth = 0.3) +
  scale_fill_manual(values = condition_colors) +
  theme_classic(base_family = "Arial") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        xis.text.x  = element_text(size = 10, colour = "black"),
        xis.text.y  = element_text(size = 10, colour = "black"),
        legend.position = "none") +
  labs(y = y_lab) +
  ggpubr::stat_compare_means(comparisons = list(c("ND","HFD")),
                             method = "wilcox.test", label = "p.format", size = 3, vjust = -0.5) 
Cited1

# ===== Ghr Stroma =====
df <- combine_diet_df(cell_type = "Stroma", gene = "Ghr", slot = slot_used)
df <- na.omit(df)
wt <- wilcox.test(Expression ~ Condition, data = df, exact = FALSE)
wt$p.value
condition_colors <- c("ND" = "#74c5be", "HFD" = "#e95503")
y_lab <- "Expression"
Ghr <- ggplot(df, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(scale = "width", trim = F, adjust = 1, width = 0.8, color = NA) +
  geom_boxplot(width = 0.12, fill = "white", alpha = 1,
               outlier.size = 0.8, outlier.color = "black",
               outlier.shape = 19, linewidth = 0.3) +
  scale_fill_manual(values = condition_colors) +
  theme_classic(base_family = "Arial") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        xis.text.x  = element_text(size = 10, colour = "black"),
        xis.text.y  = element_text(size = 10, colour = "black"),
        legend.position = "none") +
  labs(y = y_lab) +
  ggpubr::stat_compare_means(comparisons = list(c("ND","HFD")),
                             method = "wilcox.test", label = "p.format", size = 3, vjust = -0.5) 
Ghr

# ===== Ptch1 HormSens =====
df <- combine_diet_df(cell_type = "HormSens", gene = "Ptch1", slot = slot_used)
df <- na.omit(df)
wt <- wilcox.test(Expression ~ Condition, data = df, exact = FALSE)
wt$p.value
condition_colors <- c("ND" = "#74c5be", "HFD" = "#e95503")
y_lab <- "Expression"
Ptch1 <- ggplot(df, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(scale = "width", trim = F, adjust = 1, width = 0.8, color = NA) +
  geom_boxplot(width = 0.12, fill = "white", alpha = 1,
               outlier.size = 0.8, outlier.color = "black",
               outlier.shape = 19, linewidth = 0.3) +
  scale_fill_manual(values = condition_colors) +
  theme_classic(base_family = "Arial") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        xis.text.x  = element_text(size = 10, colour = "black"),
        xis.text.y  = element_text(size = 10, colour = "black"),
        legend.position = "none") +
  labs(y = y_lab) +
  ggpubr::stat_compare_means(comparisons = list(c("ND","HFD")),
                             method = "wilcox.test", label = "p.format", size = 3, vjust = -0.5) 
Ptch1

# ===== Gli2 HormSens =====
df <- combine_diet_df(cell_type = "HormSens", gene = "Gli2", slot = slot_used)
df <- na.omit(df)
wt <- wilcox.test(Expression ~ Condition, data = df, exact = FALSE)
wt$p.value
condition_colors <- c("ND" = "#74c5be", "HFD" = "#e95503")
y_lab <- "Expression"
Gli2 <- ggplot(df, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(scale = "width", trim = F, adjust = 1, width = 0.8, color = NA) +
  geom_boxplot(width = 0.12, fill = "white", alpha = 1,
               outlier.size = 0.8, outlier.color = "black",
               outlier.shape = 19, linewidth = 0.3) +
  scale_fill_manual(values = condition_colors) +
  theme_classic(base_family = "Arial") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        xis.text.x  = element_text(size = 10, colour = "black"),
        xis.text.y  = element_text(size = 10, colour = "black"),
        legend.position = "none") +
  labs(y = y_lab) +
  ggpubr::stat_compare_means(comparisons = list(c("ND","HFD")),
                             method = "wilcox.test", label = "p.format", size = 3, vjust = -0.5) 
Gli2

# ===== Tnfrsf21 HormSens =====
df <- combine_diet_df(cell_type = "HormSens", gene = "Tnfrsf21", slot = slot_used)
df <- na.omit(df)
wt <- wilcox.test(Expression ~ Condition, data = df, exact = FALSE)
wt$p.value
condition_colors <- c("ND" = "#74c5be", "HFD" = "#e95503")
y_lab <- "Expression"
Tnfrsf21 <- ggplot(df, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(scale = "width", trim = F, adjust = 1, width = 0.8, color = NA) +
  geom_boxplot(width = 0.12, fill = "white", alpha = 1,
               outlier.size = 0.8, outlier.color = "black",
               outlier.shape = 19, linewidth = 0.3) +
  scale_fill_manual(values = condition_colors) +
  theme_classic(base_family = "Arial") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        xis.text.x  = element_text(size = 10, colour = "black"),
        xis.text.y  = element_text(size = 10, colour = "black"),
        legend.position = "none") +
  labs(y = y_lab) +
  ggpubr::stat_compare_means(comparisons = list(c("ND","HFD")),
                             method = "wilcox.test", label = "p.format", size = 3, vjust = -0.5) 
Tnfrsf21

# # ===== Tgfb3 Stroma =====
# df <- combine_diet_df(cell_type = "Stroma", gene = "Tgfb3", slot = slot_used)
# df <- na.omit(df)
# wt <- wilcox.test(Expression ~ Condition, data = df, exact = FALSE)
# wt$p.value
# condition_colors <- c("ND" = "#74c5be", "HFD" = "#e95503")
# y_lab <- "Expression"
# Tgfb3 <- ggplot(df, aes(x = Condition, y = Expression, fill = Condition)) +
#   geom_violin(scale = "width", trim = F, adjust = 1, width = 0.8, color = NA) +
#   geom_boxplot(width = 0.12, fill = "white", alpha = 1,
#                outlier.size = 0.8, outlier.color = "black",
#                outlier.shape = 19, linewidth = 0.3) +
#   scale_fill_manual(values = condition_colors) +
#   theme_classic(base_family = "Arial") +
#   theme(axis.title.x = element_blank(),
#         axis.title.y = element_text(size = 12),
#         xis.text.x  = element_text(size = 10, colour = "black"),
#         xis.text.y  = element_text(size = 10, colour = "black"),
#         legend.position = "none") +
#   labs(y = y_lab) +
#   ggpubr::stat_compare_means(comparisons = list(c("ND","HFD")),
#                              method = "wilcox.test", label = "p.format", size = 3, vjust = -0.5) 
# Tgfb3
# 
# ===== Tgfb1 Stroma =====
df <- combine_diet_df(cell_type = "HormSens", gene = "Tgfb1", slot = slot_used)
df <- na.omit(df)
wt <- wilcox.test(Expression ~ Condition, data = df, exact = FALSE)
wt$p.value
condition_colors <- c("ND" = "#74c5be", "HFD" = "#e95503")
y_lab <- "Expression"
Tgfb2 <- ggplot(df, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_violin(scale = "width", trim = F, adjust = 1, width = 0.8, color = NA) +
  geom_boxplot(width = 0.12, fill = "white", alpha = 1,
               outlier.size = 0.8, outlier.color = "black",
               outlier.shape = 19, linewidth = 0.3) +
  scale_fill_manual(values = condition_colors) +
  theme_classic(base_family = "Arial") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        xis.text.x  = element_text(size = 10, colour = "black"),
        xis.text.y  = element_text(size = 10, colour = "black"),
        legend.position = "none") +
  labs(y = y_lab) +
  ggpubr::stat_compare_means(comparisons = list(c("ND","HFD")),
                             method = "wilcox.test", label = "p.format", size = 3, vjust = -0.5) 
Tgfb2

# ==== patch together ====
library(patchwork)
Areg | Esr1 | Adam17 | Timp3 | Cited1 | Egfr | Ghr | Ptch1 | Gli2
Tgfb2 | Tgfb3
