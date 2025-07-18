### PREPARATION ####

setwd("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/1composition")

library(Seurat)
library(dplyr)
library(tidyverse)
library(ggalluvial)
library(ggplot2)
library(RColorBrewer)

All <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_all.rds")

fixed_colors <- brewer.pal(7, "Paired")

### SORT SC METADATA FOR PLOTTING ###

All_info <- as.data.frame(All@meta.data)

type_dotplot_df <- All_info %>%
  group_by(orig.ident, cell_type) %>%
  tally() %>%
  group_by(orig.ident) %>%
  mutate(percentage = n / sum(n))

group_dotplot_df_all <- type_dotplot_df

group_dotplot_df_all$orig.ident <- factor(group_dotplot_df_all$orig.ident, levels = c("ND", "HFD"))

### CENTRAL DOT PLOT ###

central <- ggplot(group_dotplot_df_all, aes(x = orig.ident, y = cell_type, size = percentage, color = cell_type)) +
  geom_point(alpha = 0.85) +
  scale_color_manual(values = fixed_colors) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  ) +
  labs(y = "Cell proportion", x = NULL, size = "Cell proportion", color = "Clusters")

central

### SIGNIFICANT? ###

chisq_results <- group_dotplot_df_all %>%
  select(orig.ident, cell_type, n) %>%
  tidyr::pivot_wider(names_from = orig.ident, values_from = n, values_fill = 0) %>%
  rowwise() %>%
  mutate(
    chisq_p = chisq.test(matrix(c(ND, HFD), nrow = 2))$p.value
  ) %>%
  ungroup() %>%
  mutate(
    p_adj = p.adjust(chisq_p, method = "fdr"),
    label = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      TRUE ~ ""
    )
  ) %>%
  as.data.frame()

### UP BARPLOT ###

barplot_x <- ggplot(group_dotplot_df_all, aes(x = orig.ident, y = n, fill = cell_type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = fixed_colors) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "none"
  ) +
  labs(y = "Number of cells", x = NULL)

barplot_x

### RIGHT BARPLOT ###

celltype_sum <- group_dotplot_df_all %>%
  group_by(cell_type) %>%
  summarise(total_n = sum(n), .groups = "drop")

barplot_y <- ggplot(celltype_sum, aes(x = cell_type, y = total_n, fill = cell_type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = fixed_colors) +
  coord_flip() +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "right"
  ) +
  labs(y = "Number of cells", x = NULL)

barplot_y