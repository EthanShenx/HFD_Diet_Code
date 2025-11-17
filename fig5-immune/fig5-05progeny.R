# Load required libraries
library(Seurat)
library(progeny)
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(ggpubr)

# Read and prepare data
obj_path <- "/Users/renyixiang/Desktop/data/23BMI/ND_HFD_MG_snRNAseq/Figure7_immune/All_for_CCC.rds"
obj <- readRDS(obj_path)
subcol <- "celltype_state"
obj$combo <- paste0(obj$orig.ident, "_", obj@meta.data[[subcol]])
Idents(obj) <- "celltype_state"
obj <- subset(obj, idents = c("MacS01", "MacS02", "MacS03", "MacS04"))

# Run Progeny analysis
expr_matrix <- as.matrix(GetAssayData(obj, slot = "data", assay = "RNA"))
pathway_scores <- progeny(expr_matrix, scale = TRUE, organism = "Mouse", top = 500)

# Add pathway scores to metadata
pathway_df <- as.data.frame(pathway_scores)
pathway_df <- pathway_df[colnames(obj), , drop = FALSE]
for (pathway in colnames(pathway_df)) {
  obj@meta.data[[pathway]] <- pathway_df[[pathway]]
}

# Prepare data for plotting
pathways <- colnames(pathway_df)
combo_order <- as.vector(rbind(paste0("ND_MacS0", 1:4),
                               paste0("HFD_MacS0", 1:4)))


plot_data <- obj@meta.data %>%
  dplyr::select(combo, all_of(pathways)) %>%
  mutate(
    Condition = ifelse(grepl("^ND_", combo), "ND", "HFD"),
    CellType = gsub("^(ND|HFD)_", "", combo)
  ) %>%
  pivot_longer(cols = all_of(pathways), 
               names_to = "Pathway", 
               values_to = "Score")

# Calculate mean scores for bar plot
plot_data_mean <- plot_data %>%
  group_by(combo, Pathway, Condition, CellType) %>%
  summarise(Score = mean(Score, na.rm = TRUE), .groups = "drop") %>%
  mutate(combo = factor(combo, levels = combo_order))

# Prepare comparisons for each cell type
comps_list <- lapply(c("MacS01", "MacS02", "MacS03", "MacS04"), function(s) {
  c(paste0("ND_", s), paste0("HFD_", s))
})

# Create color palette
colors_paired <- brewer.pal(8, "Paired")
color_mapping <- setNames(colors_paired[1:8], combo_order)

# Create the base plot with mean values
p <- ggplot(plot_data_mean, aes(x = combo, y = Score, fill = combo)) +
  geom_bar(stat = "identity", width = 0.7) +
  facet_wrap(~ Pathway, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = color_mapping) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    axis.title = element_text(size = 11, face = "bold"),
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    panel.grid.minor = element_blank()
  ) +
  labs(x = "Cell Type & Condition", 
       y = "Pathway Activity Score",
       fill = "Group",
       title = "Progeny Pathway Activity Analysis",
       subtitle = "Comparison between ND and HFD conditions across MacS subtypes")

# Add significance using original cell-level data
p <- p + 
  ggpubr::stat_compare_means(
    data = plot_data %>% mutate(combo = factor(combo, levels = combo_order)),
    aes(x = combo, y = Score),
    comparisons   = comps_list,
    method        = "wilcox.test",
    label         = "p.signif",
    hide.ns       = TRUE,
    bracket.size  = 0.25,
    size          = 3,
    tip.length    = 0.01,
    step.increase = 0.08
  )
p
# Save plots
ggsave("Progeny_Pathway_Analysis.pdf", plot = p, width = 16, height = 12, dpi = 300)
ggsave("Progeny_Pathway_Analysis.png", plot = p, width = 16, height = 12, dpi = 300)

# Print plot
print(p)