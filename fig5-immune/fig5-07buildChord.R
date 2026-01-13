library(CellChat)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

setwd("/Users/renyixiang/Desktop/data/23BMI/ND_HFD_MG_snRNAseq/Figure7_immune/CellChat")

# 读取数据
HFD <- readRDS("cellchat_HFD_mac_state.rds")
ND <- readRDS("cellchat_ND_mac_state.rds")

# 定义来源和靶标
mac_sources <- c("MacS01", "MacS03")
epithelial_targets <- c("Basal", "LumProg", "HormSens")

# 分别提取ND和HFD的通讯
net_ND <- subsetCommunication(ND, 
                              sources.use = mac_sources,
                              targets.use = epithelial_targets)
net_HFD <- subsetCommunication(HFD,
                               sources.use = mac_sources,
                               targets.use = epithelial_targets)

# 创建LR pair标识
net_ND$LR_pair <- paste(net_ND$ligand, net_ND$receptor, sep = "-")
net_HFD$LR_pair <- paste(net_HFD$ligand, net_HFD$receptor, sep = "-")

# 合并并计算差异
diff_LR <- full_join(
  net_ND %>% select(source, target, LR_pair, ligand, receptor, pathway_name, prob, pval, interaction_name),
  net_HFD %>% select(source, target, LR_pair, ligand, receptor, pathway_name, prob, pval, interaction_name),
  by = c("source", "target", "LR_pair", "ligand", "receptor", "pathway_name", "interaction_name"),
  suffix = c("_ND", "_HFD")
)

# 填充NA
diff_LR$prob_ND[is.na(diff_LR$prob_ND)] <- 0
diff_LR$prob_HFD[is.na(diff_LR$prob_HFD)] <- 0
diff_LR$pval_ND[is.na(diff_LR$pval_ND)] <- 1
diff_LR$pval_HFD[is.na(diff_LR$pval_HFD)] <- 1

# 计算差异指标
diff_LR <- diff_LR %>%
  mutate(
    prob_diff = prob_HFD - prob_ND,
    log2FC = log2((prob_HFD) / (prob_ND)),
    pct_change = (prob_HFD - prob_ND) / (prob_ND) * 100,
    max_prob = pmax(prob_ND, prob_HFD),
    significant = (pval_ND < 0.05) | (pval_HFD < 0.05)
  )

# 筛选有意义的变化
diff_LR_sig <- diff_LR %>%
  filter(abs(prob_diff) > 0.03) %>%
  arrange(desc(abs(prob_diff)))

# 上调和下调
up_LR <- diff_LR_sig %>% filter(prob_diff > 0)
down_LR <- diff_LR_sig %>% filter(prob_diff < 0)

# 保存结果
write.csv(diff_LR_sig, file = "LRpairs.csv", row.names = FALSE)

# 获取所有显著变化的pathways
all_pathways <- unique(diff_LR_sig$pathway_name)
# ===== 1. 绘制ND条件下的chord图 =====
pdf("ND_chord.pdf", width = 12, height = 12)

netVisual_chord_gene(
  ND, 
  sources.use = mac_sources,
  targets.use = epithelial_targets,
  signaling = all_pathways,
  lab.cex = 0.7,
  title.name = "Macrophage-Epithelial Communication in ND",
  legend.pos.x = 10,
  legend.pos.y = 20,
  small.gap = 1,
  big.gap = 10
)

dev.off()

# ===== 2. 绘制HFD条件下的chord图 =====
pdf("HFD_chord.pdf", width = 12, height = 12)

netVisual_chord_gene(
  HFD, 
  sources.use = mac_sources,
  targets.use = epithelial_targets,
  signaling = all_pathways,
  lab.cex = 0.7,
  title.name = "Macrophage-Epithelial Communication in HFD",
  legend.pos.x = 10,
  legend.pos.y = 20,
  small.gap = 1,
  big.gap = 10
)

dev.off()

# ===== 3. Modified Dot plot with categorical regulation =====
pdf("LR_regulation_dotplot.pdf", width = 7, height = 10)

# Prepare data with categorical regulation
plot_data <- diff_LR_sig %>%
  mutate(
    cell_pair = paste(source, "→", target),
    LR_label = paste(ligand, receptor, sep = " - "),
    
    # Categorize regulation based on prob values
    regulation = case_when(
      prob_ND == 0 & prob_HFD > 0 ~ "Only in HFD",
      prob_HFD == 0 & prob_ND > 0 ~ "Only in ND",
      prob_ND > 0 & prob_HFD > 0 & prob_HFD > prob_ND ~ "HFD > ND",
      prob_ND > 0 & prob_HFD > 0 & prob_ND > prob_HFD ~ "ND > HFD",
      TRUE ~ "No change"
    ),
    
    # Calculate fold change (avoid division by zero)
    fold_change = ifelse(prob_ND == 0, 
                         prob_HFD,  # Use HFD prob directly if ND is 0
                         ifelse(prob_HFD == 0,
                                -prob_ND,  # Use negative ND prob if HFD is 0
                                prob_HFD / prob_ND)),  # Normal ratio
    
    abs_prob_diff = abs(prob_diff)
  ) %>%
  arrange(desc(abs_prob_diff))

# Limit to top LR pairs
top_n_pairs <- 30
if(nrow(plot_data) > top_n_pairs) {
  plot_data <- plot_data %>% slice(1:top_n_pairs)
}

# Set factor levels for regulation (for plotting order)
plot_data$regulation <- factor(plot_data$regulation, 
                               levels = c("Only in HFD", "HFD > ND", 
                                          "ND > HFD", "Only in ND"))

# Nature/Science style color palette
# Red shades for HFD-dominant, Blue shades for ND-dominant
regulation_colors <- c(
  "Only in HFD" = "#D62728",      # Dark red
  "HFD > ND" = "#FF7F0E",          # Orange-red
  "ND > HFD" = "#1F77B4",          # Medium blue
  "Only in ND" = "#0D47A1"         # Dark blue
)

p1 <- ggplot(plot_data, aes(x = cell_pair, y = reorder(LR_label, abs_prob_diff))) +
  geom_point(aes(size = abs_prob_diff, color = regulation)) +
  scale_color_manual(values = regulation_colors, 
                     name = "Regulation",
                     drop = FALSE) +
  scale_size_continuous(range = c(3, 10), name = "abs(prob_diff)") +
  facet_grid(regulation ~ ., scales = "free_y", space = "free_y") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 8),
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_rect(fill = "gray95", color = "gray70"),
    panel.grid.major = element_line(color = "gray90"),
    panel.border = element_rect(color = "gray70", fill = NA),
    legend.position = "right"
  ) +
  labs(
    title = "Differential L-R pairs between ND and HFD",
    x = "Cell-cell communication",
    y = "Ligand - Receptor pair"
  )

print(p1)

dev.off()
cat("✓ Generated: LR_regulation_dotplot.pdf\n")
