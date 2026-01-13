library(Seurat)
library(ggplot2)
library(dplyr)
library(ggalluvial)

# Read the data
obj <- readRDS("/Users/renyixiang/Desktop/data/23BMI/ND_HFD_MG_snRNAseq/Figure7_immune/All_for_CCC.rds")
Idents(obj) <- "celltype_state"

# Filter for the specific macrophage subtypes
mac_subtypes <- c("MacS01", "MacS02", "MacS03", "MacS04")
obj_sub <- subset(obj, idents = mac_subtypes)

# Extract metadata
metadata <- obj_sub@meta.data %>%
  select(orig.ident, celltype_state) %>%
  mutate(celltype_state = factor(celltype_state, levels = mac_subtypes))

# Create frequency table of macrophage subtypes (by diet condition)
freq_table <- table(metadata$orig.ident, metadata$celltype_state)
print("Frequency table:")
print(freq_table)

# Perform Chi-square test on the contingency table between diet conditions and macrophage subtypes
chi_test <- chisq.test(freq_table)
print("\nChi-square test results:")
print(chi_test)

# Assign significance level based on the p-value of the Chi-square test for the entire contingency table
p_value <- chi_test$p.value
sig_label <- ifelse(p_value < 0.001, "***", 
                    ifelse(p_value < 0.01, "**", 
                           ifelse(p_value < 0.05, "*", "ns")))

library(scales)  # 用来把 x 轴显示成百分比

## 1. 统计 ND / HFD 中各 MacS0* 的比例 --------------------------
comp_data <- metadata %>%
  group_by(orig.ident, celltype_state) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  mutate(freq = count / sum(count))  # 每个组内标准化到 0–1

# 看一眼数据
print(comp_data)

## 2. 画 100% 堆积条形图（两条横杠：ND 和 HFD） -----------------
# 颜色仍然用你之前的 NPG 调色板
cell_colors <- c("MacS01" = "#E64B35",  # Red
                 "MacS02" = "#4DBBD5",  # Blue
                 "MacS03" = "#00A087",  # Teal/Green
                 "MacS04" = "#F39B7F")  # Orange

p <- ggplot(comp_data,
            aes(y = orig.ident,   # ND / HFD 放在纵轴
                x = freq,         # 横轴 0–100%
                fill = celltype_state)) +
  geom_col(width = 0.6,
           color = "black",
           linewidth = 0.3) +
  scale_x_continuous(labels = percent_format(accuracy = 1),
                     expand = c(0, 0)) +
  scale_fill_manual(values = cell_colors, name = "Mac subtype") +
  labs(
    x = NULL,
    y = NULL,
    title = "Macrophage subtype composition",
    subtitle = paste0(
      "Chi-square test: p ",
      ifelse(p_value < 0.001, "< 0.001",
             ifelse(p_value < 0.01, paste("= ", round(p_value, 3)),
                    paste("= ", round(p_value, 4))))
    )
  ) +
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    legend.position = "right",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10)
  )

## 3. 加整体显著性标记（放在 ND / HFD 中间上方） ----------------
p <- p +
  annotate(
    "text",
    x = 0.5,       # 横轴中间（50% 的位置）
    y = 1.5,       # 介于第一行和第二行之间
    label = sig_label,
    size = 6,
    fontface = "bold",
    color = ifelse(sig_label == "ns", "grey30", "red")
  )

print(p)

# Save the plot
#ggsave("Sankey_MacSubtypes_ChiSquare_Modified.pdf", 
#       plot = p, 
#       width = 12, 
#       height = 8, 
#       dpi = 300)

#ggsave("Sankey_MacSubtypes_ChiSquare_Modified.png", 
#       plot = p, 
#       width = 12, 
#       height = 8, 
#       dpi = 300)
save.image(file = "Sankey.RData")
