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

# Prepare data for Sankey plot
sankey_data <- metadata %>%
  group_by(orig.ident, celltype_state) %>%
  summarise(count = n(), .groups = 'drop')

# NPG (Nature) color palette for 4 cell types
cell_colors <- c("MacS01" = "#E64B35",  # Red
                 "MacS02" = "#4DBBD5",  # Blue
                 "MacS03" = "#00A087",  # Teal/Green
                 "MacS04" = "#F39B7F")  # Orange

# Colors for orig.ident (ND vs HFD)
sample_colors <- c("#3C5488", "#DC0000")  # Navy blue for ND, Red for HFD

# Create Sankey plot
p <- ggplot(sankey_data,
            aes(axis1 = orig.ident, 
                axis2 = celltype_state, 
                y = count)) +
  geom_alluvium(aes(fill = celltype_state), 
                width = 1/12, 
                alpha = 0.8, 
                curve_type = "sigmoid") +
  geom_stratum(width = 1/12, 
               fill = "white", 
               color = "black", 
               linewidth = 0.5) +
  geom_text(stat = "stratum", 
            aes(label = after_stat(stratum)), 
            size = 4,
            fontface = "bold") +
  scale_x_discrete(limits = c("Sample", "Cell Type"),
                   expand = c(0.15, 0.05)) +
  scale_fill_manual(values = cell_colors) +
  labs(title = "Macrophage Subtype Composition Change",
       subtitle = paste0("Chi-square test: p ", 
                         ifelse(p_value < 0.001, "< 0.001", 
                                ifelse(p_value < 0.01, paste("=", round(p_value, 3)),
                                       paste("=", round(p_value, 4))))),
       fill = "Cell Type",
       y = "Cell Count") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 13, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 11),
        plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 11),
        legend.position = "right",
        legend.title = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 10))

# Add significance annotation for the entire diet vs. subtype difference
max_y <- sum(sankey_data$count) / length(unique(sankey_data$orig.ident))
p <- p + 
  annotate("text", 
           x = 1.5,  # Position between the two axes
           y = max_y * 1.15,  # Above the flows
           label = sig_label, 
           size = 12,  # Larger size for visibility
           color = ifelse(sig_label == "ns", "grey30", "red"),
           fontface = "bold")

# Print the plot
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
