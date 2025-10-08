library(Seurat)
library(ggplot2)
library(ggpubr)

gene   <- "Gas6"
df.exp <- FetchData(adipo, vars = c(gene, "orig.ident")) %>% 
  as.data.frame() %>% 
  mutate(orig.ident = factor(orig.ident, levels = c("ND", "HFD")))

p_vio <- ggplot(df.exp,
                aes(x = orig.ident,
                    y = .data[[gene]],
                    fill = orig.ident)) +
  geom_violin(trim = FALSE, alpha = 0.7, color = "black", linewidth = .2) +
  geom_boxplot(width = 0.15, fill = "white", color = "black",
               linewidth = .2, outlier.shape = NA) +
  
  stat_compare_means(method = "wilcox.test",
                     comparisons = list(c("ND", "HFD")),
                     label = "p.format",
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                        symbols = c("***", "**", "*", "ns")),
                     vjust = 1.5) +
  
  labs(x = NULL,
       y = " expression level",
       title = gene) +
  
  scale_fill_manual(values = c("ND" = "#00BFC4", "HFD" = "#F8766D")) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank())

print(p_vio)
