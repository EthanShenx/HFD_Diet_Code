All <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_All_sub.rds")

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)

genes <- c("Areg","Sp1","Egr1","Esr1","Cited1")

seu <- subset(
  All,
  subset = subcluster == "HormSens" 
  & 
  orig.ident %in% c("ND","HFD")
)


DefaultAssay(seu) <- "RNA"
slot_to_use <- "data"

genes <- intersect(genes, rownames(seu))
stopifnot(length(genes) > 0)

df <- FetchData(seu, vars = c("orig.ident", genes), slot = slot_to_use) |>
  rownames_to_column("cell") |>
  pivot_longer(cols = all_of(genes), names_to = "gene", values_to = "expr") |>
  mutate(
    orig.ident = factor(orig.ident, levels = c("ND","HFD")),
    gene       = factor(gene, levels = rev(genes))
  )

xlim <- range(df$expr, finite = TRUE)
pad  <- diff(xlim) * 0.05
xlim <- c(xlim[1] - pad, xlim[2] + pad)

xlim <- c(-2, 6)

p <- ggplot(df, aes(x = expr, y = gene, fill = orig.ident)) +
  
  geom_violin(scale = "width", 
              trim = F,
              color = NA) +
  
  geom_boxplot(width = 0.25, 
               fill = "white",
               outlier.size = 0.5, 
               outlier.color = "black",
               outlier.shape = 19, 
               linewidth = 0.3) +
  
  facet_grid(. ~ orig.ident, scales = "fixed") +
  
  scale_fill_manual(values = c(ND = "#74c5be", HFD = "#e95503"), guide = "none") +
  
  scale_x_continuous(limits = xlim, expand = expansion(mult = 0)) +
  
  labs(x = "Expression level", y = NULL) +
  
  theme_classic(base_size = 12) +
  
  theme(
    strip.text = element_text(),
    axis.text.y = element_text(face = "italic", color = "black"),
    panel.spacing.x = unit(1.2, "lines")
  )

print(p)

