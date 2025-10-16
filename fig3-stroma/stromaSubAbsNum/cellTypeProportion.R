setwd("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/fig3-stroma/stromaSubProp")

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggalluvial)
library(RColorBrewer)
library(scales)

All <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_All_Stroma_sub.rds")

meta <- All@meta.data

meta$cell_type <- NULL

ct_col <- "subcluster"

tmp2 <- meta %>%
  dplyr::count(orig.ident, subcluster, name = "count") %>%
  dplyr::group_by(orig.ident) %>%
  dplyr::mutate(total = sum(count),
         rate  = count / total) %>%
  dplyr::ungroup()

tmp2 <- tmp2 %>%
  dplyr::mutate(
    orig.ident = factor(orig.ident,
                        levels = c("ND", "HFD"),
                        labels = c("ND", "HFD"))
  )

n_ct <- length(unique(tmp2$subcluster))
base_n <- min(max(3, n_ct), 12)
pal <- colorRampPalette(brewer.pal(base_n, "Paired"))(n_ct)
colors <- setNames(pal, sort(unique(tmp2$subcluster)))

ggplot(
  data = tmp2,
  aes(x = orig.ident,
      stratum  = subcluster,
      alluvium = subcluster,
      y = rate,
      fill = subcluster,
      label = subcluster)
) +
  geom_flow(alpha = 0.55) +
  geom_stratum(width = 0.4,
               color = "white",
               linewidth = 0) +
  scale_fill_discrete() +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = "Diet type",
    y = "Cell type proportion",
    title = "Cell composition"
  ) +
  theme_classic(base_size = 14) +
  theme(
    text = element_text(family = "Arial"),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5),
    legend.position = "right",
    legend.title = element_blank(),
    legend.key.size = unit(0.8, "lines"),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3)
  ) +
  geom_col(width = 0.4, color = "black", linewidth = 0.2)

ggplot(tmp2, aes(x = orig.ident, y = count, fill = subcluster)) +
  geom_col(width = 0.6, color = "black", linewidth = 0.2) +
  scale_fill_manual(values = colors) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Diet type", y = "Cell count", title = "Cell composition (absolute)") +
  theme_classic(base_size = 14) +
  theme(
    text = element_text(family = "Arial"),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5),
    legend.position = "right",
    legend.title = element_blank()
  )