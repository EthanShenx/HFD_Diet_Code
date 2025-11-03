#——————————————————————————————————————————————————————
# Complete R script: Top25 Ligand–Receptor dot‐plot
#  • X‐axis sorted by receptor then ligand
#  • Ligand/receptor ticks with connecting black lines
# Cell type order: "Adipo","LumProg","Endo","HormSens","Stroma","Immune","Basal"
# Color palette: RColorBrewer “Paired”
#——————————————————————————————————————————————————————

# 0. Load libraries
library(tidyverse)
library(data.table)
library(viridis)
library(RColorBrewer)

# 1. Read in CellPhoneDB outputs (HFD)
means <- fread("D:/data/23BMI/ND_HFD_MG_snRNAseq/cellphonedb/result/statistical_analysis_HFD/statistical_analysis_means_10_09_2025_111445.txt")
pvals <- fread("D:/data/23BMI/ND_HFD_MG_snRNAseq/cellphonedb/result/statistical_analysis_HFD/statistical_analysis_pvalues_10_09_2025_111445.txt")
sig   <- fread("D:/data/23BMI/ND_HFD_MG_snRNAseq/cellphonedb/result/statistical_analysis_HFD/statistical_analysis_significant_means_10_09_2025_111445.txt")

# 2. 筛出 Ligand→Receptor
means_lr <- means %>% filter(directionality=="Ligand-Receptor", secreted, receptor_b)
pvals_lr <- pvals %>% filter(directionality=="Ligand-Receptor", secreted, receptor_b)
sig_lr   <- sig   %>% filter(directionality=="Ligand-Receptor", secreted, receptor_b)

# 3. pivot significant_means → long
sig_long <- sig_lr %>%
  pivot_longer(matches("\\|"), names_to="cluster_pair", values_to="mean_sig") %>%
  filter(mean_sig>0)

# 4. pivot pvalues → long
pvals_long <- pvals_lr %>%
  pivot_longer(matches("\\|"), names_to="cluster_pair", values_to="pvalue")

# 5. join mean_sig & pvalue
df <- sig_long %>%
  inner_join(pvals_long,
             by = c("interacting_pair","gene_a","gene_b","cluster_pair"))

# 6. 选 Top25 interacting_pair
top25_pairs <- df %>%
  group_by(interacting_pair) %>%
  summarise(max_mean = max(mean_sig, na.rm=TRUE)) %>%
  arrange(desc(max_mean)) %>%
  slice_head(n=25) %>%
  pull(interacting_pair)

df25 <- df %>% filter(interacting_pair %in% top25_pairs)

# 7. 拆分 sender/receiver（保留原 cluster_pair）
df25 <- df25 %>%
  separate(cluster_pair, into=c("sender","receiver"), sep="\\|", remove=FALSE)

# 8. 定义顺序 & 调色板
cell_order  <- c("Adipo","LumProg","Endo","HormSens","Stroma","Immune","Basal")
paired_cols <- brewer.pal(7, "Paired")
names(paired_cols) <- cell_order

# 9. 因子化并按 receiver→sender 排序
df25 <- df25 %>%
  mutate(sender=factor(sender, levels=cell_order),
         receiver=factor(receiver, levels=cell_order)) %>%
  arrange(receiver, sender) %>%
  mutate(cluster_pair=factor(cluster_pair, levels=unique(cluster_pair)),
         xpos=as.numeric(cluster_pair))

# 10. 计算每个 receptor 组的 xmin/xmax/xmid，用于连线和 receptor tick
receptor_groups <- df25 %>%
  group_by(receiver) %>%
  summarise(xmin=min(xpos)-0.4, xmax=max(xpos)+0.4) %>%
  mutate(xmid=(xmin+xmax)/2)

# 11. 开始绘主图
p <- ggplot(df25, aes(
  x     = cluster_pair,
  y     = interacting_pair,
  size  = mean_sig,
  color = pvalue
)) +
  geom_point() +
  scale_size_continuous(range=c(1,6), name="Avg. expr.") +
  scale_color_viridis_c(direction=-1, option="D",
                        name="p-value", limits=c(0,0.05),
                        oob=scales::squish) +
  theme_minimal(base_size=12) +
  theme(
    axis.text.x        = element_blank(),
    axis.ticks.x       = element_blank(),
    axis.title.x       = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.margin        = margin(t=5,r=5,b=100,l=5)
  ) +
  labs(y="Top 25 Ligand–Receptor pairs")

# 12. 黑线连接同一 receiver 下的所有 cluster_pairs
p <- p +
  geom_segment(data=receptor_groups,
               aes(x=xmin, xend=xmax, y=-1.05, yend=-1.05),
               color="black", inherit.aes=FALSE)

# 13. Ligand ticks + 数字
sender_df <- df25 %>% distinct(cluster_pair, sender)
p <- p +
  geom_point(data=sender_df,
             aes(x=cluster_pair, y=-0.5, fill=sender),
             shape=21, size=3, stroke=0.3, inherit.aes=FALSE) +
  geom_text(data=sender_df %>% mutate(code=as.numeric(sender)),
            aes(x=cluster_pair, y=-0.5, label=code),
            color="black", size=2, inherit.aes=FALSE)

# 14. Receptor ticks + 数字
p <- p +
  geom_point(data=receptor_groups,
             aes(x=xmid, y=-1.3, fill=receiver),
             shape=21, size=3, stroke=0.3, inherit.aes=FALSE) +
  geom_text(data=receptor_groups %>% mutate(code=as.numeric(receiver)),
            aes(x=xmid, y=-1.3, label=code),
            color="black", size=2, inherit.aes=FALSE)

# 15. 最终 fill scale & 显示面板外元素
p +
  scale_fill_manual(name="Cell type",
                    values=paired_cols,
                    breaks=cell_order) +
  coord_cartesian(clip="off")
