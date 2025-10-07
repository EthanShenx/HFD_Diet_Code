library(ggplot2)
library(dplyr)

## 1. 数据整理（同上） -------------------------------------------------------
meta <- adipo@meta.data %>% 
  as.data.frame() %>% 
  select(orig.ident, subcluster) %>% 
  mutate(orig.ident = factor(orig.ident, levels = c("ND", "HFD")))

total_df <- dplyr::count(meta, orig.ident, name = "total")

count_tbl <- meta %>% 
  dplyr::count(orig.ident, subcluster) %>% 
  left_join(total_df, by = "orig.ident") %>% 
  mutate(prop   = n / total,
         pct_lbl = sprintf("%.1f%%", prop * 100)) %>% 
  arrange(orig.ident, prop) %>% 
  group_by(orig.ident) %>% 
  mutate(y_mid = cumsum(prop) - prop / 2)

## 2. 画图 -------------------------------------------------------------------
p <- ggplot(count_tbl,
            aes(x = orig.ident,
                y = prop,
                fill = subcluster)) +
  # 2.1 柱子
  geom_col(width = 0.7, color = "black", linewidth = 0.1) +   # linewidth 替代 size
  # 2.3 柱顶总细胞数 —— 切断继承，避免去找 subcluster
  geom_text(data = total_df,
            aes(x = orig.ident, y = 1.05, label = sprintf("n = %d", total)),
            inherit.aes = FALSE, vjust = 0, size = 3.5) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.25),      # 固定 0、25%、50%、75%、100%
    labels = scales::percent,
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(x = NULL,
       y = "Proportion within diet group",
       title = "Adipocyte sub-cluster composition") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank())

print(p)
