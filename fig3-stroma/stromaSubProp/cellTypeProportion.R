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

################################
########## Statistics ##########
################################

library(dplyr)
library(tidyr)
library(purrr)
library(stringr)

counts_by <- tmp2 %>%
  dplyr::select(orig.ident, subcluster, count) %>%
  tidyr::complete(orig.ident = factor(c("ND","HFD"), levels = c("ND","HFD")),
                  subcluster,
                  fill = list(count = 0))

totals <- counts_by %>%
  dplyr::group_by(orig.ident) %>%
  dplyr::summarise(total = sum(count), .groups = "drop") %>%
  tibble::deframe()  # named numeric: c(ND = ..., HFD = ...)

wide_ct <- counts_by %>%
  tidyr::pivot_wider(names_from = orig.ident, values_from = count, values_fill = 0) %>%
  dplyr::arrange(subcluster)

eps <- 1e-6
test_list <- purrr::pmap(
  list(wide_ct$subcluster, wide_ct$ND, wide_ct$HFD),
  function(ct, nd_ct, hfd_ct) {

    nd_other  <- totals["ND"]  - nd_ct
    hfd_other <- totals["HFD"] - hfd_ct

    mat <- matrix(c(hfd_ct, hfd_other,
                    nd_ct,  nd_other),
                  nrow = 2, byrow = TRUE,
                  dimnames = list(Group = c("HFD","ND"),
                                  Category = c(ct, paste0("non-", ct))))

    ft <- fisher.test(mat, alternative = "two.sided")

    # Haldane-Anscombe 校正后的 OR，避免 0 导致 Inf
    or_ha <- ((hfd_ct + 0.5)/(hfd_other + 0.5)) /
             ((nd_ct  + 0.5)/(nd_other  + 0.5))

    nd_rate  <- nd_ct  / totals["ND"]
    hfd_rate <- hfd_ct / totals["HFD"]

    tibble::tibble(
      subcluster = ct,
      ND_count   = nd_ct,
      HFD_count  = hfd_ct,
      ND_total   = unname(totals["ND"]),
      HFD_total  = unname(totals["HFD"]),
      ND_rate    = nd_rate,
      HFD_rate   = hfd_rate,
      rate_diff  = hfd_rate - nd_rate,
      log2FC_prop = log2((hfd_rate + eps)/(nd_rate + eps)),
      odds_ratio = unname(or_ha),
      p_value    = ft$p.value,
      conf_low   = if (!is.null(ft$conf.int)) ft$conf.int[1] else NA_real_,
      conf_high  = if (!is.null(ft$conf.int)) ft$conf.int[2] else NA_real_
    )
  }
)

prop_test_df <- dplyr::bind_rows(test_list) %>%
  dplyr::mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  dplyr::arrange(p_adj)

write.csv(prop_test_df, "cell_prop_fisher_res.csv", row.names = FALSE)

##########

setwd("/Users/coellearth/Desktop/HFD_Paper/Fig1_snRNA-seqProfiling/cellTypeProportion")

library(Seurat)
library(dplyr)
library(ggplot2)
library(scales)

All <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_All_sub.rds")
meta <- All@meta.data
meta$cell_type <- NULL

ct_col <- "subcluster"

tmp2 <- meta %>%
  dplyr::count(orig.ident, subcluster, name = "count") %>%
  dplyr::group_by(orig.ident) %>%
  dplyr::mutate(total = sum(count),
                rate  = count / total) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    orig.ident = factor(orig.ident,
                        levels = c("ND", "HFD"),
                        labels = c("ND", "HFD"))
  )

ggplot(
  data = tmp2,
  aes(x = subcluster, 
      y = rate, 
      fill = orig.ident)
) +
  geom_col(position = position_dodge(width = 0),  # bar 紧贴并置
           width = 0.9) +
  scale_fill_discrete() +  # 使用 ggplot2 默认色板
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = "Cell type",
    y = "Proportion",
    fill = "Diet",
    title = "Cell composition"
  ) +
  coord_flip() +  # 横向摆放
  theme_classic(base_size = 14) +
  theme(
    text = element_text(family = "Arial"),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5),
    legend.position = "right",
    legend.title = element_blank(),
    legend.key.size = unit(0.8, "lines"),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3)
  )

## ---------------- 改动要点 ----------------
## 1) 不再使用 ggalluvial / RColorBrewer，也不再自定义调色
## 2) 用 geom_col() 画两个并置的 stacked bar
## 3) 每个堆叠段加黑色描边 color="black"
## 4) 使用 ggplot2 默认离散色板（不设 scale_fill_manual）
## 5) 两个柱之间留很小的间隔：width=0.95
## 6) 仍然显示比例（rate）与百分号刻度
## -----------------------------------------

# (可选) 保证 subcluster 的顺序稳定
tmp2 <- tmp2 %>%
  dplyr::mutate(subcluster = factor(subcluster, levels = sort(unique(subcluster))))

p <- ggplot(
  data = tmp2,
  aes(x = orig.ident, y = rate, fill = subcluster)
) +
  # 两个并置堆叠柱；黑线勾勒每个 segment
  geom_col(width = 0.95, color = "black", linewidth = 0.2) +

  # 使用 ggplot2 默认离散色板（不写 scale_fill_manual）
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
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
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5),
    legend.position = "right",
    legend.title = element_blank(),
    legend.key.size = unit(0.8, "lines"),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3)
  )

print(p)

p + coord_flip() + geom_flow(alpha = 0.55)
