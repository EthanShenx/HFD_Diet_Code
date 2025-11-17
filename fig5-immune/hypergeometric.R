# =========================
# Macrophage 功能富集判定脚本
# =========================
# 依赖：tidyverse, readr, stringr, ggplot2
suppressPackageStartupMessages({
  library(tidyverse)
})

# -------- 1) 读入 gene/state 表 --------
infile <- "/Users/renyixiang/Desktop/data/23BMI/ND_HFD_MG_snRNAseq/Figure7_immune/ecotyperResult/ND_HFD/gene_info.txt"

# 自动识别分隔符
first_kb <- readChar(infile, nchars = 1000, useBytes = TRUE)
sep <- if (str_detect(first_kb, "\t") && str_count(first_kb, "\t") >= str_count(first_kb, ",")) "\t" else
  if (str_detect(first_kb, ",")) "," else "\\s+"

df <- read.table(infile, header = TRUE, sep = sep, check.names = FALSE, quote = "", comment.char = "", fill = TRUE)
names(df) <- tolower(trimws(names(df)))

# 尝试识别 gene / state 列
gene_col  <- c("gene","genes","symbol","gene_symbol","geneid","gene_id")
state_col <- c("state","cluster","macstate","mac_state","group","label")

gcol <- intersect(gene_col,  names(df))[1]
scol <- intersect(state_col, names(df))[1]
if (is.na(gcol) || is.na(scol)) stop("未找到 gene/state 两列，请检查文件列名。")

df <- df %>%
  transmute(
    gene  = toupper(trimws(.data[[gcol]])),
    state = as.character(trimws(.data[[scol]]))
  ) %>%
  filter(!is.na(gene), gene != "", !is.na(state), state != "")

states <- sort(unique(df$state))
universe <- sort(unique(df$gene))
M <- length(universe)

# -------- 2) 功能 marker 集（可按需增删） --------
marker_sets <- list(
  LAM_TREM2_SPP1 = c("TREM2","SPP1","LPL","LIPA","CD36","FABP4","FABP5","LGALS3","CTSB","CTSL","GPNMB","CST7","LILRB4A","APOE","APOC1"),
  RESIDENT_REPAIR_LYVE1_FOLR2 = c("LYVE1","FOLR2","SEPP1","MRC1","MERTK","STAB1","CD163","MAF","IGF1","VSIG4","MS4A7","MS4A6A","SELENOP","SLCO2B1","PDGFB"),
  APC_MHCII = c("H2-AA","H2-AB1","H2-EB1","CIITA","CD74","CD80","CD86","IRF1","HLA-DRA","HLA-DRB1","HLA-DPA1","HLA-DPB1"),
  IFN_ISG = c("IFIT1","IFIT2","IFIT3","ISG15","RSAD2","OASL1","OASL2","OAS1A","OAS2","LY6E","IFI44","IFI44L","IRF7","MX1","MX2","BST2"),
  LXR_PPAR_METABOLIC = c("PPARG","PPARA","RXRA","ABCA1","ABCG1","APOE","APOA1","PLIN2","NR1H3","NR1H2","CPT1A","CD36","LPL"),
  PROLIFERATION_CYCLE = c("MKI67","TOP2A","PCNA","TK1","AURKA","AURKB","CDK1","CCNB1","CCNB2","HMGB2","BIRC5","UBE2C","MYC","E2F1","MCM2","MCM4","MCM6")
)

# 将 marker 集限制到 universe，并去重
marker_sets <- lapply(marker_sets, function(v) unique(intersect(toupper(v), universe)))

# -------- 3) 逐 state 计算超几何富集 --------
phyper_sf <- function(k, n, N, M) { # P[X >= k]
  # X ~ Hypergeom(M, n, N)
  if (k <= 0) return(1)
  # 使用下尾概率转换：sf = phyper(k-1, n, M-n, N, lower.tail = FALSE)
  stats::phyper(q = k - 1, m = n, n = M - n, k = N, lower.tail = FALSE)
}

res <- purrr::imap_dfr(marker_sets, function(markers, cat){
  if (length(markers) == 0) {
    tibble(category = cat, state = states, state_size = NA_integer_, marker_in_universe = 0,
           overlap = 0, odds_ratio = NA_real_, p_enrich = NA_real_)
  } else {
    n <- length(markers)
    purrr::map_dfr(states, function(s){
      genes_s <- unique(df$gene[df$state == s])
      N <- length(genes_s)
      k <- length(intersect(genes_s, markers))
      # 2x2 表 + 0.5 连续性校正
      a <- k
      b <- n - k
      c <- N - k
      d <- M - n - c
      odds <- ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
      
      pval <- phyper_sf(k, n, N, M)
      tibble(category = cat, state = s, state_size = N, marker_in_universe = n,
             overlap = k, odds_ratio = odds, p_enrich = pval)
    })
  }
}) %>%
  mutate(neglog10_p = if_else(is.finite(-log10(p_enrich)), -log10(p_enrich), 0))

# -------- 4) 每个 state 的最佳功能标签 --------
best_calls <- res %>%
  group_by(state) %>%
  arrange(desc(neglog10_p), desc(odds_ratio), desc(overlap), .by_group = TRUE) %>%
  slice(1L) %>%
  ungroup() %>%
  transmute(
    state,
    best_category = category,
    overlap = overlap,
    neglog10_p = round(neglog10_p, 3),
    odds_ratio = round(odds_ratio, 3)
  )

print(best_calls)

# -------- 5) 输出文件 --------
outdir <- "/Users/renyixiang/Desktop/data/23BMI/ND_HFD_MG_snRNAseq/Figure7_immune/ecotyperResult/ND_HFD/mac_function_outputs"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

#write.csv(res, file.path(outdir, "mac_state_function_enrichment_summary.csv"), row.names = FALSE)
# 重叠计数矩阵
count_mat <- res %>%
  select(category, state, overlap) %>%
  pivot_wider(names_from = state, values_from = overlap)
#write.csv(count_mat, file.path(outdir, "mac_state_marker_overlap_counts.csv"), row.names = FALSE)
# 最佳判定
#write.csv(best_calls, file.path(outdir, "mac_state_best_function_call.csv"), row.names = FALSE)

# -------- 6) 画热图（-log10 p） --------
heat <- res %>%
  select(category, state, neglog10_p) %>%
  pivot_wider(names_from = state, values_from = neglog10_p) %>%
  arrange(category)

heat_long <- res %>%
  select(category, state, neglog10_p)

p <- ggplot(heat_long, aes(x = state, y = category, fill = neglog10_p)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.1f", neglog10_p)), size = 3) +
  RotatedAxis() + 
  scale_fill_distiller(palette = "Spectral", direction = -1, name = "-log10(p)")+
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(title = "Functional marker enrichment by state",
       x = "State", y = "Marker category")
#ggsave(file.path(outdir, "mac_state_function_enrichment_heatmap.pdf"),
       #p, width = max(6, 0.9*length(states)), height = max(3.5, 0.6*length(marker_sets)), dpi = 220)

# -------- 7) 导出 GMT（可直接用于 AUCell/UCell 等） --------
gmt_path <- file.path(outdir, "mac_function_markers.gmt")
con <- file(gmt_path, open = "wt", encoding = "UTF-8")
for (nm in names(marker_sets)) {
  line <- paste(c(nm, "generated_by_script", marker_sets[[nm]]), collapse = "\t")
  writeLines(line, con)
}
close(con)

message("完成。输出目录： ", outdir)
save.image(file = "hyper_geomertic.RData")
