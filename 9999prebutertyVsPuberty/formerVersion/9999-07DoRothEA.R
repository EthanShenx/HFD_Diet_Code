library(Seurat)
library(dorothea)
library(decoupleR)
library(dplyr)
library(tidyr)
library(purrr)
library(pheatmap)
library(igraph)

# cds_obj <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/9999prebutertyVsPuberty/formerVersion/trajectoryOutput/cds_obj.rds")

All_stage <- readRDS(file = "/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/9999prebutertyVsPuberty/formerVersion/trajectoryOutput/All_stages_trajectory.rds")

load("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/9999prebutertyVsPuberty/formerVersion/proteinNetwork/up_data.rda")
load("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/9999prebutertyVsPuberty/formerVersion/proteinNetwork/df_up_common.rda")
load("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/9999prebutertyVsPuberty/formerVersion/proteinNetwork/orthologs.rda")

DefaultAssay(All_stage) <- "RNA"
expr <- as.matrix(GetAssayData(All_stage, slot = "data"))

## ===== construct DoRothEA network (mouse) =====
data(dorothea_mm)
data(dorothea_hs)

regulons_abc <- dorothea_mm %>%
  dplyr::filter(confidence %in% c("A","B","C")) %>%
  dplyr::transmute(tf = tf, target = target, mor = mor, confidence = confidence)

reg_hs <- dorothea_hs %>%
  dplyr::filter(confidence %in% c("A","B","C","D","E")) %>%
  dplyr::select(tf, target, mor, confidence)

up_hs <- orthologs %>% dplyr::filter(mouse_symbol %in% up_mouse) %>% dplyr::pull(human_symbol) %>% unique()
bg_hs <- orthologs %>% dplyr::filter(mouse_symbol %in% rownames(expr)) %>% dplyr::pull(human_symbol) %>% unique()

reg_hs_pos <- reg_hs %>% dplyr::filter(mor > 0, target %in% bg_hs)

common_genes <- intersect(rownames(expr), unique(regulons_abc$target))
expr <- expr[common_genes, , drop = FALSE]
regulons_abc <- regulons_abc %>% dplyr::filter(target %in% common_genes)

# --- collapse duplicate TFs with identical target(+mor) signatures, prefer high-confidence
sig_tbl <- regulons_abc %>%
  arrange(tf, target) %>%
  group_by(tf) %>%
  summarise(sig = paste(target, mor, collapse="|"), .groups = "drop")

dup_groups <- sig_tbl %>% group_by(sig) %>% filter(n() > 1)

if (nrow(dup_groups) > 0) {
  conf_w <- c(A=3,B=2,C=1)
  keep_dup <- dup_groups %>%
    left_join(regulons_abc, by="tf") %>%
    group_by(sig, tf) %>%
    summarise(A_targets = sum(confidence=="A"),
              total_targets = n(),
              score = sum(conf_w[confidence]), .groups="drop") %>%
    group_by(sig) %>%
    slice_max(order_by = A_targets, n=1, with_ties = FALSE) %>%
    pull(tf)
  drop_dup <- setdiff(dup_groups$tf, keep_dup)
  regulons_abc <- regulons_abc %>% filter(!tf %in% drop_dup)
}

# --- drop highly collinear TFs (identical/near-identical signatures)
cc <- decoupleR::check_corr(
  regulons_abc, 
  .source="tf", 
  .target="target", 
  .mor="mor"
)
bad <- cc %>% dplyr::filter(abs(correlation) >= 0.9)

if (nrow(bad) > 0) {
  g <- graph_from_data_frame(bad %>% transmute(from=source, to=source.2), directed=FALSE)
  comps <- components(g)$membership
  groups <- split(names(comps), comps)
  conf_w <- c(A=3,B=2,C=1)
  keep <- purrr::map_chr(groups, function(gs){
    regulons_abc %>% dplyr::filter(tf %in% gs) %>%
      group_by(tf) %>%
      summarise(A_targets = sum(confidence=="A"),
                total_targets = n(),
                score = sum(conf_w[confidence]), .groups="drop") %>%
      arrange(desc(A_targets), desc(score), desc(total_targets)) %>%
      slice(1) %>% pull(tf)
  })
  drop <- setdiff(unique(c(bad$source, bad$source.2)), keep)
  message("Dropping ", length(drop), " highly collinear TFs: ", paste(head(drop, 10), collapse=", "),
          ifelse(length(drop)>10, " ...", ""))
  regulons_abc <- regulons_abc %>% dplyr::filter(!tf %in% drop)
}

up_mouse <- df_up_common %>%
  dplyr::transmute(gene = as.character(gene)) %>%
  dplyr::pull(gene) %>% unique()

bg_mouse <- rownames(expr)

reg_pos <- regulons_abc %>%
  dplyr::filter(mor > 0, target %in% bg_mouse)

minSize <- 10L
tf_sets <- reg_pos %>%
  dplyr::group_by(tf) %>%
  dplyr::summarise(targets = list(unique(target)),
                   n_target = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n_target >= minSize)

up_in_bg <- intersect(up_mouse, bg_mouse)
N <- length(bg_mouse); n <- length(up_in_bg)

dorothea_up_common_TF_ora <- purrr::pmap_dfr(
  list(tf_sets$tf, tf_sets$targets, tf_sets$n_target),
  function(tf, targets, M){
    k <- sum(up_in_bg %in% targets)
    a <- k; b <- M - k
    c <- n - k; d <- N - M - (n - k)
    p <- stats::fisher.test(matrix(c(a,b,c,d), nrow = 2),
                            alternative = "greater")$p.value
    tibble(
      tf = tf,
      k = k, M = M, n = n, N = N,
      overlap_ratio = ifelse(M > 0, k/M, NA_real_),
      overlap_genes = paste(intersect(up_in_bg, targets), collapse = ";"),
      pvalue = p
    )
  }
) %>%
  dplyr::mutate(padj = p.adjust(pvalue, method = "BH")) %>%
  dplyr::arrange(padj, desc(k), desc(overlap_ratio))

# 输出结果文件（可按需修改路径）
# readr::write_csv(dorothea_up_common_TF_ora,
#                  "/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/9999prebutertyVsPuberty/formerVersion/proteinNetwork/dorothea_up_common_TF_ora.csv")
# save(dorothea_up_common_TF_ora,
#      file = "/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/9999prebutertyVsPuberty/formerVersion/proteinNetwork/dorothea_up_common_TF_ora.rda")

## ===== decoupleR::run_mlm to estimate TF activity (per cell) =====

mlm_res <- run_mlm(
  mat       = expr,
  network   = regulons_abc,
  .source   = "tf",
  .target   = "target",
  .mor      = "mor",
  minsize   = 5
)

mlm_mat <- mlm_res %>%
  dplyr::select(source, condition, score) %>%
  dplyr::rename(tf = source, cell = condition, score = score) %>%
  tidyr::pivot_wider(names_from = cell, values_from = score) %>%
  tibble::column_to_rownames("tf") %>%
  as.matrix()

All_stage[["dorothea_mlm"]] <- CreateAssayObject(data = mlm_mat)
