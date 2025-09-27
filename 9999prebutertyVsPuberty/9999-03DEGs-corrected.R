# ============ 每细胞 MAST（非伪bulk）差异分析脚本 ============
# 工作目录按需修改
setwd("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/9999prebutertyVsPuberty")

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(stringr)
  library(rlang)
  library(MAST)
})

# ---------------- 载入对象与元信息 ----------------
allStages <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/9999Stages/allStages.rds")

meta <- allStages@meta.data

# 你的批次映射（保持不变）
meta$gse_id <- as.character(meta$batch)
meta$gse_id[meta$batch %in% c("GSM2759555_pub.rds", "GSM2759554_pub.rds")] <- "GSE_group1"
meta$gse_id[meta$batch %in% c("GSM4994963_pre.rds", "GSM4994964_pub.rds", "GSM4994965_pub.rds")] <- "GSE_group2"

condition_col <- "orig.ident"  # prepuberty / puberty
sample_col    <- "batch"       # 将作为随机效应(若可用)
batch_col     <- "gse_id"      # 固定批次协变量
celltype_col  <- "cell_type"

# 设定参考水平为 prepuberty
meta[[condition_col]] <- as.factor(meta[[condition_col]])
if (all(c("prepuberty","puberty") %in% levels(meta[[condition_col]]))) {
  meta[[condition_col]] <- relevel(meta[[condition_col]], ref = "prepuberty")
}
allStages@meta.data <- meta

# ---------------- 每细胞 MAST 主函数 ----------------
scde_mast_by_celltype <- function(seu,
                                  celltype_col = "cell_type",
                                  out_prefix   = "allStages_CellTypes",
                                  min_cells_per_type = 50,    # 每个细胞类型至少多少细胞才分析
                                  min_detect_rate    = 0.05,  # 在该细胞类型中 ≥5% 细胞检出的基因才纳入
                                  min_samples_per_cond = 2     # 每条件至少多少样本（batch）出现
                                  ) {

  # 输出目录
  out_dir <- file.path("DEG_results", paste0(out_prefix, "_MAST_cell"))
  if (!dir.exists("DEG_results")) dir.create("DEG_results")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  # 取用 log-normalized 表达（Seurat 默认 "data" 为 log1p 归一化）
  expr_log <- Seurat::GetAssayData(seu, slot = "data")
  meta_all <- seu@meta.data

  cell_types <- sort(unique(meta_all[[celltype_col]]))
  deg_list <- list()

  for (cell in cell_types) {
    message(sprintf("[MAST] Processing cell type: %s", cell))

    # 子集该细胞类型
    keep_cells <- rownames(meta_all)[meta_all[[celltype_col]] == cell]
    if (length(keep_cells) < min_cells_per_type) {
      message(sprintf("  Skip %s: too few cells (%d).", cell, length(keep_cells)))
      next
    }
    sub_meta <- meta_all[keep_cells, , drop = FALSE]

    # 条件/样本出现性检查（每条件至少 min_samples_per_cond 个样本）
    tab_cs <- table(sub_meta[[condition_col]], sub_meta[[sample_col]])
    if (any(rowSums(tab_cs) == 0)) {
      message(sprintf("  Skip %s: only one condition present.", cell)); next
    }
    samps_per_cond <- rowSums(tab_cs > 0)
    if (any(samps_per_cond < min_samples_per_cond)) {
      message(sprintf("  Skip %s: < %d biological replicates in at least one condition. (%s)",
                      cell, min_samples_per_cond,
                      paste(names(samps_per_cond), samps_per_cond, collapse = ", ")))
      next
    }

    # 基因筛选：该细胞类型中检出率 ≥ min_detect_rate
    detectable <- Matrix::rowSums(expr_log[, keep_cells, drop = FALSE] > 0) / length(keep_cells) >= min_detect_rate
    if (sum(detectable) < 50) {
      message(sprintf("  Skip %s: too few detectable genes after filter.", cell)); next
    }

    expr_sub <- as.matrix(expr_log[detectable, keep_cells, drop = FALSE])

    # 构建 SingleCellAssay（MAST 输入）
    cData <- data.frame(sub_meta, row.names = rownames(sub_meta))
    fData <- data.frame(gene = rownames(expr_sub), row.names = rownames(expr_sub))
    sca <- MAST::FromMatrix(exprsArray = expr_sub, cData = cData, fData = fData)

    # 细胞检测率（CDR）作为协变量
    colData(sca)$ngeneson <- Matrix::colSums(expr_sub > 0)
    colData(sca)$cngeneson <- scale(colData(sca)$ngeneson)

    # 条件/批次因子确保水平
    colData(sca)[[condition_col]] <- factor(colData(sca)[[condition_col]])
    if (all(c("prepuberty","puberty") %in% levels(colData(sca)[[condition_col]]))) {
      colData(sca)[[condition_col]] <- relevel(colData(sca)[[condition_col]], "prepuberty")
    }
    colData(sca)[[batch_col]]  <- factor(colData(sca)[[batch_col]])
    colData(sca)[[sample_col]] <- factor(colData(sca)[[sample_col]])

    # --- 拟合：优先用混合效应（随机截距=样本），否则回退固定效应 ---
    frm_mixed <- as.formula(sprintf("~ %s + %s + cngeneson + (1|%s)", condition_col, batch_col, sample_col))
    frm_fixed <- as.formula(sprintf("~ %s + %s + cngeneson",       condition_col, batch_col))

    fit <- try({
      zlm(formula = frm_mixed, sca = sca, method = "glmer")  # 需要 lme4 支持
    }, silent = TRUE)

    used_mixed <- !inherits(fit, "try-error")
    if (!used_mixed) {
      message(sprintf("  [%s] Mixed-effects not available; falling back to fixed-effects.", cell))
      fit <- zlm(formula = frm_fixed, sca = sca, method = "glm")
    }

    # LRT：检验 puberty 相对 prepuberty
    contrast_name <- paste0(condition_col, "puberty")  # MAST 命名通常为 'orig.identpuberty'
    s <- summary(fit, doLRT = contrast_name)
    dt <- s$datatable

    # 提取 Hurdle 组合检验的 p 值，并用连续部分(C)的系数作为 logFC 近似
    dt_p <- dt[dt$contrast == contrast_name & dt$component == "H", c("primerid", "Pr(>Chisq)")]
    dt_c <- dt[dt$contrast == contrast_name & dt$component == "C", c("primerid", "coef")]
    colnames(dt_p) <- c("gene", "pvalue")
    colnames(dt_c) <- c("gene", "avg_log2FC")  # 方向：puberty 相对 prepuberty

    tt <- dplyr::left_join(dt_c, dt_p, by = "gene")
    tt$FDR <- p.adjust(tt$pvalue, method = "fdr")
    tt$cell_type <- cell
    tt$used_mixed_model <- used_mixed

    # 上下调标签（与你原脚本一致的阈值）
    tt$regulation <- "NotSig"
    tt$regulation[tt$FDR < 0.05 & tt$avg_log2FC >  1.5] <- "Up"
    tt$regulation[tt$FDR < 0.05 & tt$avg_log2FC < -1.5] <- "Down"

    # 排序
    tt <- tt %>%
      mutate(order_group = case_when(regulation == "Up" ~ 1,
                                     regulation == "Down" ~ 2,
                                     TRUE ~ 3),
             sort_value = case_when(regulation == "Up" ~  avg_log2FC,
                                    regulation == "Down" ~ abs(avg_log2FC),
                                    TRUE ~ NA_real_)) %>%
      arrange(order_group, desc(sort_value)) %>%
      dplyr::select(-order_group, -sort_value)

    # 保存每个细胞类型的结果
    out_file <- file.path(out_dir, paste0(gsub("[/\\s]", "_", cell), "_DEGs_MAST_cell.csv"))
    write.csv(tt, file = out_file, row.names = FALSE)

    deg_list[[cell]] <- tt
  }

  # 合并所有细胞类型
  if (length(deg_list) > 0) {
    all_degs <- dplyr::bind_rows(deg_list, .id = "cell_type_name")
    write.csv(all_degs,
              file = file.path(out_dir, paste0(out_prefix, "_Combined_DEGs_MAST_cell.csv")),
              row.names = FALSE)
  }

  assign("deg_list_mast", deg_list, envir = .GlobalEnv)
  save(deg_list, cell_types, file = file.path(out_dir, paste0(out_prefix, "-envVariables_MAST_cell.RData")))
}

scde_mast_by_celltype(
  allStages,
  celltype_col = "cell_type",
  out_prefix   = "allStages_CellTypes",
  min_cells_per_type   = 50,
  min_detect_rate      = 0.05,
  min_samples_per_cond = 2
)

check_gene_celltype <- function(results_list, gene = "Esr1", celltype = "HormSens") {
  if (!celltype %in% names(results_list)) return(NULL)
  df <- results_list[[celltype]]
  df[df$gene == gene, c("gene","avg_log2FC","pvalue","FDR","regulation","used_mixed_model")]
}
esr1_hormsens <- check_gene_celltype(get("deg_list_mast"), gene = "Esr1", celltype = "HormSens")
print(esr1_hormsens)
# ===============================================================
