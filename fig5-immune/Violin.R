# 需要的包
if (!requireNamespace("gghalves", quietly = TRUE)) install.packages("gghalves")
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(gghalves)
library(stringr)

# ---- 固定分组顺序 ----
SENDER_LEVELS   <- c("MacS01", "MacS03")
RECEIVER_LEVELS <- c("LumProg", "Basal", "HormSens")

# ---- 读对象（保持你的路径）----
setwd("/Users/renyixiang/Desktop/data/23BMI/ND_HFD_MG_snRNAseq/Figure7_immune")
obj <- readRDS("/Users/renyixiang/Desktop/data/23BMI/ND_HFD_MG_snRNAseq/Figure7_immune/All_for_CCC.rds")
DefaultAssay(obj) <- "RNA"
Idents(obj) <- "celltype_state"

# ---- 明确 sender / receiver 并设定因子顺序 ----
senders   <- subset(obj, idents = SENDER_LEVELS)
receivers <- subset(obj, idents = RECEIVER_LEVELS)
senders$celltype_state   <- factor(senders$celltype_state,   levels = SENDER_LEVELS)
receivers$celltype_state <- factor(receivers$celltype_state, levels = RECEIVER_LEVELS)

# ---- 生成 Condition（ND/HFD） ----
add_condition <- function(object) {
  if (!"Condition" %in% colnames(object@meta.data)) {
    object$Condition <- ifelse(grepl("HFD", object$orig.ident, ignore.case = TRUE), "HFD", "ND")
  }
  object$Condition <- factor(object$Condition, levels = c("ND", "HFD"))
  object
}
senders   <- add_condition(senders)
receivers <- add_condition(receivers)

# ---- 基因存在性检查 ----
gene_present <- function(seurat_obj, candidates) {
  m <- intersect(candidates, rownames(seurat_obj))
  if (length(m) == 0) return(NULL)
  m[1]
}

# ---- 取作图数据（按照给定的分组顺序）----
prep_violin_df <- function(object, feature, group_levels) {
  df <- FetchData(object, vars = c(feature, "celltype_state", "Condition"))
  colnames(df) <- c("expr", "group", "Condition")
  df$group <- factor(df$group, levels = group_levels)
  df
}

# ---- 半小提琴（ND 左半，HFD 右半）----
half_violin_plot <- function(df, title, ylab) {
  ggplot(df, aes(x = group, y = expr)) +
    gghalves::geom_half_violin(
      data = dplyr::filter(df, Condition == "ND"),
      side = "l", aes(fill = Condition), width = 0.9, trim = TRUE
    ) +
    gghalves::geom_half_violin(
      data = dplyr::filter(df, Condition == "HFD"),
      side = "r", aes(fill = Condition), width = 0.9, trim = TRUE
    ) +
    stat_summary(fun = median, geom = "point", size = 1.2, color = "black") +
    scale_fill_manual(values = c("ND" = "#71BDB6", "HFD" = "#DD5511")) +
    labs(title = title, y = ylab, x = NULL) +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = "right")
}

# ---- 主函数：一个 LR 对 → Sender/Receiver 上下拼接 ----
plot_lr_half_violin <- function(tag, lig_vec, rec_vec,
                                send_obj = senders, recv_obj = receivers,
                                outdir = "./LR_violin_half") {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  lig_gene <- gene_present(send_obj, lig_vec)
  rec_gene <- gene_present(recv_obj, rec_vec)
  
  if (is.null(lig_gene) && is.null(rec_gene)) {
    message("Skip [", tag, "]: neither ligand nor receptor present.")
    return(invisible(NULL))
  }
  
  p_sender <- if (!is.null(lig_gene)) {
    df <- prep_violin_df(send_obj, lig_gene, SENDER_LEVELS)
    half_violin_plot(df, paste0(tag, " | Sender (", lig_gene, ")"),
                     paste0(lig_gene, " expression"))
  } else NULL
  
  p_receiver <- if (!is.null(rec_gene)) {
    df <- prep_violin_df(recv_obj, rec_gene, RECEIVER_LEVELS)
    half_violin_plot(df, paste0(tag, " | Receiver (", rec_gene, ")"),
                     paste0(rec_gene, " expression"))
  } else NULL
  
  combo <- if (!is.null(p_sender) && !is.null(p_receiver)) {
    (p_sender / p_receiver) + patchwork::plot_annotation(
      title = paste0("Ligand–Receptor: ", tag)
    )
  } else if (!is.null(p_sender)) {
    p_sender + ggtitle(paste0("Ligand–Receptor: ", tag, " (Sender only)"))
  } else {
    p_receiver + ggtitle(paste0("Ligand–Receptor: ", tag, " (Receiver only)"))
  }
  
  safe_tag <- gsub("[^A-Za-z0-9]+", "_", tag); safe_tag <- gsub("^_+|_+$", "", safe_tag)
  outfile <- file.path(outdir, paste0(safe_tag, "_half.pdf"))
  ggsave(outfile, combo, width = 12, height = 8, dpi = 300)
  message("Saved: ", outfile)
}

# ---- 你的 LR 列表（示例；保持你已有的）----
lr_pairs <- list(
  "TENM4-ADGRL3" = c("Tenm4", "Adgrl3"),
  "APP-TNFRSF21" = c("App", "Tnfrsf21"),
  "LIPA-RORA"    = c("Lipa", "Rora"),
  "NRG2-ERBB4"   = c("Nrg2", "Erbb4")
)

# ---- 正确的遍历（按名字取 tag）----
for (tag in names(lr_pairs)) {
  pair <- lr_pairs[[tag]]
  plot_lr_half_violin(tag, lig_vec = pair[1], rec_vec = pair[2])
}

# ---- 仅受体端（方法2，最小改动：在末尾追加两行）----
plot_lr_half_violin("ESR1_only_receiver", lig_vec = character(0), rec_vec = "Esr1")
plot_lr_half_violin("ADGRL2_only_receiver", lig_vec = character(0), rec_vec = "Adgrl2")
