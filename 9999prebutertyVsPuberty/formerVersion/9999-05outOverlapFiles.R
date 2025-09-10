library(Seurat)
library(dplyr)
library(org.Mm.eg.db)
library(enrichplot)
library(ggplot2)
library(ReactomePA)
library(clusterProfiler)
library(tibble)

# ===== Read in and process data str =====
cell_types <- c("Basal", "LumProg", "HormSens")

setwd("D:/data/23BMI/ND_HFD_MG_snRNAseq/formerVersion/DEG_results/allStages_CellTypes_pseudoBulk")
prepub_files <- list.files(pattern = "*_DEGs_pseudobulk\\.csv$")
prepub_deg <- lapply(prepub_files, function(f) {
  df <- read.csv(f, row.names = 1, check.names = FALSE)
  df
})
names(prepub_deg) <- gsub("^prepub_|_DEGs_pseudobulk\\.csv$", "", prepub_files)

setwd("D:/data/23BMI/ND_HFD_MG_snRNAseq/formerVersion/DEG_results")
hfd_files <- list.files(pattern = "*_DEGs\\.csv$")
hfd_deg <- lapply(hfd_files, function(f) {
  df <- read.csv(f, row.names = 1, check.names = FALSE)
  df
})
names(hfd_deg) <- gsub("^HFD_|_DEGs\\.csv$", "", hfd_files)

# ===== Get dysregulated genes =====
# Note: I did not set a bar to filter p.adjust for 2 reasons:
# 1. Areg (prior knowledge)
# 2. Different level of adjustment for false positive will cause # different level of false negative

process_deg <- function(deg_list, suffix) {
  lapply(deg_list, function(df) {
    df <- as.data.frame(df)
    colnames(df) <- trimws(colnames(df))
    df$gene <- rownames(df)
    list(
      up = df %>%
        dplyr::filter(.data$avg_log2FC > 1.5) %>%
        rename_with(~ paste0(.x, "_", suffix), .cols = !any_of("gene")),
      down = df %>%
        dplyr::filter(.data$avg_log2FC < -1.5) %>%
        rename_with(~ paste0(.x, "_", suffix), .cols = !any_of("gene"))
    )
  })
}


merge_opposites <- function(pub, diet) {
  mapply(function(p, d) {
    list(
      up_down = merge(p$up, d$down, by = "gene") %>%
        dplyr::select(!matches("^cell_type|^pct")),
      down_up = merge(p$down, d$up, by = "gene") %>%
        dplyr::select(!matches("^cell_type|^pct"))
    )
  }, pub, diet, SIMPLIFY = FALSE)
}

cell_types <- c("Basal", "LumProg", "HormSens")
pub_deg <- process_deg(prepub_deg[cell_types], "pub")
diet_deg <- process_deg(hfd_deg[cell_types], "diet")
results <- merge_opposites(pub_deg, diet_deg)

Basal_Meant_Up_but_Down <- results$Basal$up_down
Basal_Meant_Down_but_Up <- results$Basal$down_up
LumProg_Meant_Up_but_Down <- results$LumProg$up_down
LumProg_Meant_Down_but_Up <- results$LumProg$down_up
HormSens_Meant_Up_but_Down <- results$HormSens$up_down
HormSens_Meant_Down_but_Up <- results$HormSens$down_up

# ===== Visualization 1: upset plot =====

library(UpSetR)

gene_lists <- list(
  Basal_UpDown   = Basal_Meant_Up_but_Down$gene,
  Basal_DownUp   = Basal_Meant_Down_but_Up$gene,
  LumProg_UpDown = LumProg_Meant_Up_but_Down$gene,
  LumProg_DownUp = LumProg_Meant_Down_but_Up$gene,
  HormSens_UpDown = HormSens_Meant_Up_but_Down$gene,
  HormSens_DownUp = HormSens_Meant_Down_but_Up$gene
)

UpSetR::upset(
  fromList(gene_lists),
  nsets   = 6,
  nintersects = NA,
  order.by    = "freq",
  mb.ratio    = c(0.6, 0.4),
  shade.color = "#ecf0f1"
)

outdir <- "D:/data/23BMI/ND_HFD_MG_snRNAseq/formerVersion/DEG_results"

# 建立列表统一管理
gene_results <- list(
  Basal_UpDown   = Basal_Meant_Up_but_Down,
  Basal_DownUp   = Basal_Meant_Down_but_Up,
  LumProg_UpDown = LumProg_Meant_Up_but_Down,
  LumProg_DownUp = LumProg_Meant_Down_but_Up,
  HormSens_UpDown = HormSens_Meant_Up_but_Down,
  HormSens_DownUp = HormSens_Meant_Down_but_Up
)

# 循环写出 CSV（无任何处理）
for (nm in names(gene_results)) {
  outfile <- file.path(outdir, paste0(nm, ".csv"))
  write.csv(gene_results[[nm]], outfile, row.names = FALSE)
  message("已保存: ", outfile)
}

