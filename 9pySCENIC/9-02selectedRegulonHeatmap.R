library(tidyverse)
library(pheatmap)
library(RColorBrewer)

regulons <- read.csv("regulons.csv", check.names = TRUE, stringsAsFactors = FALSE)
head(regulons)

genes <- c("Areg", "Lpl", "Pparg", "Dpp4", "Pdgfra", "Pi16", "Fap", "Lama2", "Mmp3", "Mmp9")

regulons_of_interest <- regulons %>%
  filter(str_detect(TargetGenes, str_c("\\b(", str_c(genes, collapse = "|"), ")\\b"))) %>%
  select(TF)

regulon_TFs <- regulons_of_interest$TF

pattern <- paste(regulon_TFs, collapse = "|")

areg_rows_HFD <- rownames(avg_HFD)[str_detect(rownames(avg_HFD), pattern)]
areg_rows_ND  <- rownames(avg_ND)[str_detect(rownames(avg_ND), pattern)]

# _sub
avg_HFD_sub <- avg_HFD[areg_rows_HFD, , drop = FALSE]
avg_ND_sub  <- avg_ND[areg_rows_ND, , drop = FALSE]

# # average same TF
# avg_HFD_sub <- avg_HFD_sub %>%
#   as.data.frame() %>%
#   rownames_to_column("TF") %>%
#   group_by(TF) %>%
#   summarise(across(everything(), mean)) %>%
#   column_to_rownames("TF")
# 
# avg_ND_sub <- avg_ND_sub %>%
#   as.data.frame() %>%
#   rownames_to_column("TF") %>%
#   group_by(TF) %>%
#   summarise(across(everything(), mean)) %>%
#   column_to_rownames("TF")

pheatmap(avg_HFD_sub,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
         scale = "row",
         main = "HFD - Areg-related regulons")

pheatmap(avg_ND_sub,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
         scale = "row",
         main = "ND - Areg-related regulons")

