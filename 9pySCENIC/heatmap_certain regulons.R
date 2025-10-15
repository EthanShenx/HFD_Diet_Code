library(tidyverse)
library(pheatmap)
library(RColorBrewer)

regulons <- read.csv("regulons.csv", check.names = TRUE, stringsAsFactors = FALSE)
head(regulons)

areg_regulons <- regulons %>%
  filter(str_detect(TargetGenes, "\\bAreg\\b")) %>%
  select(TF)

# TF name
areg_TFs <- areg_regulons$TF
cat(length(areg_TFs), "Areg related regulon：\n")
print(areg_TFs)

areg_related_keywords <- c("Areg", "Egfr", "Egr1", "Fos", "Jun", "Atf3", "Stat3", "Nr4a1", "Elf3")

# regulons
areg_rows_HFD <- rownames(avg_HFD)[str_detect(rownames(avg_HFD), paste(areg_related_keywords, collapse = "|"))]
areg_rows_ND  <- rownames(avg_ND)[str_detect(rownames(avg_ND), paste(areg_related_keywords, collapse = "|"))]

# _sub
avg_HFD_sub <- avg_HFD[areg_rows_HFD, , drop = FALSE]
avg_ND_sub  <- avg_ND[areg_rows_ND, , drop = FALSE]

# average same TF
avg_HFD_sub <- avg_HFD_sub %>%
  as.data.frame() %>%
  rownames_to_column("TF") %>%
  group_by(TF) %>%
  summarise(across(everything(), mean)) %>%
  column_to_rownames("TF")

avg_ND_sub <- avg_ND_sub %>%
  as.data.frame() %>%
  rownames_to_column("TF") %>%
  group_by(TF) %>%
  summarise(across(everything(), mean)) %>%
  column_to_rownames("TF")

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
