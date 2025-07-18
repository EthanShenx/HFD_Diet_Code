library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(ggplot2)
library(ReactomePA)

setwd("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/3Enrichment")
load("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/2DEG/02-01-envVariables.RData")
All <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_all.rds")

###### GSEA core ######

mm_gmt <- read.gmt("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/3Enrichment/mh.all.v2025.1.Mm.entrez.gmt")

gsea_results_list <- list()
gsea_res_df_list <- list()

ribo_pattern <- "^Rps|^Rpl"

cell_types <- All@meta.data$cell_type %>%
  unique()

for (cell in cell_types) {
  deg <- deg_list[[cell]]
  gene_vector <- deg$avg_log2FC
  gene_vector <- gene_vector[is.finite(gene_vector)]
  names(gene_vector) <- rownames(deg)
  gene_vector <- sort(gene_vector, decreasing = TRUE)
  gene_vector <- gene_vector[!grepl(ribo_pattern, names(gene_vector))]
  
  gene_df <- bitr(
    names(gene_vector),
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = org.Mm.eg.db
  )
  
  gene_vector <- gene_vector[gene_df$SYMBOL]
  names(gene_vector) <- gene_df$ENTREZID
  gene_vector <- sort(gene_vector, decreasing = TRUE)
  
  gsea_go <- gseGO(
    geneList      = gene_vector,
    OrgDb         = org.Mm.eg.db,
    ont           = "BP",
    keyType       = "ENTREZID",
    pAdjustMethod = "BH",
    verbose       = FALSE
  )
  
  gsea_go_df <- 
    as.data.frame(gsea_go@result)
  
  suffixes <- c("GO")
  results <- list(gsea_go)
  df_results <- list(gsea_go_df)
  
  for (i in seq_along(suffixes)) {
    name <- paste0(cell, "_", suffixes[i])
    name_df <- paste0(cell, "_", suffixes[i], "_df")
    gsea_results_list[[name]] <- results[[i]]
    gsea_res_df_list[[name_df]] <- df_results[[i]]
    assign(name, results[[i]], envir = .GlobalEnv)
    assign(name_df, df_results[[i]], envir = .GlobalEnv)
  }
}

rm(results)
rm(df_results)
rm(All)
rm(mm_gmt)
rm(ribo_pattern)
rm(gsea_go)