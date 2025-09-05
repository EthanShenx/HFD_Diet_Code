setwd("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/9999prebutertyVsPuberty/formerVersion")

library(Seurat)
library(dplyr)
library(org.Mm.eg.db)
library(enrichplot)
library(ggplot2)
library(ReactomePA)
library(clusterProfiler)
library(tibble)

load("./9999-02-envVariables.RData")

###### ORA ######

ora_results_list <- list()

for (cell in cell_types) {
  deg <- deg_list[[cell]]

  up_genes <- rownames(deg[deg$regulation == "Up", ])
  down_genes <- rownames(deg[deg$regulation == "Down", ])
  
  up_genes  <- setdiff(up_genes, grep("^Rps|^Rpl", up_genes, value=TRUE))
  down_genes  <- setdiff(down_genes, grep("^Rps|^Rpl", down_genes, value=TRUE))

  up_entrez <- bitr(up_genes,
                    fromType = "SYMBOL",
                    toType   = "ENTREZID",
                    OrgDb    = org.Mm.eg.db)$ENTREZID
  down_entrez <- bitr(down_genes,
                      fromType = "SYMBOL",
                      toType   = "ENTREZID",
                      OrgDb    = org.Mm.eg.db)$ENTREZID

  # ---------------- GOBP ----------------
  up_go <- enrichGO(
    gene          = up_entrez,
    OrgDb         = org.Mm.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    readable      = TRUE
  )
  down_go <- enrichGO(
    gene          = down_entrez,
    OrgDb         = org.Mm.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    readable      = TRUE
  )

  # ---------------- KEGG ----------------
  up_kegg <- enrichKEGG(
    gene          = up_entrez,
    organism      = "mmu",
    pAdjustMethod = "BH"
  )
  down_kegg <- enrichKEGG(
    gene          = down_entrez,
    organism      = "mmu",
    pAdjustMethod = "BH"
  )

  # ---------------- Reactome ----------------
  up_reactome <- enrichPathway(
    gene          = up_entrez,
    organism      = "mouse",
    pAdjustMethod = "BH",
    readable      = TRUE
  )
  down_reactome <- enrichPathway(
    gene          = down_entrez,
    organism      = "mouse",
    pAdjustMethod = "BH",
    readable      = TRUE
  )

  # ---------------- Save ----------------
  write.csv(as.data.frame(up_go),
            paste0("ORA_results/", cell, "_Up_ORA_GO.csv"),
            row.names = FALSE)
  write.csv(as.data.frame(down_go),
            paste0("ORA_results/", cell, "_Down_ORA_GO.csv"),
            row.names = FALSE)

  write.csv(as.data.frame(up_kegg),
            paste0("ORA_results/", cell, "_Up_ORA_KEGG.csv"),
            row.names = FALSE)
  write.csv(as.data.frame(down_kegg),
            paste0("ORA_results/", cell, "_Down_ORA_KEGG.csv"),
            row.names = FALSE)

  write.csv(as.data.frame(up_reactome),
            paste0("ORA_results/", cell, "_Up_ORA_Reactome.csv"),
            row.names = FALSE)
  write.csv(as.data.frame(down_reactome),
            paste0("ORA_results/", cell, "_Down_ORA_Reactome.csv"),
            row.names = FALSE)

  ora_results_list[[paste0(cell, "_Up_GO")]]        <- up_go
  ora_results_list[[paste0(cell, "_Down_GO")]]      <- down_go
  ora_results_list[[paste0(cell, "_Up_KEGG")]]      <- up_kegg
  ora_results_list[[paste0(cell, "_Down_KEGG")]]    <- down_kegg
  ora_results_list[[paste0(cell, "_Up_Reactome")]]  <- up_reactome
  ora_results_list[[paste0(cell, "_Down_Reactome")]]<- down_reactome
}

###### GSEA ######

gsea_results_list <- list()
mm_gmt <- read.gmt("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/3Enrichment/mh.all.v2025.1.Mm.entrez.gmt")


for (cell in cell_types) {
  deg <- deg_list[[cell]]
  
  gene_vector <- deg$avg_log2FC
  gene_vector <- gene_vector[is.finite(gene_vector)]
  gene_vector <- sort(gene_vector, decreasing = TRUE)
  names(gene_vector) <- rownames(deg)
  

  gene_df <- bitr(names(gene_vector),
                  fromType = "SYMBOL",
                  toType   = "ENTREZID",
                  OrgDb    = org.Mm.eg.db)

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

  gsea_kegg <- gseKEGG(
    geneList      = gene_vector,
    organism      = "mmu",
    pAdjustMethod = "BH",
    verbose       = FALSE
  )

  gsea_reactome <- gsePathway(
    geneList      = gene_vector,
    organism      = "mouse",
    pAdjustMethod = "BH",
    verbose       = FALSE
  )

  gsea_hallmark <- GSEA(
    geneList      = gene_vector,
    TERM2GENE     = mm_gmt,
    verbose       = FALSE
  )
  
  gsea_hallmark@result$count <-
    purrr::map_int(
      gsea_hallmark@result$core_enrichment,
      function(x) {
        length(stringr::str_split(x, pattern = "/")[[1]])
      }
    )

  write.csv(as.data.frame(gsea_go),
            file = paste0("GSEA_results/", cell, "_GSEA_GO.csv"),
            row.names = FALSE)
  write.csv(as.data.frame(gsea_kegg),
            file = paste0("GSEA_results/", cell, "_GSEA_KEGG.csv"),
            row.names = FALSE)
  write.csv(as.data.frame(gsea_reactome),
            file = paste0("GSEA_results/", cell, "_GSEA_Reactome.csv"),
            row.names = FALSE)
  write.csv(as.data.frame(gsea_hallmark),
            file = paste0("GSEA_results/", cell, "_GSEA_Hallmark.csv"),
            row.names = FALSE)

  gsea_results_list[[paste0(cell, "_GO")]]       <- gsea_go
  gsea_results_list[[paste0(cell, "_KEGG")]]     <- gsea_kegg
  gsea_results_list[[paste0(cell, "_Reactome")]] <- gsea_reactome
  gsea_results_list[[paste0(cell, "_Hallmark")]] <- gsea_hallmark
}
