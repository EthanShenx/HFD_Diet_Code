library(clusterProfiler)
library(org.Mm.eg.db)
library(ReactomePA)
library(dplyr)

# ===== set path =====
indir  <- "D:/data/23BMI/ND_HFD_MG_snRNAseq/formerVersion/DEG_results"
outdir <- file.path(indir, "ORA_results_meant")
dir.create(outdir, showWarnings = FALSE)

# ===== find all *_UpDown.csv / *_DownUp.csv files =====
deg_files <- list.files(indir, pattern = ".*_(UpDown|DownUp)\\.csv$", full.names = TRUE)

# ===== 批量 ORA =====
for (f in deg_files) {
  nm <- gsub("\\.csv$", "", basename(f))   # Basal_UpDown 之类
  
  df <- read.csv(f)
  if (!"gene" %in% colnames(df)) {
    warning("lack of gene column: ", f)
    next
  }
  
  genes <- unique(df$gene)
  
  # remove ribosome genes
  genes <- setdiff(genes, grep("^Rps|^Rpl", genes, value = TRUE))
  
  # SYMBOL -> ENTREZ
  entrez <- bitr(genes,
                 fromType = "SYMBOL",
                 toType   = "ENTREZID",
                 OrgDb    = org.Mm.eg.db)$ENTREZID
  
  if (length(entrez) == 0) {
    warning("failed to map gene: ", nm)
    next
  }
  
  # --- GOBP ---
  go_res <- enrichGO(entrez, OrgDb = org.Mm.eg.db, ont = "BP",
                     pAdjustMethod = "BH", readable = TRUE)
  
  # --- KEGG ---
  kegg_res <- enrichKEGG(entrez, organism = "mmu", pAdjustMethod = "BH")
  
  # --- Reactome ---
  reactome_res <- enrichPathway(entrez, organism = "mouse",
                                pAdjustMethod = "BH", readable = TRUE)
  
  # save
  write.csv(as.data.frame(go_res),
            file.path(outdir, paste0(nm, "_GO.csv")), row.names = FALSE)
  write.csv(as.data.frame(kegg_res),
            file.path(outdir, paste0(nm, "_KEGG.csv")), row.names = FALSE)
  write.csv(as.data.frame(reactome_res),
            file.path(outdir, paste0(nm, "_Reactome.csv")), row.names = FALSE)
  
  message("ORA completed: ", nm)
}
