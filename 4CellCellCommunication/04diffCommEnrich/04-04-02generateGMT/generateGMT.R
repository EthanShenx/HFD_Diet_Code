library(GSEABase)

data_dir <- "/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/4CellCellCommunication/04diffCommEnrich"

m3_gmt <- file.path(data_dir, "m3.all.v2025.1.Mm.symbols.gmt")
m5_gmt <- file.path(data_dir, "m5.go.bp.v2025.1.Mm.symbols.gmt")
output_m3_csv <- file.path(data_dir, "m3_mouse_intersection.csv")
output_m5_csv <- file.path(data_dir, "m5_mouse_intersection.csv")
output_m3_gmt <- file.path(data_dir, "m3_mouse_intersected.gmt")
output_m5_gmt <- file.path(data_dir, "m5_mouse_intersected.gmt")

data_dir <- "/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/4CellCellCommunication/04diffCommEnrich/mouse"

csv_files <- list.files(data_dir, pattern = "_mouse\\.csv$", full.names = TRUE)
all_genes <- unique(unlist(lapply(csv_files, function(f) {
  df <- read.csv(f)
  genes <- c(df$source_genesymbol, df$target_genesymbol)
  genes[!is.na(genes)]
})))

m3 <- getGmt(m3_gmt)
m5 <- getGmt(m5_gmt)

trim_pathways_by_genes <- function(gmt_obj, gene_pool) {
  trimmed <- lapply(gmt_obj, function(g) {
    genes <- intersect(geneIds(g), gene_pool)
    GeneSet(genes,
            setName = setName(g),
            shortDescription = description(g),
            geneIdType = geneIdType(g),
            collectionType = collectionType(g),
            organism = organism(g))
  })
  GeneSetCollection(trimmed)
}

m3_trimmed <- trim_pathways_by_genes(m3, all_genes)
m5_trimmed <- trim_pathways_by_genes(m5, all_genes)

writeGMT <- function(gsc, file) {
  con <- file(file, "w")
  for (gs in gsc) {
    line <- paste(setName(gs), description(gs), paste(geneIds(gs), collapse = "\t"), sep = "\t")
    writeLines(line, con)
  }
  close(con)
  cat("Saved：", file, "\n")
}

data_dir <- "/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/4CellCellCommunication/04diffCommEnrich/processedResGMT"

writeGMT(m3_trimmed, file.path(data_dir, "m3_trimmed.gmt"))
writeGMT(m5_trimmed, file.path(data_dir, "m5_trimmed.gmt"))
