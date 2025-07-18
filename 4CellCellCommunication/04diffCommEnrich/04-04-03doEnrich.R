setwd("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/4CellCellCommunication/04diffCommEnrich")

library(dplyr)
library(CellChat)
library(clusterProfiler)
library(GSEABase)

ND <- read.csv("/Users/coellearth/Desktop/HFD_Paper/Fig2_PhenontypicOverview/diffAllLRPairLIANA/LIANARawData/ND_LIANA_res.csv")
HFD <- read.csv("/Users/coellearth/Desktop/HFD_Paper/Fig2_PhenontypicOverview/diffAllLRPairLIANA/LIANARawData/HFD_LIANA_res.csv")

dat <- inner_join(
  ND %>%
    dplyr::select(
      source,
      target,
      ligand.complex,
      receptor.complex,
      mean_rank
    ) %>%
    dplyr::rename(mean_rank_ND = mean_rank),
  HFD %>%
    select(
      source,
      target,
      ligand.complex,
      receptor.complex,
      mean_rank
    ) %>%
    dplyr::rename(mean_rank_HFD = mean_rank),
  by = c("source", "target", "ligand.complex", "receptor.complex")
) %>%
  mutate(
    score = mean_rank_HFD - mean_rank_ND
  ) %>%
  dplyr::rename(
    ligand = ligand.complex,
    receptor = receptor.complex
  ) %>%
  dplyr::select(
    - mean_rank_ND,
    - mean_rank_HFD
  )

Up_reg <- dat %>%
  filter(score > 0) %>%
  arrange(desc(score))

Down_reg <- dat %>%
  filter(score < 0) %>%
  mutate(score = abs(score)) %>%
  arrange(desc(score))

up_genes <- unique(Up_reg$ligand, Up_reg$receptor)
down_genes <- unique(Down_reg$ligand, Down_reg$receptor)

gmt_to_df <- function(gmt_path) {
  gmt <- getGmt(gmt_path)
  gmt <- gmt[sapply(gmt, 
                    function(gs) 
                      length(geneIds(gs)) > 0)]
  term2gene <- do.call(rbind, lapply(gmt, function(gs) {
    data.frame(term = setName(gs), 
               gene = geneIds(gs), 
               stringsAsFactors = FALSE)
  }))
  
  return(term2gene)
}

m3_term2gene <- gmt_to_df("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/4CellCellCommunication/04diffCommEnrich/processedResGMT/m3_trimmed.gmt")
m5_term2gene <- gmt_to_df("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/4CellCellCommunication/04diffCommEnrich/processedResGMT/m5_trimmed.gmt")

Up_m3_enrich <- enricher(up_genes, 
                         TERM2GENE = m3_term2gene)
Up_m3_enrich_result <- as.data.frame(Up_m3_enrich@result)
Up_m3_enrich_result <- Up_m3_enrich_result |>
  filter(p.adjust < 0.05) |>
  arrange(desc(FoldEnrichment))

Up_m5_enrich <- enricher(up_genes, 
                         TERM2GENE = m5_term2gene)
Up_m5_enrich_result <- as.data.frame(Up_m5_enrich@result)
Up_m5_enrich_result <- Up_m5_enrich_result |>
  filter(p.adjust < 0.05) |>
  arrange(desc(FoldEnrichment))

