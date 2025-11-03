library(mapa)
library(org.Mm.eg.db)
library(dplyr)
library(purrr)
library(tibble)
load("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/9999prebutertyVsPuberty/formerVersion/9999-02-envVariables.RData")

deg_list$Basal$symbol <- rownames(deg_list$Basal)
deg_list$HormSens$symbol <- rownames(deg_list$HormSens)
deg_list$LumProg$symbol <- rownames(deg_list$LumProg)

to_variable_info <- function(df, id_type = "symbol", order_col = "logFC") {
  if (!(id_type %in% names(df))) {
    df <- tibble::rownames_to_column(df, var = id_type)
  }

  names(df) <- tolower(names(df))

  names(df) <- sub("^pct\\.1$", "pct_1", names(df))
  names(df) <- sub("^pct\\.2$", "pct_2", names(df))

  df <- df %>% mutate(variable_id = .data[[id_type]])

  vi <- convert_id(
    data = df %>% dplyr::select(variable_id, all_of(id_type)),
    query_type   = "gene",
    from_id_type = id_type,
    organism     = "org.Mm.eg.db"
  )

  keep_back <- c("variable_id", tolower(order_col), intersect(c("fdr","regulation"), names(df)))
  vi <- vi %>% left_join(dplyr::select(df, any_of(keep_back)), by = "variable_id")

  vi
}

variable_info_list <- imap(deg_list, 
                           ~to_variable_info(.x, 
                                             id_type = "symbol", 
                                             order_col = "logFC"))

basal_down_info <- variable_info_list$Basal %>%
  filter(fdr < 0.05) %>%
  filter(logfc < -0.5)

hs_down_info <- variable_info_list$HormSens %>%
  filter(fdr < 0.05) %>%
  filter(logfc < -0.5)

lp_down_info <- variable_info_list$LumProg %>%
  filter(fdr < 0.05) %>%
  filter(logfc < -0.5)

basal_up_info <- variable_info_list$Basal %>%
  filter(fdr < 0.05) %>%
  filter(logfc > 0.5)

hs_up_info <- variable_info_list$HormSens %>%
  filter(fdr < 0.05) %>%
  filter(logfc > 0.5)

lp_up_info <- variable_info_list$LumProg %>%
  filter(fdr < 0.05) %>%
  filter(logfc > 0.5)

basal_ora_down <- enrich_pathway(
  variable_info   = basal_down_info,
  query_type      = "gene",
  database        = c("go", "kegg", "reactome"),
  go.orgdb        = org.Mm.eg.db,
  go.keytype      = "ENTREZID",
  kegg.organism   = "mmu",
  kegg.keytype    = "kegg",
  reactome.organism = "mouse"
)

hs_ora_down <- enrich_pathway(
  variable_info   = hs_down_info,
  query_type      = "gene",
  database        = c("go", "kegg", "reactome"),
  go.orgdb        = org.Mm.eg.db,
  go.keytype      = "ENTREZID",
  kegg.organism   = "mmu",
  kegg.keytype    = "kegg",
  reactome.organism = "mouse"
)
lp_ora_down <- enrich_pathway(
  variable_info   = lp_down_info,
  query_type      = "gene",
  database        = c("go", "kegg", "reactome"),
  go.orgdb        = org.Mm.eg.db,
  go.keytype      = "ENTREZID",
  kegg.organism   = "mmu",
  kegg.keytype    = "kegg",
  reactome.organism = "mouse"
)
basal_ora_up <- enrich_pathway(
  variable_info   = basal_up_info,
  query_type      = "gene",
  database        = c("go", "kegg", "reactome"),
  go.orgdb        = org.Mm.eg.db,
  go.keytype      = "ENTREZID",
  kegg.organism   = "mmu",
  kegg.keytype    = "kegg",
  reactome.organism = "mouse"
)
hs_ora_up <- enrich_pathway(
  variable_info   = hs_up_info,
  query_type      = "gene",
  database        = c("go", "kegg", "reactome"),
  go.orgdb        = org.Mm.eg.db,
  go.keytype      = "ENTREZID",
  kegg.organism   = "mmu",
  kegg.keytype    = "kegg",
  reactome.organism = "mouse"
)
lp_ora_up <- enrich_pathway(
  variable_info   = lp_up_info,
  query_type      = "gene",
  database        = c("go", "kegg", "reactome"),
  go.orgdb        = org.Mm.eg.db,
  go.keytype      = "ENTREZID",
  kegg.organism   = "mmu",
  kegg.keytype    = "kegg",
  reactome.organism = "mouse"
)

basal_down_similarity_result <- 
  merge_pathways(
    object = basal_ora_down,
    database = c("go", "kegg", "reactome"),
    p.adjust.cutoff.go = 0.05,
    p.adjust.cutoff.kegg = 0.05,
    p.adjust.cutoff.reactome = 0.05,
    count.cutoff.go = 5,
    count.cutoff.kegg = 5,
    count.cutoff.reactome = 5,
    measure.method.go = "Sim_XGraSM_2013",
    go.orgdb = "org.Mm.eg.db",               # Required for GO analysis
    measure.method.kegg = "jaccard",        # Gene overlap similarity
    measure.method.reactome = "jaccard"     # Gene overlap similarity
  )

basal_up_similarity_result <- 
  merge_pathways(
    object = basal_ora_up,
    database = c("go", "kegg", "reactome"),
    p.adjust.cutoff.go = 0.05,
    p.adjust.cutoff.kegg = 0.05,
    p.adjust.cutoff.reactome = 0.05,
    count.cutoff.go = 5,
    count.cutoff.kegg = 5,
    count.cutoff.reactome = 5,
    measure.method.go = "Sim_XGraSM_2013",  # GO semantic similarity
    go.orgdb = "org.Mm.eg.db",               # Required for GO analysis
    measure.method.kegg = "jaccard",        # Gene overlap similarity
    measure.method.reactome = "jaccard"     # Gene overlap similarity
  )

hs_down_similarity_result <-
  merge_pathways(
    object = hs_ora_down,
    database = c("go", "kegg", "reactome"),
    p.adjust.cutoff.go = 0.05,
    p.adjust.cutoff.kegg = 0.05,
    p.adjust.cutoff.reactome = 0.05,
    count.cutoff.go = 5,
    count.cutoff.kegg = 5,
    count.cutoff.reactome = 5,
    measure.method.go = "Sim_XGraSM_2013",  # GO semantic similarity
    go.orgdb = "org.Mm.eg.db",               # Required for GO analysis
    measure.method.kegg = "jaccard",        # Gene overlap similarity
    measure.method.reactome = "jaccard"     # Gene overlap similarity
  )

hs_up_similarity_result <-
  merge_pathways(
    object = hs_ora_up,
    database = c("go", "kegg", "reactome"),
    p.adjust.cutoff.go = 0.05,
    p.adjust.cutoff.kegg = 0.05,
    p.adjust.cutoff.reactome = 0.05,
    count.cutoff.go = 5,
    count.cutoff.kegg = 5,
    count.cutoff.reactome = 5,
    measure.method.go = "Sim_XGraSM_2013",  # GO semantic similarity
    go.orgdb = "org.Mm.eg.db",               # Required for GO analysis
    measure.method.kegg = "jaccard",        # Gene overlap similarity
    measure.method.reactome = "jaccard"     # Gene overlap similarity
  )

lp_down_similarity_result <-
  merge_pathways(
    object = lp_ora_down,
    database = c("go", "kegg", "reactome"),
    p.adjust.cutoff.go = 0.05,
    p.adjust.cutoff.kegg = 0.05,
    p.adjust.cutoff.reactome = 0.05,
    count.cutoff.go = 5,
    count.cutoff.kegg = 5,
    count.cutoff.reactome = 5,
    measure.method.go = "Sim_XGraSM_2013",  # GO semantic similarity
    go.orgdb = "org.Mm.eg.db",               # Required for GO analysis
    measure.method.kegg = "jaccard",        # Gene overlap similarity
    measure.method.reactome = "jaccard"     # Gene overlap similarity
  )

lp_up_similarity_result <-
  merge_pathways(
    object = lp_ora_up,
    database = c("go", "kegg", "reactome"),
    p.adjust.cutoff.go = 0.05,
    p.adjust.cutoff.kegg = 0.05,
    p.adjust.cutoff.reactome = 0.05,
    count.cutoff.go = 5,
    count.cutoff.kegg = 5,
    count.cutoff.reactome = 5,
    measure.method.go = "Sim_XGraSM_2013",  # GO semantic similarity
    go.orgdb = "org.Mm.eg.db",               # Required for GO analysis
    measure.method.kegg = "jaccard",        # Gene overlap similarity
    measure.method.reactome = "jaccard"     # Gene overlap similarity
  )

basal_up_functional_modules <-
  get_functional_modules(
    object = basal_up_similarity_result,
    sim.cutoff = 0.45,
    cluster_method = "louvain"
  )

hs_up_functional_modules <-
  get_functional_modules(
    object = hs_up_similarity_result,
    sim.cutoff = 0.45,
    cluster_method = "louvain"
  )

lp_up_functional_modules <-
  get_functional_modules(
    object = lp_up_similarity_result,
    sim.cutoff = 0.45,
    cluster_method = "louvain"
  )

plot_similarity_network(
  object = basal_up_functional_modules,
  level = "functional_module",
  degree_cutoff = 3,
  text = TRUE
)

 plot_similarity_network(
  object = lp_up_functional_modules,
  level = "functional_module",
  degree_cutoff = 4,
  text = TRUE
)

plot_similarity_network(
  object = hs_up_functional_modules,
  level = "functional_module",
  degree_cutoff = 4,
  text = TRUE
)

write.csv(basal_up_functional_modules@merged_module$functional_module_result,
          file = "/Users/coellearth/Desktop/HFD_Paper/Yixiang/Basal_SI.csv")
write.csv(lp_up_functional_modules@merged_module$functional_module_result,
          file = "/Users/coellearth/Desktop/HFD_Paper/Yixiang/LumProg_SI.csv")
write.csv(hs_up_functional_modules@merged_module$functional_module_result,
          file = "/Users/coellearth/Desktop/HFD_Paper/Yixiang/HormSens_SI.csv")

