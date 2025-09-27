library(mapa)
library(org.Mm.eg.db)
library(dplyr)
library(purrr)
library(tibble)
load("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/2DEG/02-01-envVariables.RData")

deg_list$Basal$symbol <- rownames(deg_list$Basal)
deg_list$HormSens$symbol <- rownames(deg_list$HormSens)
deg_list$LumProg$symbol <- rownames(deg_list$LumProg)

to_variable_info <- function(
  df,
  id_type   = "symbol",
  order_col = "logFC",

  extra_cols = c("p.adjust", "p_val_adj", "avg_log2FC", "avg_log2fc", "pvalue", "p_val", "padj")
) {

  if (!(id_type %in% names(df))) {
    df <- tibble::rownames_to_column(df, var = id_type)
  }

  names(df) <- tolower(names(df))
  names(df) <- sub("^pct\\.1$", "pct_1", names(df))
  names(df) <- sub("^pct\\.2$", "pct_2", names(df))

  df <- df %>% dplyr::mutate(variable_id = .data[[tolower(id_type)]])

  vi <- convert_id(
    data = df %>% dplyr::select(variable_id, dplyr::all_of(tolower(id_type))),
    query_type   = "gene",
    from_id_type = tolower(id_type),
    organism     = "org.Mm.eg.db"
  )

  extra_cols <- tolower(extra_cols)
  keep_back <- c(
    "variable_id",
    tolower(order_col),
    intersect(c("fdr", "regulation", extra_cols), names(df))
  )

  vi <- vi %>% dplyr::left_join(dplyr::select(df, dplyr::any_of(keep_back)), by = "variable_id")

  vi
}

variable_info_list <- imap(deg_list, 
                           ~to_variable_info(.x, 
                                             id_type = "symbol", 
                                             order_col = "logFC"))

preprocess4mapa <- function(upordown, 
                            variable_info_list, 
                            celltype,
                            sim_cutoff){
  
  sim_cutoff <- sim_cutoff

  if (upordown == "up"){
    variable_info <- variable_info_list[[celltype]] %>%
      filter(avg_log2fc > 0.5) %>%
      filter(p_val < 0.05)
  } else if (upordown == "down"){
    variable_info <- variable_info_list[[celltype]] %>%
      filter(avg_log2fc < -0.5) %>%
      filter(p_val < 0.05)
  } else {
    stop("Argument 'upordown' must be 'up' or 'down'")
  }

  ora_result <- enrich_pathway(
    variable_info   = variable_info,
    query_type      = "gene",
    database        = c("go", "kegg", "reactome"),
    go.orgdb        = org.Mm.eg.db,
    go.keytype      = "ENTREZID",
    kegg.organism   = "mmu",
    kegg.keytype    = "kegg",
    reactome.organism = "mouse"
  )

  similarity_result <- merge_pathways(
    object = ora_result,
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

  functional_modules <- get_functional_modules(
    object = similarity_result,
    sim.cutoff = 0.4,
    cluster_method = "louvain"
  )

  return(functional_modules  = functional_modules)
}

basal_up_functional_modules <- 
  preprocess4mapa("up", 
                  variable_info_list, 
                  "Basal",
                  sim_cutoff = 0.4)

basal_down_functional_modules <- 
  preprocess4mapa("down", 
                  variable_info_list, 
                  "Basal",
                  sim_cutoff = 0.4)

lumprog_up_functional_modules <- 
  preprocess4mapa("up", 
                  variable_info_list, 
                  "LumProg",
                  sim_cutoff = 0.5)

lumprog_down_functional_modules <-
  preprocess4mapa("down", 
                  variable_info_list, 
                  "LumProg",
                  sim_cutoff = 0.4)

hormsens_up_functional_modules <-
  preprocess4mapa("up", 
                  variable_info_list, 
                  "HormSens",
                  sim_cutoff = 0.6)

hormsens_down_functional_modules <-
  preprocess4mapa("down", 
                  variable_info_list, 
                  "HormSens",
                  sim_cutoff = 0.5)

plot_similarity_network(
  object = basal_up_functional_modules,
  level = "functional_module",
  degree_cutoff = 3,
  text = TRUE
)

plot_similarity_network(
  object = basal_down_functional_modules,
  level = "functional_module",
  degree_cutoff = 1,
  text = TRUE
)

plot_similarity_network(
  object = lumprog_up_functional_modules,
  level = "functional_module",
  degree_cutoff = 3,
  text = TRUE
)

plot_similarity_network(
  object = lumprog_down_functional_modules,
  level = "functional_module",
  degree_cutoff = 1,
  text = TRUE
)

plot_similarity_network(
  object = hormsens_up_functional_modules,
  level = "functional_module",
  degree_cutoff = 3,
  text = TRUE
)

plot_similarity_network(
  object = hormsens_down_functional_modules,
  level = "functional_module",
  degree_cutoff = 3,
  text = TRUE
)
