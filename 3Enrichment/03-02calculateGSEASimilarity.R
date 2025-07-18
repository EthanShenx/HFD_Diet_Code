setwd("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/3Enrichment")
source("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/3Enrichment/03-99GSEAforSimilarity.R")

library(simplifyEnrichment)
library(igraph)
library(ggraph)
library(ggdendro)
library(RColorBrewer)
library(pheatmap)
library(DOSE)
library(dplyr)
library(clusterProfiler)
library(GOSemSim)
library(RCy3)

mmGO <- godata("org.Mm.eg.db", ont = "BP", computeIC = TRUE)
sim_res_list <- list()

for (cell in cell_types) {
  # =========== Simplify GO pathways ===========

  sim_res <- clusterProfiler::simplify(
    gsea_results_list[[paste0(cell, "_GO")]],
    cutoff = 0.7,
    by = "p.adjust",
    select_fun = min,
    measure = "Wang",
    semData = NULL
  )

  gosim_res <- data.frame(sim_res@result)
  sim_res_list[[paste0(cell, "_GO_sim")]] <- sim_res

  assign(paste0(cell, "_GO_sim"), gosim_res, envir = .GlobalEnv)

  setwd("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/3Enrichment/afterGSEASimplification")

  write.table(gosim_res,
    file = paste0(cell, "GSEAsim.txt"),
    sep = "\t",
    row.names = F,
    quote = F
  )

  # =========== Calculate similarity: Jaccard ===========

  pathway_genes <- lapply(sim_res@result$ID, function(go_id) {
    genes <- strsplit(sim_res@result[sim_res@result$ID == go_id, "core_enrichment"], "/")[[1]]
    unique(genes)
  })

  names(pathway_genes) <- sim_res@result$ID
  jaccard_matrix <- matrix(
    nrow = length(pathway_genes),
    ncol = length(pathway_genes)
  )
  rownames(jaccard_matrix) <- names(pathway_genes)
  colnames(jaccard_matrix) <- names(pathway_genes)

  for (i in 1:(length(pathway_genes) - 1)) {
    for (j in (i + 1):length(pathway_genes)) {
      a <- length(intersect(pathway_genes[[i]], pathway_genes[[j]]))
      b <- length(union(pathway_genes[[i]], pathway_genes[[j]]))
      jaccard_matrix[i, j] <- a / b
    }
  }

  edges <- as.data.frame(
    which(
      jaccard_matrix > 0.1 &
        upper.tri(jaccard_matrix),
      arr.ind = TRUE
    )
  ) %>%
    mutate(
      from = rownames(jaccard_matrix)[row],
      to = colnames(jaccard_matrix)[col],
      weight = jaccard_matrix[cbind(row, col)]
    ) %>%
    select(from, to, weight) %>%
    filter(weight > 0.25)

  nodes <- data.frame(ID = unique(c(edges$from, edges$to)))

  nodes <- merge(
    nodes,
    sim_res@result[, c("ID", 
                       "Description", 
                       "p.adjust",
                       "NES")],
    by    = "ID",
    all.x = TRUE
  )

  nodes <- nodes |>
    filter(p.adjust < 0.05)
  nodes$negLogP <- -log10(nodes$p.adjust)

  colnames(nodes)[colnames(nodes) == "ID"] <- "id"
  colnames(edges)[colnames(edges) == "from"] <- "source"
  colnames(edges)[colnames(edges) == "to"] <- "target"

  assign(paste0(cell, "JaccardEdges"),
    edges,
    envir = .GlobalEnv
  )

  assign(paste0(cell, "JaccardNodes"),
    nodes,
    envir = .GlobalEnv
  )

  setwd("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/3Enrichment/computedJaccardSim")

  write.table(edges,
    file = paste0(cell, "NetworkEdges.txt"),
    sep = "\t",
    row.names = F,
    quote = F
  )

  write.table(nodes,
    file = paste0(cell, "NetworkNodes.txt"),
    sep = "\t",
    row.names = F,
    quote = F
  )

  # ============ Save to Cytoscape ==============

  createNetworkFromDataFrames(
    nodes = nodes,
    edges = edges,
    title = paste0(cell, "_Jaccard_network"),
    collection = "GSEA_similarity",
    node.id.list = "id",
    source.id.list = "source",
    target.id.list = "target"
  )

  # ============ Set Cytoscape style =============
  cat(paste0(
    "=== Now setting Cytoscape style for ",
    cell,
    "... ==="
  ))

  style.name <- paste0(cell, "_Jaccard_Style")
  
  neg.range <- range(nodes$negLogP, na.rm = TRUE)
  wt.range <- range(edges$weight, na.rm = TRUE)
  nes.vec <- nodes$NES
  max.abs.nes <- max(abs(nes.vec), na.rm = TRUE)
  
  nodeSizeMap <- mapVisualProperty(
    visual.prop         = "NODE_SIZE",
    table.column        = "negLogP",
    mapping.type        = "c",
    table.column.values = neg.range,
    visual.prop.values  = c(10, 30)
  )

  edgeWidthMap <- mapVisualProperty(
    visual.prop         = "EDGE_WIDTH",
    table.column        = "weight",
    mapping.type        = "c",
    table.column.values = wt.range,
    visual.prop.values  = c(1, 5)
  )
  
  colorMap <- mapVisualProperty(
    visual.prop         = "NODE_FILL_COLOR",
    table.column        = "NES",
    mapping.type        = "c",
    table.column.values = c(-max.abs.nes, 
                            0, 
                            max.abs.nes),
    visual.prop.values  = c("#276419", 
                            "#F7F7F7", 
                            "#C51B7D")
  )

  createVisualStyle(
    style.name = style.name,
    defaults   = list(NODE_SHAPE  = "ELLIPSE",
                      EDGE_WIDTH  = 1),
    mappings   = list(nodeSizeMap,
                      edgeWidthMap,
                      colorMap)
  )

  setVisualStyle(style.name)

  # # =========== Calculate similarity: GOSemSim ===========
  #
  # sem_sim <- mgoSim(sim_res@result$ID,
  #                   sim_res@result$ID,
  #                   semData=mmGO,
  #                   measure="Wang",
  #                   combine=NULL)
  #
  # edges_sem <- as.data.frame(as.table(sem_sim)) %>%
  #   filter(Var1 != Var2) %>%
  #   rename(from=Var1, to=Var2, weight=Freq) %>%
  #   filter(weight > 0.3)
  #
  # assign(paste0(cell, "SemEdges"), edges_sem, envir = .GlobalEnv)
  #
  # setwd("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/3Enrichment/computedSemSim")
  #
  # write.table(edges_sem,
  #   file = paste0(cell, "NetworkEdges.txt"),
  #   sep = "\t",
  #   row.names = F,
  #   quote = F
  # )
}

setwd("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/3Enrichment")
