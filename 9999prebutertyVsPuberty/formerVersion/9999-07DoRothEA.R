library(Seurat)
library(dorothea)
library(decoupleR)
library(dplyr)
library(tidyr)
library(purrr)
library(pheatmap)

# cds_obj <- readRDS("/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/9999prebutertyVsPuberty/formerVersion/trajectoryOutput/cds_obj.rds")

All_stage <- readRDS(file = "/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/9999prebutertyVsPuberty/formerVersion/trajectoryOutput/All_stages_trajectory.rds")

DefaultAssay(All_stage) <- "RNA"
expr <- as.matrix(GetAssayData(All_stage, slot = "data"))

## ===== construct DoRothEA network =====
data(dorothea_mm)
regulons_abc <- dorothea_mm %>%
  filter(confidence %in% c("A","B","C")) %>%
  transmute(tf = tf, 
            target = target, 
            mor = mor, 
            confidence = confidence)

common_genes <- intersect(rownames(expr), unique(regulons_abc$target))
expr <- expr[common_genes, , drop = FALSE]
regulons_abc <- regulons_abc %>% filter(target %in% common_genes)

sig_tbl <- regulons_abc %>%
  arrange(tf, target) %>%
  group_by(tf) %>%
  summarise(sig = paste(target, mor, collapse="|"), .groups = "drop")

dup_groups <- sig_tbl %>% group_by(sig) %>% filter(n() > 1)

if (nrow(dup_groups) > 0) {
  conf_w <- c(A=3,B=2,C=1)
  keep_dup <- dup_groups %>%
    left_join(regulons_abc, by="tf") %>%
    group_by(sig, tf) %>%
    summarise(A_targets = sum(confidence=="A"),
              total_targets = n(),
              score = sum(conf_w[confidence]), .groups="drop") %>%
    group_by(sig) %>%
    slice_max(order_by = A_targets, n=1, with_ties = FALSE) %>%
    pull(tf)
  drop_dup <- setdiff(dup_groups$tf, keep_dup)
  regulons_abc <- regulons_abc %>% filter(!tf %in% drop_dup)
}

cc <- decoupleR::check_corr(
  regulons_abc, 
  .source="tf", 
  .target="target", 
  .mor="mor"
)

bad <- cc %>% filter(abs(correlation) >= 0.9)

if (nrow(bad) > 0) {
  g <- graph_from_data_frame(bad %>% transmute(from=source, to=source.2), directed=FALSE)
  comps <- components(g)$membership
  groups <- split(names(comps), comps)
  conf_w <- c(A=3,B=2,C=1)
  keep <- map_chr(groups, function(gs){
    regulons_abc %>% filter(tf %in% gs) %>%
      group_by(tf) %>%
      summarise(A_targets = sum(confidence=="A"),
                total_targets = n(),
                score = sum(conf_w[confidence]), .groups="drop") %>%
      arrange(desc(A_targets), desc(score), desc(total_targets)) %>%
      slice(1) %>% pull(tf)
  })
  drop <- setdiff(unique(c(bad$source, bad$source.2)), keep)
  message("Dropping ", length(drop), " highly collinear TFs: ", paste(head(drop, 10), collapse=", "),
          ifelse(length(drop)>10, " ...", ""))
  regulons_abc <- regulons_abc %>% filter(!tf %in% drop)
}


## ===== decoupleR::run_mlm to estimate TF activity =====
mlm_res <- run_mlm(
  mat       = expr,
  network   = regulons_abc,
  .source   = "tf",
  .target   = "target",
  .mor      = "mor",
  minsize   = 5
)

mlm_mat <- mlm_res %>%
  select(source, condition, score) %>%
  rename(tf = source, cell = condition, score = score) %>%
  pivot_wider(names_from = cell, values_from = score) %>%
  tibble::column_to_rownames("tf") %>%
  as.matrix()

All_stage[["dorothea_mlm"]] <- CreateAssayObject(data = mlm_mat)

# Split into time boxes according to pseudotime

All_stage <- FindClusters(All_stage, resolution = 0.2)

DimPlot(All_stage, 
        reduction = "umap", 
        group.by = "seurat_clusters", 
        label = TRUE, 
        repel = TRUE)

K <- 15
All_stage$pt_bin <- ggplot2::cut_number(All_stage$monocle3_pseudotime, 
                               n = K, 
                               labels = paste0("B", 1:K))

DimPlot(All_stage, 
        reduction = "umap", 
        group.by = "pt_bin", 
        label = TRUE, 
        repel = TRUE)

DefaultAssay(All_stage) <- "dorothea_mlm"

bins <- paste0("B", 1:K)
transitions_B <- tibble(from = bins[-length(bins)], to = bins[-1])

DE_list_B <- pmap(transitions_B, function(from, to){
  
  keep <- which(All_stage$pt_bin %in% c(from, to))
  
  sub  <- subset(All_stage, cells = colnames(All_stage)[keep])
  
  Idents(sub) <- 
    factor(
      ifelse(sub$pt_bin == to, "to", "from"), 
      levels = c("from","to"))
  
  res <- FindMarkers(sub, 
                     ident.1 = "to", 
                     ident.2 = "from",
                     logfc.threshold = 0, 
                     min.pct = 0,
                     test.use = "t")
  
  res$tf <- rownames(res)
  
  res$from <- from
  res$to <- to
  
  res <- res[order(res$avg_log2FC, decreasing = TRUE), ]
  up   <- head(res, 5)
  down <- head(res[order(res$avg_log2FC, decreasing = FALSE), ], 5)
  list(all = res, up = up, down = down)
})

top_tbl_B <- map2_df(DE_list_B, seq_along(DE_list_B), function(x, i){
  bind_rows(
    mutate(x$up,   direction = "up"),
    mutate(x$down, direction = "down")
  ) %>% mutate(step = i) %>% select(step, from, to, tf, avg_log2FC, p_val_adj, direction)
})

####################################
####################################
####################################

pga <- monocle3::principal_graph_aux(cds_obj)
closest_idx <- as.integer(pga$UMAP$pr_graph_cell_proj_closest_vertex[, 1])
node_names  <- colnames(pga$UMAP$dp_mst)
closest_node <- node_names[closest_idx]
names(closest_node) <- rownames(pga$UMAP$pr_graph_cell_proj_closest_vertex)

pt <- monocle3::pseudotime(cds_obj)
All_stage$monocle3_pseudotime <- pt[Seurat::Cells(All_stage)]
All_stage$node <- closest_node[Seurat::Cells(All_stage)]

node_order <- as.data.frame(Seurat::FetchData(All_stage, vars = c("node","monocle3_pseudotime"))) %>%
  group_by(node) %>%
  summarise(med_pt = median(monocle3_pseudotime, na.rm = TRUE), .groups = "drop") %>%
  arrange(med_pt) %>% pull(node)

DimPlot(All_stage, 
        reduction = "umap", 
        group.by = "node", 
        label = F, 
        repel = F)
