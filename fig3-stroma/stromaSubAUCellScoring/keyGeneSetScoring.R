# =========== Prep ===========
obj_path <- "/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/*originaldata/Harmony/harmony_All_Stroma_sub.rds"
gs_dir   <- "/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/fig3-stroma/stromaSubAUCellScoring/geneSets"
out_dir  <- "/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/fig3-stroma/stromaSubGSVAScoring"

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(GSVA)
  library(GSEABase)
  library(ggsci)
})

# =========== Load Seurat object ===========
obj <- readRDS(obj_path)
subcol <- "subcluster"
obj$combo <- paste0(obj$orig.ident, "_", obj@meta.data[[subcol]])

# =========== Expression matrix for scoring ===========
DefaultAssay(obj) <- "RNA"
expr <- as.matrix(GetAssayData(obj, slot = "data"))

dup <- duplicated(rownames(expr))
if (any(dup)) expr <- expr[!dup, , drop = FALSE]
rownames(expr) <- toupper(rownames(expr))

# =========== Read gene sets (GMT) ===========
gmt_files <- list.files(gs_dir, pattern = "\\.gmt$", full.names = TRUE)

read_gmt_as_list <- function(gmt_path) {
  gsc <- getGmt(gmt_path)
  lst <- lapply(gsc, function(gs) toupper(geneIds(gs)))
  names(lst) <- paste0(basename(gmt_path), "::", sapply(gsc, setName))
  lst
}
gs_list <- unlist(lapply(gmt_files, read_gmt_as_list), recursive = FALSE)

min_sz <- 10
max_sz <- 500
gs_list <- lapply(gs_list, function(v) intersect(v, rownames(expr)))
gs_list <- gs_list[sapply(gs_list, function(v) length(v) >= min_sz & length(v) <= max_sz)]
if (length(gs_list) == 0) stop("Stop")

# =========== Run GSVA ===========
set.seed(1234)
gp <- GSVA::gsvaParam(
  exprData = expr,
  geneSets = gs_list,
  kcdf     = "Gaussian",
  minSize  = min_sz,
  maxSize  = max_sz,
  maxDiff  = TRUE
)

gsva_scores <- GSVA::gsva(gp)

saveRDS(gsva_scores, file = file.path(out_dir, "GSVA_scores_by_cell.rds"))

# =========== Long format + stats ===========
scores_df <- as.data.frame(t(gsva_scores))
scores_df$cell  <- rownames(scores_df)
scores_df$combo <- obj$combo[rownames(scores_df)]

long_df <- scores_df |>
  tidyr::pivot_longer(cols = -c(cell, combo),
                      names_to = "gene_set",
                      values_to = "score") |>
  dplyr::filter(!is.na(combo))

long_df_norm <- long_df |>
  group_by(gene_set) |>
  mutate(
    mu    = mean(score, na.rm = TRUE),
    sigma = sd(score,  na.rm = TRUE),
    score_z = ifelse(is.finite(sigma) & sigma > 0, (score - mu) / sigma, 0),
    smin  = min(score, na.rm = TRUE),
    smax  = max(score, na.rm = TRUE),
    score_minmax = ifelse(smax > smin, (score - smin) / (smax - smin), 0)
  ) |>
  ungroup() |>
  dplyr::select(-mu, -sigma, -smin, -smax)

summary_df <- long_df_norm |>
  group_by(gene_set, combo) |>
  summarise(
    n = dplyr::n(),
    mean_raw   = mean(score,        na.rm = TRUE),
    median_raw = median(score,      na.rm = TRUE),
    mean_z     = mean(score_z,      na.rm = TRUE),
    median_z   = median(score_z,    na.rm = TRUE),
    mean_mm    = mean(score_minmax, na.rm = TRUE),
    median_mm  = median(score_minmax, na.rm = TRUE),
    .groups = "drop"
  )

p_z <- ggplot(long_df_norm, aes(x = combo, y = score_z, fill = combo)) +
  geom_boxplot(outlier.size = 0.4, width = 0.7) +
  facet_wrap(~ gene_set, scales = "free_y", ncol = 2) +
  scale_fill_npg(guide = "none") +
  labs(title = "GSVA (z-score) by combo", x = "Combo", y = "GSVA z-score") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 9))

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr); library(tidyr); library(ggplot2)
  library(GSEABase)
  library(ggsci)
  library(AUCell)
})

# =========== Run AUCell ===========
set.seed(1234)

rankings <- AUCell_buildRankings(
  expr,
  plotStats = FALSE,
  verbose   = TRUE
)

min_sz <- 5
max_sz <- 500
gs_list <- lapply(gs_list, function(v) intersect(v, rownames(rankings)))
gs_list <- gs_list[sapply(gs_list, function(v) length(v) >= min_sz & length(v) <= max_sz)]

nGenes <- nrow(rankings)
aucMaxRank <- min(nGenes, max(50, ceiling(0.05 * nGenes)))

cellsAUC <- AUCell_calcAUC(gs_list, rankings, aucMaxRank = aucMaxRank, nCores = 1)
auc_mat  <- as.matrix(getAUC(cellsAUC))

scores_df <- as.data.frame(t(auc_mat))
scores_df$cell  <- rownames(scores_df)
scores_df$combo <- obj$combo[rownames(scores_df)]

long_df <- scores_df |>
  tidyr::pivot_longer(cols = -c(cell, combo),
                      names_to = "gene_set",
                      values_to = "score") |>
  dplyr::filter(!is.na(combo))

long_df_norm <- long_df |>
  group_by(gene_set) |>
  mutate(
    mu    = mean(score, na.rm = TRUE),
    sigma = sd(score,  na.rm = TRUE),
    score_z = ifelse(is.finite(sigma) & sigma > 0, (score - mu) / sigma, 0),
    smin  = min(score, na.rm = TRUE),
    smax  = max(score, na.rm = TRUE),
    score_minmax = ifelse(smax > smin, (score - smin) / (smax - smin), 0)
  ) |>
  ungroup() |>
  dplyr::select(-mu, -sigma, -smin, -smax)

summary_df <- long_df_norm |>
  group_by(gene_set, combo) |>
  summarise(
    n = dplyr::n(),
    mean_raw   = mean(score,        na.rm = TRUE),
    median_raw = median(score,      na.rm = TRUE),
    mean_z     = mean(score_z,      na.rm = TRUE),
    median_z   = median(score_z,    na.rm = TRUE),
    mean_mm    = mean(score_minmax, na.rm = TRUE),
    median_mm  = median(score_minmax, na.rm = TRUE),
    .groups = "drop"
  )

p_z <- 
  ggplot(
    long_df_norm, 
         aes(x = combo, y = score_z, fill = combo)
         ) +
  geom_boxplot(outlier.size = 0.35, 
               width = 0.65,
               staplewidth = 1,
               color = "black",
               size = 0.3) +
  # stat_summary(
  #   fun = median,
  #   geom = "crossbar",
  #   width = 0.65,
  #   fatten = 0,
  #   color = "white",
  #   size = 0.1
  # ) +
  facet_wrap(~ gene_set, scales = "free_y", ncol = 2) +
  scale_fill_npg(guide = "none") +
  labs(title = "AUCell AUC (z-score) by combo", x = "Combo", y = "AUC z-score") +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
        axis.text.y = element_text(colour = "black"),
        strip.text = element_text(size = 9),
        axis.line.y = element_line(color = "black", size = 0.3),
        axis.line.x = element_line(color = "black", size = 0.3),
        axis.ticks = element_line(color = "black", size = 0.3))
print(p_z)


make_gs_labels <- function(x) {
  x1 <- sub("^.*::", "", x)

  x2 <- gsub("_", " ", x1)

  x3 <- sub("^\\S+\\s*", "", x2)
  x3 <- trimws(x3)

  x3[x3 == ""] <- x2[x3 == ""]
  x3
}

gs_levels <- unique(long_df_norm$gene_set)
gs_labmap <- setNames(make_gs_labels(gs_levels), gs_levels)

p_z <- 
  ggplot(long_df_norm, aes(x = combo, y = score_z, fill = combo)) +
  geom_boxplot(
    outlier.size = 0.2, width = 0.65, staplewidth = 1,
    color = "black", size = 0.25
  ) +
  facet_wrap(
    ~ gene_set, scales = "free_y", ncol = 2,
    labeller = labeller(gene_set = gs_labmap)
  ) +
  scale_fill_npg(guide = "none") +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
    axis.text.y = element_text(colour = "black"),
    strip.text  = element_text(size = 5),
    axis.line.y = element_line(color = "black", size = 0.3),
    axis.line.x = element_line(color = "black", size = 0.3),
    axis.ticks  = element_line(color = "black", size = 0.3)
  )

print(p_z) # pdf: 2.74 x 7.33

#############################
########### col 2 ###########
#############################

obj$orig.ident <- factor(obj$orig.ident, levels = c("ND", "HFD"))

# Manual colors for ND / HFD
diet_cols <- c("ND" = "#74c5be", "HFD" = "#e95503")

# Link per-cell meta (orig.ident, subcluster) to long table if not already done
if (!all(c("orig.ident", "subcluster") %in% colnames(long_df_norm))) {
  stopifnot(subcol %in% colnames(obj@meta.data))
  cell_meta <- obj@meta.data[, c("orig.ident", subcol)]
  cell_meta$cell <- rownames(cell_meta)
  colnames(cell_meta)[2] <- "subcluster"
  long_df_norm <- dplyr::left_join(long_df_norm, cell_meta, by = "cell") |>
    dplyr::filter(!is.na(orig.ident), !is.na(subcluster))
}

# --------- spacing knobs ---------
dodge_w <- 0.70    # ① within-pair spacing: bigger = ND and HFD farther apart
box_w   <- 0.56    # ② box width: smaller = slimmer boxes (visually more gap)

# Make significance bracket span follow your dodge width automatically
sub_levels <- levels(obj$subcluster)
pos_map <- setNames(seq_along(sub_levels), sub_levels)
span <- 0.5 * dodge_w * 0.95   # cover both ND & HFD; tweak 0.95 if you want wider/narrower brackets

anno_df <- anno_df |>
  dplyr::mutate(
    xmin = pos_map[subcluster] - span,
    xmax = pos_map[subcluster] + span
  )

# --------- plot ---------
p_pairs <- ggplot(
  long_df_norm,
  aes(x = subcluster, y = score_z,
      fill = orig.ident,                               # color by condition
      group = interaction(subcluster, orig.ident))     # keep two boxes per subcluster
) +
  geom_boxplot(
    position = position_dodge(width = dodge_w),        # ① within-pair spacing
    width = box_w,                                     # ② box width
    outlier.size = 0.30,
    color = "black",
    size = 0.28
  ) +
  # significance brackets across ND-HFD per subcluster (manual mode)
  ggsignif::geom_signif(
    data = anno_df,
    aes(xmin = xmin, xmax = xmax,
        annotations = p_label, y_position = y_position),
    inherit.aes = FALSE,
    manual = TRUE,
    tip_length = 0.01,
    textsize = 2.8,
    vjust = 0
  ) +
  coord_cartesian(clip = "off") +
  facet_wrap(
    ~ gene_set, scales = "free_y", ncol = 2,
    labeller = labeller(gene_set = gs_labmap)
  ) +
  scale_fill_manual(values = diet_cols, breaks = c("ND","HFD"), name = "Diet") + 
  labs(
    x = "Cell type (subcluster)",
    y = "AUC z-score"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
    axis.text.y = element_text(colour = "black"),
    strip.text  = element_text(size = 7),
    axis.line   = element_line(color = "black", size = 0.3),
    axis.ticks  = element_line(color = "black", size = 0.3),
    panel.spacing = unit(0.9, "lines"),
    legend.position = "right"
  )

print(p_pairs)

# ===================== Build and save stats table (per gene_set × subcluster) =====================
# Requirements: long_df_norm has columns: cell, gene_set, score, score_z, orig.ident (ND/HFD), subcluster
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

# Safety checks
stopifnot(all(c("gene_set", "subcluster", "orig.ident", "score", "score_z") %in% colnames(long_df_norm)))
long_df_norm$orig.ident <- factor(long_df_norm$orig.ident, levels = c("ND", "HFD"))

# 1) Condition-level summaries (ND vs HFD) for each gene_set × subcluster
cond_summ <- long_df_norm %>%
  group_by(gene_set, subcluster, orig.ident) %>%
  summarise(
    n          = dplyr::n(),
    mean_z     = mean(score_z, na.rm = TRUE),
    median_z   = median(score_z, na.rm = TRUE),
    mean_raw   = mean(score,    na.rm = TRUE),
    median_raw = median(score,  na.rm = TRUE),
    .groups = "drop"
  )

# Wide format: columns for ND and HFD
wide <- cond_summ %>%
  tidyr::pivot_wider(
    names_from  = orig.ident,
    values_from = c(n, mean_z, median_z, mean_raw, median_raw),
    names_sep   = "."
  )

# 2) Wilcoxon rank-sum test on z-scores + Hodges–Lehmann (HL) location shift
#    Note: wilcox.test estimate = ND - HFD, so we flip sign to report HL for (HFD - ND)
test_df <- long_df_norm %>%
  group_by(gene_set, subcluster) %>%
  dplyr::group_modify(~{
    x <- .x$score_z[.x$orig.ident == "ND"]
    y <- .x$score_z[.x$orig.ident == "HFD"]
    if (length(x) > 0 && length(y) > 0) {
      wt <- suppressWarnings(wilcox.test(x, y, alternative = "two.sided", conf.int = TRUE))
      tibble::tibble(
        n_ND  = length(x),
        n_HFD = length(y),
        p_value = wt$p.value,
        # HL for ND - HFD (from wilcox.test), then convert to HFD - ND:
        HL_shift_HFD_minus_ND  = -unname(wt$estimate),
        conf_low_HFD_minus_ND  = -wt$conf.int[2],  # flip interval endpoints
        conf_high_HFD_minus_ND = -wt$conf.int[1]
      )
    } else {
      tibble::tibble(
        n_ND = length(x), n_HFD = length(y),
        p_value = NA_real_,
        HL_shift_HFD_minus_ND = NA_real_,
        conf_low_HFD_minus_ND = NA_real_,
        conf_high_HFD_minus_ND = NA_real_
      )
    }
  }) %>%
  ungroup()

stats_tbl <- wide %>%
  left_join(test_df, by = c("gene_set", "subcluster")) %>%
  mutate(
    delta_mean_z     = `mean_z.HFD`     - `mean_z.ND`,
    delta_median_z   = `median_z.HFD`   - `median_z.ND`,
    delta_mean_raw   = `mean_raw.HFD`   - `mean_raw.ND`,
    delta_median_raw = `median_raw.HFD` - `median_raw.ND`,
    p_adj_BH = p.adjust(p_value, method = "BH"),
    sig = dplyr::case_when(
      is.na(p_adj_BH)     ~ "NA",
      p_adj_BH < 0.001    ~ "***",
      p_adj_BH < 0.01     ~ "**",
      p_adj_BH < 0.05     ~ "*",
      TRUE                ~ "ns"
    )
  ) %>%
  # Reorder and round columns for a clean table
  transmute(
    gene_set, subcluster,
    n_ND = coalesce(`n.ND`, n_ND),
    n_HFD = coalesce(`n.HFD`, n_HFD),

    mean_z_ND   = `mean_z.ND`,
    mean_z_HFD  = `mean_z.HFD`,
    delta_mean_z = delta_mean_z,

    median_z_ND   = `median_z.ND`,
    median_z_HFD  = `median_z.HFD`,
    delta_median_z = delta_median_z,

    mean_raw_ND   = `mean_raw.ND`,
    mean_raw_HFD  = `mean_raw.HFD`,
    delta_mean_raw = delta_mean_raw,

    median_raw_ND   = `median_raw.ND`,
    median_raw_HFD  = `median_raw.HFD`,
    delta_median_raw = delta_median_raw,

    HL_shift_HFD_minus_ND,
    conf_low_HFD_minus_ND,
    conf_high_HFD_minus_ND,

    p_value,
    p_adj_BH,
    sig
  ) %>%
  mutate(across(where(is.numeric), ~round(., 4))) %>%
  arrange(gene_set, subcluster)

write.csv(stats_tbl, file = "/Users/coellearth/Desktop/Mammary_Gland_Diet_Project/fig3-stroma/stromaSubAUCellScoring", row.names = FALSE)

print(head(stats_tbl, 12))
message("Saved stats table to: ", out_csv)

#####################################
#####################################
#####################################
# ======================== GSVA: tidy, visualize, and stats ========================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
  library(ggsci); library(ggsignif)
})

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- Long format ---
scores_df <- as.data.frame(t(gsva_scores))
scores_df$cell  <- rownames(scores_df)
scores_df$combo <- obj$combo[rownames(scores_df)]

long_df <- scores_df |>
  tidyr::pivot_longer(cols = -c(cell, combo),
                      names_to = "gene_set",
                      values_to = "score") |>
  dplyr::filter(!is.na(combo))

# --- Normalize (z and min-max) ---
long_df_norm <- long_df |>
  group_by(gene_set) |>
  mutate(
    mu    = mean(score, na.rm = TRUE),
    sigma = sd(score,  na.rm = TRUE),
    score_z = ifelse(is.finite(sigma) & sigma > 0, (score - mu) / sigma, 0),
    smin  = min(score, na.rm = TRUE),
    smax  = max(score, na.rm = TRUE),
    score_minmax = ifelse(smax > smin, (score - smin) / (smax - smin), 0)
  ) |>
  ungroup() |>
  dplyr::select(-mu, -sigma, -smin, -smax)

# --- Combo-wise boxplots (like your first GSVA plot) ---
p_combo_z <- ggplot(long_df_norm, aes(x = combo, y = score_z, fill = combo)) +
  geom_boxplot(outlier.size = 0.35, width = 0.7, color = "black", size = 0.25) +
  facet_wrap(~ gene_set, scales = "free_y", ncol = 2) +
  scale_fill_npg(guide = "none") +
  labs(title = "GSVA (z-score) by combo", x = "Combo", y = "GSVA z-score") +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
        axis.text.y = element_text(colour = "black"),
        strip.text  = element_text(size = 8),
        axis.line   = element_line(color = "black", size = 0.3),
        axis.ticks  = element_line(color = "black", size = 0.3))
print(p_combo_z)

# Save (optional)
ggsave(file.path(out_dir, "GSVA_combo_boxplot_z.pdf"), p_combo_z, width = 7, height = 10, useDingbats = FALSE)

# --- Pretty labels like your AUCell helper ---
make_gs_labels <- function(x) {
  x1 <- sub("^.*::", "", x)
  x2 <- gsub("_", " ", x1)
  x3 <- sub("^\\S+\\s*", "", x2); x3 <- trimws(x3)
  x3[x3 == ""] <- x2[x3 == ""]
  x3
}
gs_levels <- unique(long_df_norm$gene_set)
gs_labmap <- setNames(make_gs_labels(gs_levels), gs_levels)

# ================= ND vs HFD per subcluster figure (with brackets) =================
# Ensure meta cols are present
if (!all(c("orig.ident", "subcluster") %in% colnames(long_df_norm))) {
  stopifnot(subcol %in% colnames(obj@meta.data))
  cell_meta <- obj@meta.data[, c("orig.ident", subcol)]
  cell_meta$cell <- rownames(cell_meta)
  colnames(cell_meta)[2] <- "subcluster"
  long_df_norm <- dplyr::left_join(long_df_norm, cell_meta, by = "cell") |>
    dplyr::filter(!is.na(orig.ident), !is.na(subcluster))
}

# Order and colors
obj$orig.ident <- factor(obj$orig.ident, levels = c("ND", "HFD"))
long_df_norm$orig.ident <- factor(long_df_norm$orig.ident, levels = c("ND", "HFD"))
diet_cols <- c("ND" = "#74c5be", "HFD" = "#e95503")

# Wilcoxon (z-scores) per gene_set × subcluster
test_df <- long_df_norm %>%
  group_by(gene_set, subcluster) %>%
  dplyr::group_modify(~{
    x <- .x$score_z[.x$orig.ident == "ND"]
    y <- .x$score_z[.x$orig.ident == "HFD"]
    if (length(x) > 0 && length(y) > 0) {
      wt <- suppressWarnings(wilcox.test(x, y, alternative = "two.sided", conf.int = TRUE))
      tibble::tibble(
        p_value = wt$p.value,
        HL_shift_HFD_minus_ND  = -unname(wt$estimate),
        conf_low_HFD_minus_ND  = -wt$conf.int[2],
        conf_high_HFD_minus_ND = -wt$conf.int[1]
      )
    } else {
      tibble::tibble(
        p_value = NA_real_,
        HL_shift_HFD_minus_ND  = NA_real_,
        conf_low_HFD_minus_ND  = NA_real_,
        conf_high_HFD_minus_ND = NA_real_
      )
    }
  }) %>% ungroup() %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH"),
         p_label = dplyr::case_when(
           is.na(p_adj_BH)  ~ "NA",
           p_adj_BH < 0.001 ~ "***",
           p_adj_BH < 0.01  ~ "**",
           p_adj_BH < 0.05  ~ "*",
           TRUE             ~ "ns"
         ))

# Bracket y-positions per panel: slightly above the local max
panel_max <- long_df_norm %>%
  group_by(gene_set, subcluster) %>%
  summarise(ymax = max(score_z, na.rm = TRUE), .groups = "drop")
test_df <- test_df %>%
  left_join(panel_max, by = c("gene_set", "subcluster")) %>%
  mutate(y_position = ifelse(is.finite(ymax), ymax + 0.2, 0.2))

# Compute manual x-span per subcluster (paired ND/HFD)
dodge_w <- 0.70
box_w   <- 0.56
sub_levels <- if (!is.null(levels(obj$subcluster))) levels(obj$subcluster) else sort(unique(as.character(obj$subcluster)))
long_df_norm$subcluster <- factor(long_df_norm$subcluster, levels = sub_levels)
pos_map <- setNames(seq_along(sub_levels), sub_levels)
span <- 0.5 * dodge_w * 0.95

anno_df <- test_df %>%
  mutate(
    xmin = pos_map[subcluster] - span,
    xmax = pos_map[subcluster] + span
  )

p_pairs <- ggplot(
  long_df_norm,
  aes(x = subcluster, y = score_z,
      fill = orig.ident,
      group = interaction(subcluster, orig.ident))
) +
  geom_boxplot(
    position = position_dodge(width = dodge_w),
    width = box_w,
    outlier.size = 0.30,
    color = "black",
    size = 0.28
  ) +
  ggsignif::geom_signif(
    data = anno_df,
    aes(xmin = xmin, xmax = xmax,
        annotations = p_label, y_position = y_position),
    inherit.aes = FALSE,
    manual = TRUE,
    tip_length = 0.01,
    textsize = 2.8,
    vjust = 0
  ) +
  coord_cartesian(clip = "off") +
  facet_wrap(
    ~ gene_set, scales = "free_y", ncol = 2,
    labeller = labeller(gene_set = gs_labmap)
  ) +
  scale_fill_manual(values = diet_cols, breaks = c("ND","HFD"), name = "Diet") +
  labs(
    title = "GSVA (z-score) by subcluster and diet",
    x = "Cell type (subcluster)",
    y = "GSVA z-score"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
    axis.text.y = element_text(colour = "black"),
    strip.text  = element_text(size = 7),
    axis.line   = element_line(color = "black", size = 0.3),
    axis.ticks  = element_line(color = "black", size = 0.3),
    panel.spacing = unit(0.9, "lines"),
    legend.position = "right"
  )
print(p_pairs)

ggsave(file.path(out_dir, "GSVA_by_subcluster_diet_boxplot_z.pdf"), p_pairs, width = 7.2, height = 10.5, useDingbats = FALSE)

# ===================== Significance stats table (GSVA) =====================
# Condition-level summaries and deltas (mirrors your AUCell table)
cond_summ <- long_df_norm %>%
  group_by(gene_set, subcluster, orig.ident) %>%
  summarise(
    n          = dplyr::n(),
    mean_z     = mean(score_z, na.rm = TRUE),
    median_z   = median(score_z, na.rm = TRUE),
    mean_raw   = mean(score,    na.rm = TRUE),
    median_raw = median(score,  na.rm = TRUE),
    .groups = "drop"
  )

wide <- cond_summ %>%
  tidyr::pivot_wider(
    names_from  = orig.ident,
    values_from = c(n, mean_z, median_z, mean_raw, median_raw),
    names_sep   = "."
  )

stats_tbl <- wide %>%
  left_join(test_df %>% dplyr::select(gene_set, subcluster,
                               HL_shift_HFD_minus_ND,
                               conf_low_HFD_minus_ND,
                               conf_high_HFD_minus_ND,
                               p_value, p_adj_BH),
            by = c("gene_set", "subcluster")) %>%
  mutate(
    delta_mean_z     = `mean_z.HFD`     - `mean_z.ND`,
    delta_median_z   = `median_z.HFD`   - `median_z.ND`,
    delta_mean_raw   = `mean_raw.HFD`   - `mean_raw.ND`,
    delta_median_raw = `median_raw.HFD` - `median_raw.ND`,
    sig = dplyr::case_when(
      is.na(p_adj_BH)     ~ "NA",
      p_adj_BH < 0.001    ~ "***",
      p_adj_BH < 0.01     ~ "**",
      p_adj_BH < 0.05     ~ "*",
      TRUE                ~ "ns"
    )
  ) %>%
  transmute(
    gene_set, subcluster,
    n_ND = `n.ND`, n_HFD = `n.HFD`,

    mean_z_ND = `mean_z.ND`,   mean_z_HFD = `mean_z.HFD`,   delta_mean_z,
    median_z_ND = `median_z.ND`, median_z_HFD = `median_z.HFD`, delta_median_z,

    mean_raw_ND = `mean_raw.ND`, mean_raw_HFD = `mean_raw.HFD`, delta_mean_raw,
    median_raw_ND = `median_raw.ND`, median_raw_HFD = `median_raw.HFD`, delta_median_raw,

    HL_shift_HFD_minus_ND,
    conf_low_HFD_minus_ND,
    conf_high_HFD_minus_ND,

    p_value, p_adj_BH, sig
  ) %>%
  mutate(across(where(is.numeric), ~round(., 4))) %>%
  arrange(gene_set, subcluster)

# Write out (fixed file name)
out_csv <- file.path(out_dir, "GSVA_stats_per_geneSet_subcluster.csv")
write.csv(stats_tbl, file = out_csv, row.names = FALSE)