my_mushroom_cellchat <- function(cellchat_table,
                                 source_filter = NULL,
                                 target_filter = NULL,
                                 top_n = NULL,          # NULL = 显示全部
                                 show_rankings = FALSE,
                                 ND_color = "#71BDB6",
                                 HFD_color = "#DD5511",
                                 circle_scale = 0.7,
                                 rowname_size = 8,
                                 colname_size = 10,
                                 legend_title_size = 10,
                                 legend_text_size = 8,
                                 ...)  {
  
  # 必要包（只检查，不 attach）
  requireNamespace("dplyr")
  requireNamespace("tidyr")
  requireNamespace("ggplot2")
  requireNamespace("ggforce")
  requireNamespace("shadowtext")
  
  # 1. 检查必需列 --------------------------------------------------------------
  required_cols <- c("source", "target", "LR_pair", "prob_ND", "prob_HFD")
  if (!all(required_cols %in% colnames(cellchat_table))) {
    stop(paste(
      "Missing required columns:",
      paste(required_cols[!required_cols %in% colnames(cellchat_table)],
            collapse = ", ")
    ))
  }
  
  # 2. 筛选数据 ---------------------------------------------------------------
  filtered_data <- cellchat_table
  
  if (!is.null(source_filter)) {
    filtered_data <- dplyr::filter(filtered_data, source %in% source_filter)
  }
  if (!is.null(target_filter)) {
    filtered_data <- dplyr::filter(filtered_data, target %in% target_filter)
  }
  
  if (nrow(filtered_data) == 0) {
    stop("No data after filtering")
  }
  
  # 3. 计算总概率，用于排序 ----------------------------------------------------
  filtered_data <- dplyr::mutate(filtered_data,
                                 total_prob = prob_ND + prob_HFD)
  
  # 4. 设置因子顺序 -----------------------------------------------------------
  if (!is.null(source_filter)) {
    filtered_data$source <- factor(filtered_data$source,
                                   levels = source_filter)
  } else {
    filtered_data$source <- factor(filtered_data$source,
                                   levels = unique(filtered_data$source))
  }
  
  if (!is.null(target_filter)) {
    filtered_data$target <- factor(filtered_data$target,
                                   levels = target_filter)
  } else {
    filtered_data$target <- factor(filtered_data$target,
                                   levels = unique(filtered_data$target))
  }
  
  # 5. 排序 + top_n 选择 ------------------------------------------------------
  filtered_data <- dplyr::arrange(filtered_data,
                                  source, target, dplyr::desc(total_prob))
  
  if (!is.null(top_n) && top_n < nrow(filtered_data)) {
    # 保留每个 (source, target) 的一定数量，让组合比较均匀
    n_combinations <- length(unique(paste(filtered_data$source,
                                          filtered_data$target)))
    per_group <- ceiling(top_n / n_combinations)
    
    filtered_data <- filtered_data |>
      dplyr::group_by(source, target) |>
      dplyr::slice_head(n = per_group) |>
      dplyr::ungroup() |>
      dplyr::slice_head(n = top_n)
  }
  
  if (nrow(filtered_data) == 0) {
    stop("No data to plot after applying top_n / filters")
  }
  
  # 6. 全局排名 & 行标签 ------------------------------------------------------
  filtered_data <- dplyr::mutate(filtered_data,
                                 prioritization_rank = dplyr::row_number())
  
  filtered_data <- dplyr::mutate(
    filtered_data,
    row_label = paste0(source, " ", target, " ", LR_pair)
  )
  
  # 当前要展示的行顺序
  order_rows <- filtered_data$row_label
  n_rows     <- length(order_rows)
  
  # 每个 row_label 对应的 y 坐标（第一个在最上面）
  #   order_rows[1] -> y = n_rows
  #   order_rows[n_rows] -> y = 1
  row_pos <- setNames(n_rows:1, order_rows)
  
  # 7. 宽变长：一条互作 → ND/HFD 两行 --------------------------------------
  long_data <- filtered_data |>
    dplyr::select(source, target, row_label,
                  prob_ND, prob_HFD,
                  prioritization_rank, total_prob) |>
    tidyr::pivot_longer(
      cols = c(prob_ND, prob_HFD),
      names_to = "condition",
      values_to = "probability",
      names_prefix = "prob_"
    )
  
  long_data$row_label <- factor(long_data$row_label,
                                levels = order_rows)
  long_data$source <- factor(long_data$source,
                             levels = levels(filtered_data$source))
  
  # 8. 计算坐标和半径 ---------------------------------------------------------
  max_prob <- max(c(filtered_data$prob_ND,
                    filtered_data$prob_HFD), na.rm = TRUE)
  
  long_data <- dplyr::mutate(
    long_data,
    x = as.numeric(source),                          # X: sender index
    y = row_pos[as.character(row_label)],           # Y: 行 index（已翻转）
    start = ifelse(condition == "ND", -pi, 0),      # ND 左半圆, HFD 右半圆
    size = if (max_prob > 0) {
      ifelse(probability > 0,
             sqrt(probability / max_prob), 0)
    } else {
      0
    }
  )
  
  # 9. 图形辅助数据 -----------------------------------------------------------
  scale     <- circle_scale
  n_sources <- length(levels(long_data$source))
  
  legend_df <- data.frame(
    values = c(0.25, 0.5, 0.75, 1),
    x = (n_sources + 2.5):(n_sources + 5.5),
    y = rep(floor(n_rows / 3), 4),
    start = -pi
  )
  
  axis_rect <- data.frame(
    xmin = 0, xmax = n_sources + 1,
    ymin = 0, ymax = n_rows + 1
  )
  
  panel_grid_y <- data.frame(
    x = rep(seq(from = 0.5, to = n_sources + 0.5, by = 1), each = 2),
    y = c(n_rows + 1, 0),
    group = rep(1:(n_sources + 1), each = 2)
  )
  
  panel_grid_x <- data.frame(
    y = rep(seq(from = 0.5, to = n_rows + 0.5, by = 1), each = 2),
    x = c(n_sources + 1, 0),
    group = rep(1:(n_rows + 1), each = 2)
  )
  
  # 10. 主题参数 --------------------------------------------------------------
  theme_args <- list(
    panel.grid.major = ggplot2::element_blank(),
    legend.box       = "horizontal",
    panel.background = ggplot2::element_blank(),
    axis.text.y      = ggplot2::element_text(size = rowname_size),
    axis.text.x      = ggplot2::element_text(size = colname_size),
    legend.title     = ggplot2::element_text(size = legend_title_size),
    legend.text      = ggplot2::element_text(size = legend_text_size)
  )
  theme_args[names(list(...))] <- list(...)
  
  if (!"legend.justification" %in% names(theme_args)) {
    theme_args$legend.justification <- c(1, 0.7)
  }
  if (!"legend.position" %in% names(theme_args)) {
    theme_args$legend.position <- c(1, 0.7)
  }
  
  # 11. 绘图主体 --------------------------------------------------------------
  p1 <- ggplot2::ggplot() +
    # ND 左半圆
    ggforce::geom_arc_bar(
      data = dplyr::filter(long_data, condition == "ND", probability > 0),
      ggplot2::aes(x0 = x, y0 = y,
                   r0 = 0, r = size * scale,
                   start = start, end = start + pi),
      fill  = ND_color,
      color = "white",
      linewidth = 0.5
    ) +
    # HFD 右半圆
    ggforce::geom_arc_bar(
      data = dplyr::filter(long_data, condition == "HFD", probability > 0),
      ggplot2::aes(x0 = x, y0 = y,
                   r0 = 0, r = size * scale,
                   start = start, end = start + pi),
      fill  = HFD_color,
      color = "white",
      linewidth = 0.5
    ) +
    # 大小图例
    ggforce::geom_arc_bar(
      data = legend_df,
      ggplot2::aes(x0 = x, y0 = y,
                   r0 = 0, r = sqrt(values) * scale,
                   start = start, end = start + pi),
      fill = "black"
    ) +
    ggplot2::geom_rect(
      data = legend_df,
      ggplot2::aes(xmin = x - 0.5, xmax = x + 0.5,
                   ymin = y - 0.5, ymax = y + 0.5),
      color = "gray90", fill = NA
    ) +
    ggplot2::geom_text(
      data = legend_df,
      ggplot2::aes(label = values, x = x, y = y - 0.6),
      vjust = 1, size = legend_text_size / 2.8
    ) +
    ggplot2::geom_text(
      data = data.frame(
        x = (n_sources + 4),
        y = floor(n_rows / 3) + 1,
        label = "Relative\nProbability"
      ),
      ggplot2::aes(x = x, y = y, label = label),
      size = legend_title_size / 2.8,
      vjust = 0,
      lineheight = 0.75
    ) +
    # 网格线和边框
    ggplot2::geom_line(
      data = panel_grid_y,
      ggplot2::aes(x = x, y = y, group = group),
      color = "gray90"
    ) +
    ggplot2::geom_line(
      data = panel_grid_x,
      ggplot2::aes(x = x, y = y, group = group),
      color = "gray90"
    ) +
    ggplot2::geom_rect(
      data = axis_rect,
      ggplot2::aes(xmin = xmin, ymin = ymin,
                   xmax = xmax, ymax = ymax),
      color = "black", fill = "transparent"
    ) +
    # 坐标轴
    ggplot2::scale_y_continuous(
      breaks = n_rows:1,        # y=n_rows 在最上面
      labels = order_rows,      # 与 row_pos 对齐
      expand = ggplot2::expansion(add = c(0, 0))
    ) +
    ggplot2::scale_x_continuous(
      breaks  = seq_len(n_sources),
      labels  = levels(long_data$source),
      position = "top",
      expand  = ggplot2::expansion(add = c(0, 0))
    ) +
    ggplot2::xlab("Sender cell types") +
    ggplot2::ylab("Sender Target (LR pair)") +
    ggplot2::coord_fixed() +
    do.call(ggplot2::theme, theme_args)
  
  # 12. 排名标签（可选） ------------------------------------------------------
  if (show_rankings) {
    p1 <- p1 +
      shadowtext::geom_shadowtext(
        data = dplyr::filter(long_data, condition == "ND", probability > 0),
        ggplot2::aes(x = x, y = y, label = prioritization_rank),
        size = legend_text_size / 2.8
      )
  }
  
  # 13. 左右条件说明 ----------------------------------------------------------
  p1 <- p1 +
    ggplot2::annotate(
      "text",
      x = n_sources + 4,
      y = n_rows,
      label = "Left: ND\nRight: HFD",
      size = legend_text_size / 2.8,
      hjust = 0.5,
      lineheight = 0.75
    )
  
  return(p1)
}

setwd("/Users/renyixiang/Desktop/data/23BMI/ND_HFD_MG_snRNAseq/Figure7_immune/CellChat")
your_data <- read.csv("LRpairs.csv")

#然后绘图
p <- my_mushroom_cellchat(
  cellchat_table = your_data,
  source_filter = c("MacS01", "MacS03"),
  target_filter = c("Basal", "LumProg", "HormSens"),
  top_n = NULL,
  ND_color = "#71BDB6",
  HFD_color = "#DD5511",
  circle_scale = 0.5,
  rowname_size = 7,
  colname_size = 11
)
print(p)
