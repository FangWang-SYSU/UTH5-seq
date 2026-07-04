suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(ggpubr)
  library(ggplot2)
  library(ggraph)
  library(ggrepel)
  library(igraph)
  library(readxl)
  library(scales)
})

Sys.setlocale("LC_ALL", "en_US.UTF-8")

source_data_file <- "../source_data/Figure_3_Source_Data.xlsx"
if (!file.exists(source_data_file)) {
  stop("Cannot find source data file: ", source_data_file)
}

output_dir <- "Figure3_output"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

module_colors <- c(
  Module1 = "#003C30", Module2 = "#43AC5E", Module3 = "#2D4999",
  Module4 = "#C9B2D6", Module5 = "#E6780C", Module6 = "#7E9A29",
  Module7 = "#C51E7E", Module8 = "#723F91"
)

tf_module_colors <- c(
  TFc1 = "#4EBCD5", TFc2 = "#E74C35", TFc3 = "#07A087",
  TFc4 = "#A75728", TFc5 = "#F193A8", TFc6 = "#BDBD21",
  TFc7 = "#FABE70", TFc8 = "#3D5589", TFc9 = "#B19D86"
)

read_source_data <- function(sheet_name, ...) {
  read_excel(source_data_file, sheet = sheet_name, ...) %>%
    as.data.frame(check.names = FALSE)
}

read_matrix_source <- function(sheet_name) {
  df <- read_source_data(sheet_name)
  row_ids <- df$transcription_factor
  mat <- as.matrix(df[, setdiff(names(df), "transcription_factor"), drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- row_ids
  mat
}

infer_tf_annotation <- function(tf_order, n_modules = 9) {
  module_id <- cut(seq_along(tf_order), breaks = n_modules, labels = FALSE)
  data.frame(
    transcription_factor = tf_order,
    tf_module = paste0("TFc", module_id),
    stringsAsFactors = FALSE
  )
}

plot_tf_heatmap <- function(sheet_name, output_file, legend_title, color_limits, color_values) {
  mat <- read_matrix_source(sheet_name)
  tf_annotation <- infer_tf_annotation(rownames(mat))
  ordered_tfs <- intersect(tf_annotation$transcription_factor, rownames(mat))
  mat <- mat[ordered_tfs, ordered_tfs]
  tf_module <- tf_annotation$tf_module[match(ordered_tfs, tf_annotation$transcription_factor)]

  clipped_mat <- mat
  clipped_mat[clipped_mat < min(color_limits)] <- min(color_limits)
  clipped_mat[clipped_mat > max(color_limits)] <- max(color_limits)

  heatmap_colors <- colorRamp2(color_limits, color_values)
  row_anno <- rowAnnotation(
    TF_module = tf_module,
    col = list(TF_module = tf_module_colors)
  )
  col_anno <- HeatmapAnnotation(
    TF_module = tf_module,
    col = list(TF_module = tf_module_colors)
  )

  pdf(file.path(output_dir, output_file), width = 12, height = 10)
  heatmap <- Heatmap(
    clipped_mat,
    show_row_names = FALSE,
    show_column_names = FALSE,
    show_column_dend = FALSE,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    col = heatmap_colors,
    name = legend_title,
    use_raster = TRUE,
    raster_quality = 4,
    row_split = tf_module,
    column_split = tf_module,
    top_annotation = col_anno,
    right_annotation = row_anno
  )
  draw(
    heatmap,
    heatmap_legend_side = "left",
    annotation_legend_side = "left",
    legend_grouping = "original"
  )
  dev.off()

  heatmap
}

plot_lr_bipartite <- function(df,
                              left = "TF_module",
                              right = "TF_module1",
                              edge_width = "tf_prop",
                              edge_color = "meanDiffcors",
                              node_size = 3.8,
                              left_title = "NTCs",
                              right_title = "Targeted cells") {
  df$Rewiring_Index[df$Rewiring_Index < 0.2] <- 0.2
  df$Rewiring_Index[df$Rewiring_Index > 0.6] <- 0.6
  df <- df %>% mutate(across(all_of(c(left, right)), as.character))

  edges <- df %>%
    transmute(
      from = paste0("L:", .data[[left]]),
      to = paste0("R:", .data[[right]]),
      edge_w = .data[[edge_width]],
      edge_c = .data[[edge_color]]
    )

  nodes_left <- df %>%
    distinct(label = .data[[left]]) %>%
    transmute(name = paste0("L:", label), label, side = "left", ord = row_number())
  nodes_right <- df %>%
    distinct(label = .data[[right]]) %>%
    transmute(name = paste0("R:", label), label, side = "right", ord = row_number())
  nodes <- bind_rows(nodes_left, nodes_right)
  layout <- nodes %>%
    group_by(side) %>%
    arrange(ord, .by_group = TRUE) %>%
    mutate(
      x = if_else(side == "left", 0, 1),
      y = rescale(-row_number(), to = c(1, 0))
    ) %>%
    ungroup() %>%
    select(name, x, y)

  graph <- graph_from_data_frame(edges, vertices = nodes, directed = FALSE)
  layout <- layout[match(V(graph)$name, layout$name), ]

  ggraph(graph, layout = "manual", x = layout$x, y = layout$y) +
    geom_edge_link(
      aes(edge_width = edge_w, edge_colour = edge_c),
      alpha = 0.8,
      lineend = "round"
    ) +
    scale_edge_width(
      range = c(0, 1.2),
      name = edge_width,
      breaks = c(0.2, 0.4, 0.6)
    ) +
    scale_edge_colour_gradient2(
      low = "#f1a340",
      mid = "#f7f7f7",
      high = "#6e016b",
      midpoint = 0.4,
      limits = c(0.2, 0.6),
      breaks = c(0.2, 0.4, 0.6),
      name = edge_color
    ) +
    geom_node_point(
      size = node_size,
      shape = 20,
      stroke = 0.8,
      fill = "white",
      colour = "black",
      show.legend = FALSE
    ) +
    geom_node_text(
      aes(label = label, hjust = if_else(side == "left", 1.5, -0.5)),
      size = 8
    ) +
    coord_cartesian(xlim = c(-0.25, 1.25)) +
    theme_void(base_size = 12) +
    theme(
      legend.position = "right",
      plot.margin = margin(10, 30, 10, 30)
    ) +
    annotate("text", x = -0.05, y = 1.04, label = left_title, hjust = 1, size = 4.2) +
    annotate("text", x = 1.05, y = 1.04, label = right_title, hjust = 0, size = 4.2)
}



plot_network <- function(sheet_name, output_file) {
  df <- read_source_data(sheet_name) %>%
    rename(
      TF_module = NTC,
      TF_module1 = KO,
      Rewiring_Index = `Disruption score`
    ) %>%
    mutate(
      TF_module = factor(TF_module, levels = 1:9),
      TF_module1 = factor(TF_module1, levels = 1:9),
      Rewiring_Index = as.numeric(Rewiring_Index)
    )

  p <- plot_lr_bipartite(
    df,
    left = "TF_module",
    right = "TF_module1",
    edge_width = "Rewiring_Index",
    edge_color = "Rewiring_Index",
    node_size = 4
  )

  ggsave(file.path(output_dir, output_file), p, width = 4.2, height = 5.2)

  p
}

h3b <- plot_tf_heatmap(
  sheet_name = "Fig3b",
  output_file = "Figure3b_TF_correlation_NTC.pdf",
  legend_title = "Correlation",
  color_limits = c(-0.6, -0.2, 0.2, 0.6),
  color_values = c("#114D72", "white", "white", "#78110C")
)

h3c <- plot_tf_heatmap(
  sheet_name = "Fig3c",
  output_file = "Figure3c_TF_correlation_DDX25_perturbed_cells.pdf",
  legend_title = "Correlation",
  color_limits = c(-0.6, -0.2, 0.2, 0.6),
  color_values = c("#114D72", "white", "white", "#78110C")
)

h3d <- plot_tf_heatmap(
  sheet_name = "Fig3d",
  output_file = "Figure3d_TF_correlation_change_DDX25_vs_NTC.pdf",
  legend_title = "Difference",
  color_limits = c(-0.5, 0, 0.5),
  color_values = c("#0F6056", "white", "#F25AA6")
)

figure3e <- read_source_data("Fig3e") %>%
  mutate(
    perturbation_module = factor(perturbation_module, levels = names(module_colors)),
    number_of_perturbations = as.numeric(number_of_perturbations),
    average_proportion_differential_tf_pairs = as.numeric(average_proportion_differential_tf_pairs)
  )

p3e <- ggplot(
  figure3e,
  aes(
    x = number_of_perturbations,
    y = average_proportion_differential_tf_pairs,
    color = perturbation_module,
    label = perturbation_module
  )
) +
  geom_point(size = 2.6) +
  geom_text_repel(size = 2.4, show.legend = FALSE, max.overlaps = Inf) +
  scale_color_manual(values = module_colors, drop = FALSE) +
  labs(
    x = "The number of perturbations within module",
    y = "Average proportion of significantly altered TF pairs",
    color = "Module"
  ) +
  theme_classic(base_size = 9)

ggsave(file.path(output_dir, "Figure3e_perturbation_module_tf_disruption.pdf"), p3e, width = 4.8, height = 4)


figure3f <- read_source_data("Fig3f") %>%
  mutate(
    perturbation_module = factor(perturbation_module, levels = names(module_colors)),
    proportion_altered_tfs = as.numeric(proportion_altered_tfs),
    proportion_differential_tf_pairs = as.numeric(proportion_differential_tf_pairs)
  )

figure3f_label <- figure3f %>%
  mutate(
    top_label_score = rescale(proportion_altered_tfs) +
      rescale(proportion_differential_tf_pairs)
  ) %>%
  arrange(desc(top_label_score)) %>%
  slice_head(n = 10)

p3f <- ggplot(
  figure3f,
  aes(
    x = proportion_altered_tfs,
    y = proportion_differential_tf_pairs,
    color = perturbation_module
  )
) +
  geom_point(size = 0.7, alpha = 0.9) +
  geom_text_repel(
    data = figure3f_label,
    aes(label = perturbation),
    size = 2.2,
    show.legend = FALSE,
    max.overlaps = Inf
  ) +
  scale_color_manual(values = module_colors, drop = FALSE) +
  coord_cartesian(ylim = c(0, 0.35)) +
  labs(
    x = "Proportion of altered TFs",
    y = "Proportion of differential TF pairs",
    color = "Module"
  ) +
  theme_classic(base_size = 9)

ggsave(file.path(output_dir, "Figure3f_perturbation_tf_disruption_scatter.pdf"), p3f, width = 5.2, height = 4)


network_panels <- c(
  Fig3g = "Figure3g_DDX25_TF_network_disruption.pdf",
  Fig3h = "Figure3h_IL1RL1_TF_network_disruption.pdf",
  Fig3i = "Figure3i_ASB11_TF_network_disruption.pdf",
  Fig3j = "Figure3j_CALB2_TF_network_disruption.pdf",
  Fig3k = "Figure3k_POSTN_TF_network_disruption.pdf"
)

network_plots <- lapply(
  names(network_panels),
  function(source_file) plot_network(source_file, network_panels[[source_file]])
)
names(network_plots) <- names(network_panels)

pdf("Figure3_updated_all_panels.pdf", width = 12, height = 10)
draw(h3b, heatmap_legend_side = "left", annotation_legend_side = "left", legend_grouping = "original")
draw(h3c, heatmap_legend_side = "left", annotation_legend_side = "left", legend_grouping = "original")
draw(h3d, heatmap_legend_side = "left", annotation_legend_side = "left", legend_grouping = "original")
print(ggarrange(p3e, p3f, ncol = 2, labels = c("Fig3e", "Fig3f")))
print(ggarrange(plotlist = network_plots, ncol = 3, labels = names(network_plots)))
dev.off()




