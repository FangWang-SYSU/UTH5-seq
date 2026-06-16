suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(ggplot2)
  library(ggraph)
  library(ggrepel)
  library(igraph)
  library(scales)
})

script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
script_dir <- if (length(script_path) == 0) getwd() else dirname(normalizePath(script_path))
project_dir <- normalizePath(file.path(script_dir, ".."))
data_dir <- file.path(project_dir, "Figure_dat", "Figure3")
output_dir <- file.path(data_dir, "Figure3_output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

module_colors <- c(
  Module1 = "#003C30", Module2 = "#43AC5E", Module3 = "#2D4999",
  Module4 = "#C9B2D6", Module5 = "#E6780C", Module6 = "#7E9A29",
  Module7 = "#C51E7E", Module8 = "#723F91"
)

read_matrix_source <- function(file_name) {
  df <- read.csv(file.path(data_dir, file_name), check.names = FALSE)
  row_ids <- df$transcription_factor
  mat <- as.matrix(df[, setdiff(names(df), "transcription_factor"), drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- row_ids
  mat
}

plot_tf_heatmap <- function(matrix_file, output_file, legend_title, color_limits, color_values) {
  tf_annotation <- read.csv(
    file.path(data_dir, "Figure3bcd_tf_module_annotation_source_data.csv"),
    check.names = FALSE
  )
  mat <- read_matrix_source(matrix_file)
  ordered_tfs <- intersect(tf_annotation$transcription_factor, rownames(mat))
  mat <- mat[ordered_tfs, ordered_tfs]
  tf_module <- tf_annotation$tf_module[match(ordered_tfs, tf_annotation$transcription_factor)]

  clipped_mat <- mat
  clipped_mat[clipped_mat < min(color_limits)] <- min(color_limits)
  clipped_mat[clipped_mat > max(color_limits)] <- max(color_limits)

  heatmap_colors <- colorRamp2(color_limits, color_values)
  row_anno <- rowAnnotation(TF_module = tf_module)
  col_anno <- HeatmapAnnotation(TF_module = tf_module)

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
}

plot_lr_bipartite <- function(df,
                              left = "TF_module",
                              right = "TF_module1",
                              edge_width = "tf_prop",
                              edge_color = "meanDiffcors",
                              node_size = 3.8,
                              left_title = "NTCs",
                              right_title = "Targeted cells") {
  df$Rewiring_Index[df$Rewiring_Index < 0.2] <- 0
  df$Rewiring_Index[df$Rewiring_Index > 0.7] <- 0.7
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
      range = c(0, 1.6),
      name = edge_width,
      breaks = c(0.3, 0.5, 0.7)
    ) +
    scale_edge_colour_gradient2(
      low = "#f1a340",
      mid = "#f7f7f7",
      high = "#6e016b",
      midpoint = 0.35,
      limits = c(0, 0.7),
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

plot_network <- function(source_file, output_file) {
  df <- read.csv(file.path(data_dir, source_file), check.names = FALSE) %>%
    mutate(
      TF_module = factor(TF_module, levels = 1:9),
      TF_module1 = factor(TF_module1, levels = 1:9)
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
}

plot_tf_heatmap(
  matrix_file = "Figure3b_source_data.csv",
  output_file = "Figure3b_TF_correlation_NTC.pdf",
  legend_title = "Correlation",
  color_limits = c(-0.6, -0.2, 0.2, 0.6),
  color_values = c("#114D72", "white", "white", "#78110C")
)

plot_tf_heatmap(
  matrix_file = "Figure3c_source_data.csv",
  output_file = "Figure3c_TF_correlation_DDX25_perturbed_cells.pdf",
  legend_title = "Correlation",
  color_limits = c(-0.6, -0.2, 0.2, 0.6),
  color_values = c("#114D72", "white", "white", "#78110C")
)

plot_tf_heatmap(
  matrix_file = "Figure3d_source_data.csv",
  output_file = "Figure3d_TF_correlation_change_DDX25_vs_NTC.pdf",
  legend_title = "Difference",
  color_limits = c(-0.5, 0, 0.5),
  color_values = c("#0F6056", "white", "#F25AA6")
)

figure3e <- read.csv(file.path(data_dir, "Figure3e_source_data.csv"), check.names = FALSE) %>%
  mutate(perturbation_module = factor(perturbation_module, levels = names(module_colors)))

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

figure3f <- read.csv(file.path(data_dir, "Figure3f_source_data.csv"), check.names = FALSE) %>%
  mutate(
    perturbation_module = factor(perturbation_module, levels = names(module_colors)),
    label = ifelse(is.na(label), "", label)
  )

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
    data = filter(figure3f, label != ""),
    aes(label = label),
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
  Figure3g_DDX25_tf_network_source_data.csv = "Figure3g_DDX25_TF_network_disruption.pdf",
  Figure3h_IL1RL1_tf_network_source_data.csv = "Figure3h_IL1RL1_TF_network_disruption.pdf",
  Figure3i_ASB11_tf_network_source_data.csv = "Figure3i_ASB11_TF_network_disruption.pdf",
  Figure3j_CALB2_tf_network_source_data.csv = "Figure3j_CALB2_TF_network_disruption.pdf",
  Figure3k_POSTN_tf_network_source_data.csv = "Figure3k_POSTN_TF_network_disruption.pdf"
)
main_network_genes <- c("DDX25", "IL1RL1", "ASB11", "CALB2", "POSTN")
supplementary_network_genes <- c("TSPYL6", "TREM1", "MAPK10", "AUTS2", "NPTN")

for (source_file in names(network_panels)) {
  plot_network(source_file, network_panels[[source_file]])
for (gene in main_network_genes) {
  plot_network(
    paste0("Figure3_network_", gene, "_source_data.csv"),
    paste0("Figure3_network_", gene, "_TF_network_disruption.pdf")
  )
}

for (gene in supplementary_network_genes) {
  plot_network(
    paste0("Figure3_network_", gene, "_source_data.csv"),
    paste0("Figure3_supp_network_", gene, "_TF_network_disruption.pdf")
  )
}
