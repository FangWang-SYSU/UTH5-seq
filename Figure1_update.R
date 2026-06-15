library(ggplot2)
library(ggpubr)
library(dplyr)
library(scales)

source_data_dir <- if (dir.exists("../Figure_dat_for_public/Figure1")) {
  "../Figure_dat_for_public/Figure1"
} else if (dir.exists("Figure1")) {
  "Figure1"
} else {
  "."
}

output_dir <- file.path(source_data_dir, "Figure1_output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

read_source_data <- function(file_name) {
  read.csv(file.path(source_data_dir, file_name), header = TRUE, check.names = FALSE)
}

umi_colors <- c("#ccece6", "#238b45", "#006d2c", "#00441b")

plot_genome_umi_scatter <- function(data, title, legend_title) {
  ggscatter(
    data,
    x = "mouse_umi_thousands",
    y = "human_umi_thousands",
    size = 1,
    alpha = 0.7,
    shape = 16,
    color = "total_umi",
    ylab = "Human genome (UMI counts x 1000)",
    xlab = "Mouse genome (UMI counts x 1000)",
    legend = "right"
  ) +
    ggtitle(title) +
    theme_classic() +
    scale_color_gradientn(
      colours = umi_colors,
      name = legend_title,
      oob = squish
    )
}

plot_genome_umi_scatter_pair <- function(data_1, data_2, title_1, title_2) {
  ggarrange(
    plot_genome_umi_scatter(data_1, title_1, "UMI counts"),
    plot_genome_umi_scatter(data_2, title_2, "UMI counts"),
    ncol = 2
  )
}

plot_collision_pie <- function(data, title) {
  assignment_colors <- c(
    "Collision" = "grey",
    "Single index" = "#006d2c",
    "Non-multiplet" = "#006d2c"
  )
  assignment_levels <- names(assignment_colors)[names(assignment_colors) %in% data$index_assignment]

  pie_data <- data %>%
    mutate(
      label = sprintf("%s\n%.1f%%", index_assignment, proportion * 100),
      index_assignment = factor(index_assignment, levels = assignment_levels)
    )

  ggpie(
    pie_data,
    x = "cell_count",
    label = "label",
    fill = "index_assignment",
    lab.pos = "in",
    lab.font = c(4, "bold", "black"),
    palette = unname(assignment_colors[assignment_levels]),
    title = title,
    legend = "right"
  ) +
    font("title", face = "bold")
}

plot_index_assignment_pie_pair <- function(pie_data, title_1, title_2) {
  ggarrange(
    plot_collision_pie(filter(pie_data, condition == title_1), title_1),
    plot_collision_pie(filter(pie_data, condition == title_2), title_2),
    ncol = 2
  )
}

figure1b1_data <- read_source_data("Figure1b1_source_data.csv")
figure1b2_data <- read_source_data("Figure1b2_source_data.csv") %>%
  filter(total_umi > 1000)
figure1b_scatter <- plot_genome_umi_scatter_pair(
  figure1b1_data,
  figure1b2_data,
  title_1 = "Figure1b1",
  title_2 = "Figure1b2"
)
ggsave(
  file.path(output_dir, "Figure1b_genome_umi_scatter.pdf"),
  figure1b_scatter,
  width = 8,
  height = 4
)

figure1b_pie_data <- read_source_data("Figure1b_pie_source_data.csv")
figure1b_pie <- plot_index_assignment_pie_pair(
  figure1b_pie_data,
  title_1 = "Figure1b1",
  title_2 = "Figure1b2"
)
ggsave(
  file.path(output_dir, "Figure1b_index_assignment_pie.pdf"),
  figure1b_pie,
  width = 7,
  height = 3.5
)

figure1c1_data <- read_source_data("Figure1c1_source_data.csv")
figure1c2_data <- read_source_data("Figure1c2_source_data.csv")
figure1c_scatter <- plot_genome_umi_scatter_pair(
  figure1c1_data,
  figure1c2_data,
  title_1 = "Figure1c1",
  title_2 = "Figure1c2"
)
ggsave(
  file.path(output_dir, "Figure1c_genome_umi_scatter.pdf"),
  figure1c_scatter,
  width = 8,
  height = 4
)

figure1c_pie_data <- read_source_data("Figure1c_pie_source_data.csv")
figure1c_pie <- plot_index_assignment_pie_pair(
  figure1c_pie_data,
  title_1 = "Figure1c1",
  title_2 = "Figure1c2"
)
ggsave(
  file.path(output_dir, "Figure1c_index_assignment_pie.pdf"),
  figure1c_pie,
  width = 7,
  height = 3.5
)



figure1d_data <- read_source_data("Figure1d_source_data.csv")
figure1d_genes <- ggboxplot(
  figure1d_data,
  x = "cells_per_droplet",
  y = "detected_genes",
  xlab = "Cell number per droplet",
  ylab = "Number of genes",
  bxp.errorbar = FALSE,
  outliers = FALSE,
  outlier.shape = NA,
  width = 0.3,
  legend = "right",
  color = "#3969AC",
  palette = "igv"
)

figure1d_umis <- ggboxplot(
  figure1d_data,
  x = "cells_per_droplet",
  y = "total_umis",
  xlab = "Cell number per droplet",
  ylab = "UMI counts",
  bxp.errorbar = FALSE,
  outliers = FALSE,
  outlier.shape = NA,
  width = 0.3,
  legend = "right",
  color = "#3969AC",
  palette = "igv"
)

figure1d_boxplots <- ggarrange(figure1d_genes, figure1d_umis, ncol = 2)
ggsave(
  file.path(output_dir, "Figure1d_gene_umi_boxplots.pdf"),
  figure1d_boxplots,
  width = 7,
  height = 3.5
)

figure1f_data <- read_source_data("Figure1f_source_data.csv") %>%
  mutate(cell_count_100k = cell_count / 100000)
figure1f_sgrna_capture <- ggscatter(
  figure1f_data,
  x = "cell_count_100k",
  y = "sgrna_cell_percent",
  fill = "experiment",
  color = "experiment",
  palette = "npg",
  repel = TRUE,
  legend = "right",
  label = "experiment",
  label.select = figure1f_data$experiment,
  xlab = "Number of cells (x10^5)",
  ylab = "Cells with sgRNA (>=500 genes, %)",
  size = 3
)
ggsave(
  file.path(output_dir, "Figure1f_sgrna_capture_scatter.pdf"),
  figure1f_sgrna_capture,
  width = 6,
  height = 4
)

figure1g_data <- read_source_data("Figure1g_source_data.csv")
figure1g_cell_counts <- gghistogram(
  figure1g_data,
  x = "cell_count",
  bins = 20000,
  xlim = c(0, 1200),
  xlab = "# Number of cells for each perturbation",
  ylab = "# gRNA counts",
  col = NA,
  fill = "#1E4768",
  alpha = 0.8
) +
  geom_vline(
    xintercept = 50,
    linetype = "dashed",
    color = "darkgrey"
  )
ggsave(
  file.path(output_dir, "Figure1g_perturbation_cell_count_histogram.pdf"),
  figure1g_cell_counts,
  width = 3.5,
  height = 3
)
