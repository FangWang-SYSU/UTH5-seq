library(Seurat)
library(ggpubr)
library(Nebulosa)
library(readxl)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(scales)
library(patchwork)
library(ggalluvial)
library(ggforce)

Sys.setlocale("LC_ALL", "en_US.UTF-8")

source_data_file <- "../source_data/Figure_5_Source_Data.xlsx"
if (!file.exists(source_data_file)) {
  stop("Cannot find source data file: ", source_data_file)
}

output_dir <- "Figure5_output"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

read_source <- function(sheet_name) {
  read_excel(source_data_file, sheet = sheet_name) %>%
    as.data.frame(check.names = FALSE)
}

as_numeric_column <- function(x) {
  as.numeric(as.character(x))
}

as_cluster_id <- function(x) {
  x <- as.character(x)
  ifelse(grepl("^S", x), x, paste0("S", x))
}

umapColor1 <- c("#7F3C8D", "#11A579", "#3969AC", "#E73F74", "#80BA5A",
                "#E68310", "#008695", "#CF1C90", "#f97b72", "#E7ABFD",
                "#9a221c", "#FED8A3", "#fcb93e", "#F9BEBB")
names(umapColor1) <- paste0("S", 0:13)
cluster_levels <- c("S0", "S1", "S2", "S3", "S8", "S5", "S6", "S7",
                    "S10", "S11", "S4", "S9", "S12")
umapColor2 <- umapColor1[cluster_levels]

## Figure 5a
replot_data <- read_source("Fig5a")
replot_data[["Seurat cluster"]] <- as_cluster_id(replot_data[["Seurat cluster"]])
replot_data[["Seurat cluster"]] <- factor(replot_data[["Seurat cluster"]], levels = cluster_levels)
replot_data[["UMAP 1"]] <- as_numeric_column(replot_data[["UMAP 1"]])
replot_data[["UMAP 2"]] <- as_numeric_column(replot_data[["UMAP 2"]])

figure5a_plot <- ggplot(
  replot_data,
  aes(x = `UMAP 1`, y = `UMAP 2`, color = `Seurat cluster`)
) +
  geom_point(size = 0.25, stroke = 0, alpha = 1) +
  coord_equal() +
  scale_color_manual(values = umapColor2, na.translate = FALSE) +
  labs(x = "UMAP 1", y = "UMAP 2", color = "Seurat cluster") +
  theme_classic(base_size = 10) +
  theme(
    legend.key.height = unit(0.35, "cm"),
    legend.key.width = unit(0.35, "cm"),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8)
  )

ggsave(
  file.path(output_dir, "Figure5a_umap_clusters.pdf"),
  figure5a_plot,
  width = 4.2,
  height = 3.8,
  units = "in",
  useDingbats = FALSE
)

## Figure 5b marker gene expression
sc_dat.markers <- read_source("Fig5b")
colnames(sc_dat.markers)[1:2] <- c("gene", "cluster")
sc_dat.markers$cluster <- factor(sc_dat.markers$cluster, levels = cluster_levels)
markerexpression <- sc_dat.markers[, cluster_levels, drop = FALSE]
markerexpression[] <- lapply(markerexpression, as_numeric_column)
markerexpression <- as.matrix(markerexpression)
rownames(markerexpression) <- sc_dat.markers$gene
markerexpression <- markerexpression[order(sc_dat.markers$cluster), , drop = FALSE]
splitrow <- sc_dat.markers$cluster[order(sc_dat.markers$cluster)]

rowcol <- list(Clusters = umapColor2)
rowanno <- rowAnnotation(Clusters = splitrow, col = rowcol)
color <- colorRamp2(c(-2, -0.1, 0.1, 2), c("#5ab4ac", "white", "white", "#b2182b"))

pdf(file.path(output_dir, "Figure5b_marker_gene_expression.pdf"), width = 40, height = 20)
ht <- ComplexHeatmap::Heatmap(
  markerexpression,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  row_title = NULL,
  column_title = NULL,
  row_split = splitrow,
  left_annotation = rowanno,
  col = color,
  heatmap_legend_param = list(title = "Gene expression"),
  width = unit(10, "cm"),
  height = unit(30, "cm"),
  name = "correlation",
  column_names_rot = 0
)
draw(ht,
     heatmap_legend_side = "right",
     annotation_legend_side = "right",
     legend_grouping = "original")
dev.off()

## Figure 5c
plotdat <- read_source("Fig5c")
plotdat <- plotdat %>%
  transmute(
    groups = as_cluster_id(`Seurat cluster`),
    cytotrace_score = as_numeric_column(`CytoTRACE2 score`)
  )
plotdat <- plotdat %>%
  mutate(group = forcats::fct_reorder(groups, cytotrace_score, .fun = "median", .desc = TRUE))

p <- ggviolin(plotdat,
              x = "group",
              y = "cytotrace_score",
              color = "group",
              fill = "group",
              legend = "right",
              add = "boxplot",
              palette = umapColor1,
              add.params = list(fill = "white"),
              error.plot = "errorbar",
              border = "white")
ggsave(p, file = file.path(output_dir, "Figure5c_cluster_cytotrace_score.pdf"))

## Figure 5d/e
figure5d_source_data <- read_source("Fig5d")
figure5e_source_data <- read_source("Fig5e")
figure5d_source_data[["UMAP 1"]] <- as_numeric_column(figure5d_source_data[["UMAP 1"]])
figure5d_source_data[["UMAP 2"]] <- as_numeric_column(figure5d_source_data[["UMAP 2"]])
figure5d_source_data[["PRRX1 density"]] <- as_numeric_column(figure5d_source_data[["PRRX1 density"]])
figure5e_source_data[["UMAP 1"]] <- as_numeric_column(figure5e_source_data[["UMAP 1"]])
figure5e_source_data[["UMAP 2"]] <- as_numeric_column(figure5e_source_data[["UMAP 2"]])
figure5e_source_data[["ACTA2 density"]] <- as_numeric_column(figure5e_source_data[["ACTA2 density"]])

p1 <- ggplot(figure5d_source_data, aes(`UMAP 1`, `UMAP 2`, color = `PRRX1 density`)) +
  geom_point(size = 0.1) +
  scale_color_viridis_c(option = "plasma", name = "Density") +
  ggtitle("PRRX1") +
  theme_void()
p2 <- ggplot(figure5e_source_data, aes(`UMAP 1`, `UMAP 2`, color = `ACTA2 density`)) +
  geom_point(size = 0.1) +
  scale_color_viridis_c(option = "plasma", name = "Density") +
  ggtitle("ACTA2") +
  theme_void()
p <- p1 + p2
ggsave(p, file = file.path(output_dir, "Figure5d_e_density.pdf"))

## Figure 5f marker genes bubble plot
featureplot <- read_source("Fig5f")
featureplot$anno <- factor(featureplot$anno, levels = unique(featureplot$anno))
featureplot$features.plot <- featureplot$features
featureplot$features.plot <- factor(featureplot$features.plot,
                                    levels = rev(unique(featureplot$features.plot)))
featureplot$features <- featureplot$features.plot
featureplot$cluster <- as_cluster_id(featureplot$cluster)
featureplot$cluster <- factor(featureplot$cluster, levels = cluster_levels)
featureplot$avg.exp.scaled <- as_numeric_column(featureplot$avg.exp.scaled)
featureplot$pct_scaled <- as_numeric_column(featureplot$pct_scaled)

cols <- rev(c("#b2182b", "#ef8a62", "#fddbc7", "#e0e0e0", "#999999", "#4d4d4d"))

pp <- ggscatter(featureplot,
                x = "cluster",
                y = "features",
                color = "avg.exp.scaled",
                size = "pct_scaled",
                legend = "right") +
  scale_color_gradientn(
    colours = cols,
    limits = c(-2, 2),
    oob = squish
  )
ggsave(pp, file = file.path(output_dir, "Figure5f_feature_marker_bubble.pdf"))

## Figure 5k
prop <- read_source("Fig5k")
prop$Proportion <- as_numeric_column(prop$Proportion)
prop$clusters <- factor(prop$clusters, levels = cluster_levels)
p <- ggplot(prop, aes(x = Group, y = Proportion, fill = clusters, stratum = clusters,
                      alluvium = clusters)) +
  geom_col(width = 0.5, color = "white", size = 0.5) +
  geom_flow(width = 0.5, alpha = 0.3, knot.pos = 0.2, color = "white", size = 0.5) +
  coord_flip() +
  scale_fill_manual(values = umapColor2, na.translate = FALSE) +
  scale_y_continuous(expand = c(0, 0), name = "",
                     label = c("0%", "25%", "50%", "75%", "100%")) +
  scale_x_discrete(expand = c(0, 0), name = "") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(color = "black", size = 10),
        axis.ticks.length.x = unit(0.1, "cm"),
        plot.margin = margin(10, 10, 10, 10))
ggsave(p, file = file.path(output_dir, "Figure5k_cluster_proportion.pdf"), width = 6, height = 4)

## Figure 5g
cols <- c("#ffffcc", "#ffeda0", "#fed976", "#feb24c", "#fd8d3c",
          "#fc4e2a", "#e31a1c", "#bd0026", "#800026")
vln.df <- read_source("Fig5g")
vln.df <- vln.df %>%
  transmute(
    group = factor(Clusters, levels = cluster_levels),
    exp = as_numeric_column(scores),
    mean = as_numeric_column(mean_scores)
  )
p <- vln.df %>%
  ggplot(aes(group, exp)) +
  geom_violin(aes(fill = mean), scale = "width", linewidth = 0) +
  geom_boxplot(width = 0.3, outlier.shape = NA, color = "white") +
  scale_fill_gradientn(colors = cols) +
  scale_x_discrete("") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x.bottom = element_text(angle = 0, hjust = 1, vjust = 1, size = 8)
  ) +
  theme(strip.text.y = element_text(angle = 0)) +
  theme(
    strip.text = element_text(face = "bold", size = rel(0.5)),
    strip.background = element_rect(fill = "white", colour = "black", size = 1),
    strip.placement = "outside"
  )
ggsave(p, file = file.path(output_dir, "Figure5g_score_violin.pdf"), width = 5, height = 3)

## Figure 5h pathway bubble plot
pathway_bubble <- read_source("Fig5h")
pathway_bubble <- pathway_bubble %>%
  mutate(
    cluster = factor(as_cluster_id(cluster), levels = cluster_levels),
    log2FC = as_numeric_column(log2FC),
    log2FC_plot = pmin(pmax(log2FC, 0), 1.5),
    zscore = as_numeric_column(zscore),
    pathway_label = gsub("\\.+", " ", Pathways)
  )

pathway_levels <- rev(unique(pathway_bubble$pathway_label))
pathway_bubble <- pathway_bubble %>%
  mutate(
    pathway_label = factor(pathway_label, levels = pathway_levels),
    pathway_index = as.numeric(pathway_label),
    cluster_index = as.numeric(cluster)
  )

cluster_anno <- data.frame(
  cluster = factor(cluster_levels, levels = cluster_levels),
  cluster_index = seq_along(cluster_levels),
  group = c(rep("PRRX1+", 5), rep("ACTA2+", 5), rep("Others", 3)),
  stringsAsFactors = FALSE
)
cluster_anno$group <- factor(cluster_anno$group, levels = c("PRRX1+", "ACTA2+", "Others"))
cluster_anno_label <- cluster_anno %>%
  group_by(group) %>%
  summarise(cluster_index = mean(cluster_index), .groups = "drop")

p <- ggplot(pathway_bubble, aes(x = cluster_index, y = pathway_index)) +
  geom_tile(
    data = cluster_anno,
    aes(x = cluster_index, y = length(pathway_levels) + 1, fill = group),
    inherit.aes = FALSE,
    width = 0.9,
    height = 0.55
  ) +
  geom_text(
    data = cluster_anno_label,
    aes(x = cluster_index, y = length(pathway_levels) + 1.75, label = group),
    inherit.aes = FALSE,
    size = 3.5,
    fontface = "bold"
  ) +
  geom_point(aes(size = log2FC_plot, color = zscore), alpha = 0.95) +
  scale_fill_manual(values = c("PRRX1+" = "#C87983", "ACTA2+" = "#79AEB2", "Others" = "#D9A344"),
                    guide = "none") +
  scale_color_gradient2(
    low = "#2C5AA0",
    mid = "white",
    high = "#C43C7D",
    midpoint = 0,
    limits = c(-2, 2),
    oob = squish,
    name = "Z-score"
  ) +
  scale_size_continuous(range = c(0.35, 4.2), limits = c(0, 1.5), name = "log2FC") +
  scale_x_continuous(
    breaks = seq_along(cluster_levels),
    labels = cluster_levels,
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    breaks = seq_along(pathway_levels),
    labels = pathway_levels,
    position = "right",
    expand = expansion(add = c(0.5, 2.1))
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 8) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.ticks = element_line(linewidth = 0.25),
    legend.position = "right",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    plot.margin = margin(12, 8, 8, 8)
  )
ggsave(p, file = file.path(output_dir, "Figure5h_metabolism_bubble.pdf"),
       width = 7.2, height = 9.2, useDingbats = FALSE)

## Figure 5j NTC and perturbed CAF correlation triangle heatmap
corr_mat <- read_source("Fig5j")
rownames(corr_mat) <- corr_mat$Clusters
corr_mat$Clusters <- NULL
corr_mat[] <- lapply(corr_mat, as_numeric_column)

ntc_order <- paste0("NTC_", cluster_levels)
ko_order <- paste0("KO_", cluster_levels)
corr_mat <- as.matrix(corr_mat[ntc_order, ko_order])

triangle_df <- as.data.frame(as.table(corr_mat), stringsAsFactors = FALSE)
colnames(triangle_df) <- c("NTC", "KO", "Correlation")
triangle_df <- triangle_df %>%
  mutate(
    NTC_cluster = sub("^NTC_", "", NTC),
    KO_cluster = sub("^KO_", "", KO),
    row_index = match(NTC_cluster, cluster_levels),
    col_index = match(KO_cluster, cluster_levels)
  ) %>%
  filter(row_index >= col_index) %>%
  mutate(
    x = row_index + col_index,
    y = row_index - col_index,
    cell_id = row_number()
  )

triangle_poly <- triangle_df %>%
  slice(rep(row_number(), each = 4)) %>%
  group_by(cell_id) %>%
  mutate(
    corner = row_number(),
    x_poly = x + c(0, 0.9, 0, -0.9),
    y_poly = y + c(0.9, 0, -0.9, 0)
  ) %>%
  ungroup()

bottom_labels <- data.frame(
  label = cluster_levels,
  x = 2 * seq_along(cluster_levels),
  y = -0.9
)
left_labels <- data.frame(
  label = cluster_levels,
  x = seq_along(cluster_levels) + 1,
  y = seq_along(cluster_levels) - 1
)
right_labels <- data.frame(
  label = cluster_levels,
  x = length(cluster_levels) + seq_along(cluster_levels),
  y = length(cluster_levels) - seq_along(cluster_levels)
)

p <- ggplot(triangle_poly, aes(x = x_poly, y = y_poly, group = cell_id, fill = Correlation)) +
  geom_polygon(color = NA, linewidth = 0) +
  
  geom_text(data = left_labels, aes(x = x - 0.65, y = y, label = label),
            inherit.aes = FALSE, angle = 45, size = 3, hjust = 1) +
  geom_text(data = right_labels, aes(x = x + 0.65, y = y, label = label),
            inherit.aes = FALSE, angle = -45, size = 3, hjust = 0) +
  annotate("text", x = 2.2, y = length(cluster_levels) - 1.5, label = "NTCs",
           angle = 45, size = 3.4, fontface = "bold") +
  annotate("text", x = 2 * length(cluster_levels) - 2.2, y = length(cluster_levels) - 1.5,
           label = "Perturbed CAFs", angle = -45, size = 3.4, fontface = "bold") +
  scale_fill_gradientn(
    colors = c("#4FA66A", "#F7F7F2", "#F2A13A"),
    limits = c(0.5, 0.9),
    oob = squish,
    name = "Correlation"
  ) +
  coord_fixed(clip = "off") +
  theme_void(base_size = 9) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    plot.margin = margin(12, 20, 8, 20)
  )
ggsave(p, file = file.path(output_dir, "Figure5j_NTC_KO_triangle_heatmap.pdf"),
       width = 4.6, height = 4.0, useDingbats = FALSE)
