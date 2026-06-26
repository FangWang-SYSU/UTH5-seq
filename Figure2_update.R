library(data.table)
library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(readxl)
library(grid)

ht_opt$message <- FALSE

Sys.setlocale("LC_ALL", "en_US.UTF-8")

source_data_file <- "../source_data/Figure_2_Source_Data.xlsx"
if (!file.exists(source_data_file)) {
  stop("Cannot find source data file: ", source_data_file)
}

output_dir <- "Figure2_output"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

read_source_data <- function(sheet_name, ...) {
  as.data.table(read_excel(source_data_file, sheet = sheet_name, ...))
}

highlight_perturbations <- c(
  "DSCAML1", "FOXS1", "FOXH1", "IL1RL1", "LPS1",
  "PEG10", "TRAF3", "GATA4", "ASB11", "EPB41L1",
  "MFGE8", "POSTN", "INPPL1", "MAP2K3", "C1QL1",
  "C1QL4", "PTX3"
)

highlight_genes <- c(
  "CAV1", "CCL2", "MYC", "INPPL1", "STAT3", "COL6A1",
  "PTX3", "EGF", "NID2", "COL25A1", "IL6", "LOX", "FGF7",
  "COL4A5", "IL1RL1", "SMOC2", "POSTN", "ETV1", "TNC",
  "MMP3", "FOXO1", "PRRX1", "RUNX2", "GATA4",
  "FOXB1", "F3", "COL5A2", "AREG"
)

module_colors <- c(
  Module1 = "#003C30",
  Module2 = "#43AC5E",
  Module3 = "#2D4999",
  Module4 = "#C9B2D6",
  Module5 = "#E6780C",
  Module6 = "#7E9A29",
  Module7 = "#C51E7E",
  Module8 = "#723F91"
)

program_colors <- c(
  Program1 = "#4EBCD5",
  Program2 = "#E74C35",
  Program3 = "#07A087",
  Program4 = "#A75728",
  Program5 = "#F193A8",
  Program6 = "#BDBD21",
  Program7 = "#FABE70",
  Program8 = "#3D5589",
  Program9 = "#B19D86"
)

effect_size_data <- read_source_data("Fig2b")
gene_program_data <- read_source_data("Fig2c") %>%
  mutate(
    umap_1 = as.numeric(umap_1),
    umap_2 = as.numeric(umap_2)
  ) %>%
  as.data.table()
perturbation_module_data <- read_source_data("Fig2f")

effect_matrix <- as.matrix(effect_size_data[, -"perturbation"])
storage.mode(effect_matrix) <- "numeric"
rownames(effect_matrix) <- effect_size_data$perturbation
effect_matrix_raw <- effect_matrix

missing_genes <- setdiff(colnames(effect_matrix), gene_program_data$gene_symbol)
missing_perturbations <- setdiff(rownames(effect_matrix), perturbation_module_data$perturbation)
if (length(missing_genes) > 0) {
  stop("Figure2c_source_data.csv is missing heatmap genes: ", paste(head(missing_genes, 10), collapse = ", "))
}
if (length(missing_perturbations) > 0) {
  stop("Figure2f_source_data.csv is missing heatmap perturbations: ", paste(head(missing_perturbations, 10), collapse = ", "))
}

gene_program_data <- gene_program_data[match(colnames(effect_matrix), gene_symbol)]
perturbation_module_data <- perturbation_module_data[match(rownames(effect_matrix), perturbation)]

highlight_perturbation_index <- which(rownames(effect_matrix) %in% highlight_perturbations)
highlight_gene_index <- which(colnames(effect_matrix) %in% highlight_genes)

row_annotation <- rowAnnotation(
  Module = perturbation_module_data$perturbation_module,
  col = list(Module = module_colors),
  Perturbation = anno_mark(
    at = highlight_perturbation_index,
    labels = rownames(effect_matrix)[highlight_perturbation_index],
    labels_gp = gpar(fontsize = 10),
    lines_gp = gpar()
  )
)

column_annotation <- HeatmapAnnotation(
  Genes = anno_mark(
    at = highlight_gene_index,
    labels = colnames(effect_matrix)[highlight_gene_index],
    labels_gp = gpar(fontsize = 10),
    side = "top",
    lines_gp = gpar()
  ),
  Program = gene_program_data$gene_program,
  col = list(Program = program_colors)
)

heatmap_colors <- colorRamp2(
  c(-0.2, -0.05, 0, 0.05, 0.2),
  c("#00768B", "#79B8C6", "white", "#E7A15A", "#C36600")
)
effect_matrix[effect_matrix > 0.2] <- 0.2
effect_matrix[effect_matrix < -0.2] <- -0.2

effect_heatmap <- ComplexHeatmap::Heatmap(
  effect_matrix,
  show_row_names = FALSE,
  name = "Effect size",
  show_column_names = FALSE,
  show_column_dend = FALSE,
  cluster_columns = FALSE,
  cluster_row_slices = TRUE,
  cluster_rows = FALSE,
  col = heatmap_colors,
  use_raster = TRUE,
  raster_quality = 4,
  raster_device = "png",
  row_split = perturbation_module_data$perturbation_module,
  column_split = gene_program_data$gene_program,
  top_annotation = column_annotation,
  right_annotation = row_annotation
)
draw_effect_heatmap <- function() {
  draw(
    effect_heatmap,
    show_heatmap_legend = TRUE,
    heatmap_legend_side = "left",
    annotation_legend_side = "left",
    legend_grouping = "original"
  )
}

pdf(file.path(output_dir, "Figure2b_effect_size_heatmap.pdf"), width = 25, height = 20)
draw_effect_heatmap()
dev.off()

gene_correlation_histogram_data <- read_source_data("Fig2d") %>%
  mutate(
    bin_midpoint = as.numeric(bin_midpoint),
    bin_width = as.numeric(bin_width),
    density = as.numeric(density)
  )
figure2d_correlation_histogram <- ggplot(
  gene_correlation_histogram_data,
  aes(x = bin_midpoint, y = density, fill = comparison_group, color = comparison_group)
) +
  geom_col(aes(width = bin_width), alpha = 0.55, position = "identity") +
  scale_fill_manual(values = c("inter-program" = "#2078B5", "intra-program" = "#D72829")) +
  scale_color_manual(values = c("inter-program" = "#2078B5", "intra-program" = "#D72829")) +
  labs(x = "Correlation of expression", y = "Density", fill = NULL, color = NULL) +
  theme_classic()
ggsave(
  file.path(output_dir, "Figure2d_gene_program_correlation_histogram.pdf"),
  figure2d_correlation_histogram,
  width = 4,
  height = 3
)

effect_density_data <- read_source_data("Fig2e") %>%
  mutate(absolute_effect_size = as.numeric(absolute_effect_size))
figure2e_effect_density <- ggdensity(
  effect_density_data,
  x = "absolute_effect_size",
  color = "comparison_group",
  palette = c("#91B1DE", "#F18516"),
  legend = "right",
  xlab = "Absolute effect size",
  ylab = "Density"
)
ggsave(
  file.path(output_dir, "Figure2e_effect_size_density.pdf"),
  figure2e_effect_density,
  width = 4,
  height = 3
)

show_genes <- c(
  "DHX9", "YWHAQ", "MCM6", "SF3B1", "NCL", "MLH1", "POLR2B",
  "PPP2CA", "NPM1", "HDAC2", "CUL1", "YWHAZ", "TFDP1",
  "HNRNPC", "SRSF1", "CAPZB", "LMNA", "ARF1", "DCTN1",
  "PSMB1", "CYCS", "YWHAG", "GSN", "TUBA1A", "CTSA", "ATP2B4", "LPIN1",
  "HADHA", "HADHB", "KLKB1", "HSD17B4", "TFPI2", "MYOF",
  "CD9", "SMAD3", "WNT2B", "DDR2", "TGFB2", "COL3A1", "WNT5A", "PLOD2",
  "FGF2", "LOX", "SOX4", "DAAM2", "COL12A1", "COL10A1",
  "CDC42", "LYST", "MAP4K3", "KIF13A", "PRKN", "NTRK2", "FGF7",
  "VPS13C", "VPS35", "ITGA9", "CCK", "SUCNR1", "CXCL13",
  "SCIN", "CNTNAP2", "VIPR2", "TACR2", "BDKRB1", "MAP2K1", "CLDN6",
  "ITGAM", "PVR", "PDGFB", "KIT", "PDGFC", "VEGFC", "GRIA1",
  "EGFR", "SHANK2", "SYT1", "IGF1", "NTRK3", "DLGAP1", "NEGR1",
  "ERBB4", "ADGRL3", "EREG", "PTPN13", "GRID2", "RASA1", "ELMO1",
  "DOCK1", "HIF1A", "PRKACB", "CDC42EP3", "ERLEC1", "GCC2",
  "GOLGA4", "AREG", "DKK2", "PEBP1", "PPP2R5E", "CAV1", "CCL2",
  "MYC", "INPPL1", "STAT3", "COL6A1", "PTX3", "EGF", "NID2",
  "COL25A1", "IL6", "COL4A5", "IL1RL1", "SMOC2", "POSTN",
  "ETV1", "TNC", "MMP3", "FOXO1", "PRRX1", "RUNX2", "GATA4", "FOXB1",
  "F3", "COL5A2"
)

figure2c_gene_program_umap <- ggscatter(
  gene_program_data,
  x = "umap_1",
  y = "umap_2",
  color = "gene_program",
  size = 0.3,
  label = "gene_symbol",
  label.select = show_genes,
  palette = unname(program_colors)
) +
  labs(x = "UMAP_1", y = "UMAP_2", color = NULL) +
  theme_classic()
ggsave(
  file.path(output_dir, "Figure2c_gene_program_umap.pdf"),
  figure2c_gene_program_umap,
  width = 6,
  height = 5
)

figure2f_perturbation_module_umap <- ggscatter(
  perturbation_module_data,
  x = "umap_1",
  y = "umap_2",
  color = "perturbation_module",
  size = 2,
  label = "perturbation",
  label.select = highlight_perturbations,
  palette = unname(module_colors)
) +
  labs(x = "UMAP_1", y = "UMAP_2", color = NULL) +
  theme_classic()
ggsave(
  file.path(output_dir, "Figure2f_perturbation_module_umap.pdf"),
  figure2f_perturbation_module_umap,
  width = 6,
  height = 5
)

pdf("Figure2_updated_all_panels.pdf", width = 25, height = 20)
draw_effect_heatmap()
print(
  ggarrange(
    figure2c_gene_program_umap,
    figure2d_correlation_histogram,
    figure2e_effect_density,
    figure2f_perturbation_module_umap,
    ncol = 2,
    nrow = 2,
    labels = c("Fig2c", "Fig2d", "Fig2e", "Fig2f")
  )
)
dev.off()
