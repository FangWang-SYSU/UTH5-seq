library(ggplot2)
library(data.table)
library(reshape2)
library(parallel)
library(Seurat)
library(reticulate)
library(anndata)
library(SeuratDisk)
library(readxl)
library(dplyr)
library(ggpubr)
library(patchwork)
library(fmsb)
library(pheatmap)
library(scales)

Sys.setlocale("LC_ALL", "en_US.UTF-8")

source_data_file <- "../source_data/Figure_7_Source_Data.xlsx"
if (!file.exists(source_data_file)) {
  stop("Cannot find source data file: ", source_data_file)
}

output_dir <- "Figure7_output"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

read_source <- function(sheet_name) {
  read_excel(source_data_file, sheet = sheet_name) %>%
    as.data.frame(check.names = FALSE)
}

as_numeric_column <- function(x) {
  as.numeric(as.character(x))
}

radar_source <- function(sheet_name) {
  df <- read_source(sheet_name)
  rownames(df) <- df[[1]]
  df[[1]] <- NULL
  df[] <- lapply(df, as_numeric_column)
  df
}

umapColor1 <- c("#7F3C8D", "#11A579", "#3969AC", "#E73F74", "#80BA5A",
                "#E68310", "#008695", "#CF1C90", "#f97b72", "#E7ABFD",
                "#9a221c", "#FED8A3", "#fcb93e", "#F9BEBB")
names(umapColor1) <- paste0("S", 0:13)
umapColor2 <- c("#53A85F", "#E5D2DD", "#E59CC4", "#D6E7A3", "#F1BB72", "#F3B1A0")
cols4 <- c("#7F3C8D", "#11A579", "#3969AC", "#E73F74", "#80BA5A",
           "#E68310", "#008695", "#CF1C90", "#f97b72", "#E7ABFD")

## Response clusters
perturbtions <- c("ACTG2", "MFGE8", "GPR1", "PTX4", "IL1RAPL1", "OLFM3")
for (pt in perturbtions) {
  plotdat <- read_source(paste0("Fig7a_", pt))
  plotdat$UMAP_1 <- as_numeric_column(plotdat$UMAP_1)
  plotdat$UMAP_2 <- as_numeric_column(plotdat$UMAP_2)
  plotdat$Subtype <- factor(plotdat$Subtype, levels = unique(plotdat$Subtype))
  plotdat$`Response cluster` <- factor(plotdat$`Response cluster`, levels = unique(plotdat$`Response cluster`))

  subtype_cols <- setNames(umapColor1[seq_along(levels(plotdat$Subtype))], levels(plotdat$Subtype))
  response_cols <- setNames(colorRampPalette(cols4)(length(levels(plotdat$`Response cluster`))),
                            levels(plotdat$`Response cluster`))

  p1 <- ggplot(plotdat, aes(UMAP_1, UMAP_2, color = Subtype)) +
    geom_point(alpha = 0.8, size = 0.4) +
    scale_color_manual(values = subtype_cols, na.translate = FALSE) +
    coord_equal() +
    theme_classic() +
    ggtitle(pt)
  p2 <- ggplot(plotdat, aes(UMAP_1, UMAP_2, color = `Response cluster`)) +
    geom_point(alpha = 0.8, size = 0.4) +
    scale_color_manual(values = response_cols, na.translate = FALSE) +
    coord_equal() +
    theme_classic() +
    ggtitle(pt)
  pp <- p1 / p2
  ggsave(pp, file = file.path(output_dir, paste0("Figure7a_", pt, "_scdat_mRNA.pdf")))
}

synergyScore <- read_source("Fig7b")
synergyScore$synergyScore <- as_numeric_column(synergyScore$synergyScore)
synergyScore$clusters <- factor(synergyScore$clusters, levels = paste0("S", 0:12))
synergyScore$groups <- synergyScore$Perturbations
p <- synergyScore %>%
  ggplot(aes(clusters, synergyScore)) +
  geom_violin(aes(fill = clusters),
              scale = "width",
              linewidth = 0) +
  facet_wrap(~groups,
             scales = "fixed",
             nrow = 3) +
  scale_fill_manual(values = umapColor1) +
  theme_classic2()

ggsave(p, file = file.path(output_dir, "Figure7b_synergy_score.pdf"))

genes <- c("IGFBP5", "COL8A1", "USP53", "CHSY3",
           "CITED2", "CTSB", "PRR16", "SOX5",
           "NRXN3", "PLCB4", "CHRM2")
df <- read_source("Fig7c")
df$exp <- as_numeric_column(df$expression)
df$mean <- as_numeric_column(df$mean_expression)
df$group <- factor(df$group,
                   levels = rev(c("NTC", "PTX4", "IL1RAPL1",
                                  "OLFM3", "GPR1", "OLFM3_PTX4",
                                  "IL1RAPL1_PTX4", "GPR1_PTX4")))

cols <- c("#fcc5c0", "#f768a1", "#dd3497", "#7a0177")
df$gene <- factor(df$gene, levels = genes)

p <- df %>%
  ggplot(aes(gene, exp)) +
  geom_violin(aes(fill = mean),
              scale = "width", linewidth = 0) +
  facet_wrap(~group, scales = "free_y", nrow = 10) +
  scale_fill_gradientn(colors = cols) +
  theme_classic2()
ggsave(p, file = file.path(output_dir, "Figure7c_shared_synergy_genes.pdf"))

genes <- c("GPR176", "KCNQ1OT1", "TBC1D19", "C5orf30",
           "LIMCH1", "CALM2", "MYH9", "MSN", "OXR1",
           "HEG1", "EIF4A2", "SDC2", "SKA2", "ENO1")
df <- read_source("Fig7d")
df$exp <- as_numeric_column(df$expression)
df$mean <- as_numeric_column(df$mean_expression)
df$group <- factor(df$group,
                   levels = rev(c("NTC", "PTX4", "IL1RAPL1",
                                  "OLFM3", "GPR1", "OLFM3_PTX4",
                                  "IL1RAPL1_PTX4", "GPR1_PTX4")))

cols <- c("#fcc5c0", "#f768a1", "#dd3497", "#7a0177")
df$gene <- factor(df$gene, levels = genes)

p <- df %>%
  ggplot(aes(gene, exp)) +
  geom_violin(aes(fill = mean),
              scale = "width", linewidth = 0) +
  facet_wrap(~group, scales = "free_y", nrow = 10) +
  scale_fill_gradientn(colors = cols) +
  theme_classic2()

ggsave(p, file = file.path(output_dir, "Figure7d_specific_synergy_genes.pdf"))

plotdat <- read_source("Fig7e")
plotdat$Description <- factor(plotdat$Description, levels = unique(plotdat$Description))
plotdat$qvalue <- as_numeric_column(plotdat$qvalue)
plotdat$Count <- as_numeric_column(plotdat$Count)
pp <- ggscatter(plotdat,
                x = "groups",
                y = "Description",
                color = "qvalue",
                size = "Count",
                legend = "right") +
  scale_color_gradientn(colors = cols) +
  scale_size_continuous(name = "Count", breaks = c(3, 5),
                        range = c(4, 6)) +
  rotate_x_text(45)
ggsave(pp, file = file.path(output_dir, "Figure7e_GO_BP_top10.pdf"))

df <- radar_source("Fig7f")
pdf(file.path(output_dir, "Figure7f_enrichment_KEGG_Reactome.pdf"), width = 8, height = 8)
radarchart(df[1:5, ],
           pty = c(16, 16, 16),
           axistype = 1,
           pcol = c("#4B114A", "#55998C", "#E4AD3D"),
           plwd = c(3, 3, 3),
           plty = 1,
           cglcol = "grey60",
           cglty = 1,
           cglwd = 1,
           axislabcol = "grey60",
           vlcex = 0.8,
           vlabels = colnames(df),
           caxislabels = c(0, 10, 20, 30, 40),
           calcex = 0.8)

legend(x = "topleft",
       legend = c("GPR1_PTX4", "IL1RAPL1_PTX4", "OLFM3_PTX4"),
       horiz = FALSE,
       bty = "n",
       pch = 15,
       col = c("#4B114A", "#55998C", "#E4AD3D"),
       text.col = "black",
       cex = 0.8,
       pt.cex = 1)
dev.off()

top_acts_mat <- radar_source("Fig7g")
colors <- rev(c("#8c510a", "white", "#35978f"))
colors.use <- grDevices::colorRampPalette(colors = colors)(100)
my_breaks <- c(seq(-0.3, 0, length.out = ceiling(100 / 2) + 1),
               seq(0.01, 0.3, length.out = floor(100 / 2)))

top_acts_mat <- data.frame(top_acts_mat)
top_acts_mat$Androgen <- NULL
top_acts_mat$Estrogen <- NULL
pheatmap::pheatmap(mat = as.matrix(top_acts_mat),
                   color = colors.use,
                   cluster_rows = FALSE,
                   border_color = "white",
                   breaks = my_breaks,
                   cellwidth = 15,
                   cellheight = 15,
                   treeheight_row = 20,
                   treeheight_col = 20,
                   filename = file.path(output_dir, "Figure7g_pathway_activity_heatmap.pdf"))
