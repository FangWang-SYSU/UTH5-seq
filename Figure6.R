library(data.table)
library(reshape2)
library(ggplot2)
library(parallel)
library(dplyr)
library(ggpubr)
library(fmsb)
library(readxl)
library(pheatmap)
library(scales)
library(patchwork)

Sys.setlocale("LC_ALL", "en_US.UTF-8")

source_data_file <- "../source_data/Figure_6_Source_Data.xlsx"
if (!file.exists(source_data_file)) {
  stop("Cannot find source data file: ", source_data_file)
}

output_dir <- "Figure6_output"
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

radar_source <- function(sheet_name) {
  dat <- read_source(sheet_name)
  rownames(dat) <- dat[[1]]
  dat[[1]] <- NULL
  dat[] <- lapply(dat, as_numeric_column)
  dat
}

## Figure 6a
modules <- paste0("Module", 1:8)
module_plots <- list()
scale.col <- colorRampPalette(c("#fff7ec", "#fee8c8", "#fdbb84", "#fc8d59",
                                "#ef6548", "#d7301f", "#990000"))(16)
for (i in seq_along(modules)) {
  ms <- modules[i]
  plotdat <- read_source(paste0("Fig6a_M", i))
  plotdat[["UMAP 1"]] <- as_numeric_column(plotdat[["UMAP 1"]])
  plotdat[["UMAP 2"]] <- as_numeric_column(plotdat[["UMAP 2"]])
  plotdat[["Density"]] <- as_numeric_column(plotdat[["Density"]])

  breaks <- unique(as.numeric(quantile(plotdat[["Density"]],
                                       probs = seq(0.5, 1, length.out = length(scale.col) + 1),
                                       na.rm = TRUE)))
  pp <- ggplot(plotdat, aes(x = `UMAP 1`, y = `UMAP 2`)) +
    geom_contour_filled(
      aes(z = Density, fill = after_stat(level)),
      breaks = breaks,
      colour = "ivory"
    ) +
    scale_fill_manual(values = scale.col, guide = "none") +
    coord_equal() +
    ggtitle(ms) +
    theme_classic()

  module_plots[[ms]] <- pp
  ggsave(pp, file = file.path(output_dir, paste0("Figure6a_", ms, "_density.pdf")))
}
ggsave(wrap_plots(module_plots, ncol = 4),
       file = file.path(output_dir, "Figure6a_module_density_all.pdf"),
       width = 12,
       height = 6)

# Heatmap of module enrichment
moduleEnrich <- read_source("Fig6b")
rownames(moduleEnrich) <- moduleEnrich$Modules
moduleEnrich$Modules <- NULL
moduleEnrich[] <- lapply(moduleEnrich, as_numeric_column)
moduleEnrich <- t(as.matrix(moduleEnrich))
moduleEnrich <- moduleEnrich[c("S5", "S6", "S10", "S12", "S11", "S0", "S1", "S3", "S7", "S8", "S9", "S2", "S4"),
                             c("Module1", "Module2", "Module4", "Module3", "Module5", "Module6", "Module7", "Module8")]

moduleEnrich[moduleEnrich > 1.5] <- 1.5
moduleEnrich[moduleEnrich < (-1.5)] <- -1.5
pheatmap::pheatmap(moduleEnrich,
                   scale = "none",
                   cluster_rows = FALSE,
                   cluster_cols = FALSE,
                   border_color = "white",
                   color = colorRampPalette(c("#33ABC1", "white", "#B11927"))(20),
                   filename = file.path(output_dir, "Figure6b_module_enrichment_heatmap.pdf"))

## Top10 perturbation bubble plot
plotdat <- read_source("Fig6c")
plotdat$Module <- gsub("M_", "Module", plotdat$Module)
plotdat$log2ratio <- as_numeric_column(plotdat$log2ratio)
plotdat$log10pval <- as_numeric_column(plotdat$log10pval)
plotdat$Module <- factor(plotdat$Module, levels = rev(c("Module1", "Module2", "Module4", "Module3", "Module5", "Module6", "Module7", "Module8")))
plotdat <- plotdat[order(plotdat$Module), ]
plotdat$Perturbation <- factor(plotdat$Perturbation, levels = unique(plotdat$Perturbation))
plotdat$cluster <- as_cluster_id(plotdat$cluster)
plotdat$cluster <- factor(plotdat$cluster, levels = c("S5", "S6", "S10", "S12", "S11", "S0", "S1", "S3", "S7", "S8", "S9", "S2", "S4"))

cols <- c("#313772", "#2c4ca0", "#326db6", "#478ecc", "#75b5dc", "white", "#c44438", "#b7282e")
pp <- ggscatter(plotdat,
                x = "cluster",
                y = "Perturbation",
                color = "log2ratio",
                size = "log10pval",
                legend = "right") +
  scale_color_gradientn(
    colours = cols,
    limits = c(-6, 2),
    oob = squish
  )

ggsave(pp, file = file.path(output_dir, "Figure6c_top10_perturbation_bubble.pdf"), width = 6, height = 17)

### ACTG2 and MEGF8 function
enrichDat <- read_source("Fig6e")
enrichDat$Description <- factor(enrichDat$Description, levels = unique(enrichDat$Description))
enrichDat$log10qvalue <- as_numeric_column(enrichDat$log10qvalue)
p <- ggscatter(enrichDat,
               x = "groupby",
               y = "Description",
               color = "groups",
               size = "log10qvalue",
               legend = "right",
               palette = c("#45174E", "#D86B33")) +
  rotate_x_text(45)
ggsave(p, file = file.path(output_dir, "Figure6e_function_enrichment.pdf"))

plotdat <- radar_source("Fig6f")
pdf(file.path(output_dir, "Figure6f_MFGE8_up_radar.pdf"), width = 8, height = 8)
radarchart(plotdat,
           pty = c(16, 16, 32),
           axistype = 1,
           pcol = c("#008A89", "#D92927", "#3272AF"),
           plwd = c(3, 3, 3),
           plty = 1,
           cglcol = "grey60",
           cglty = 1,
           cglwd = 1,
           axislabcol = "grey60",
           vlcex = 0.8,
           vlabels = colnames(plotdat),
           caxislabels = c(0, 10, 20, 30, 40),
           calcex = 0.8)
legend(x = "bottomright",
       legend = c("S0", "S10", "S3"),
       horiz = FALSE,
       bty = "n",
       pch = 15,
       col = c("#008A89", "#D92927", "#3272AF"),
       text.col = "black",
       cex = 1, pt.cex = 1.5)
dev.off()

plotdat <- radar_source("Fig6g")
pdf(file.path(output_dir, "Figure6g_MFGE8_down_radar.pdf"), width = 8, height = 8)
radarchart(plotdat,
           pty = c(16, 16, 32),
           axistype = 1,
           pcol = c("#74C9BC", "#F79C7B", "#8FB4DC"),
           plwd = c(3, 3, 3),
           plty = 1,
           cglcol = "grey60",
           cglty = 1,
           cglwd = 1,
           axislabcol = "grey60",
           vlcex = 0.8,
           vlabels = colnames(plotdat),
           caxislabels = c(0, 10, 20, 30, 40),
           calcex = 0.8)
legend(x = "bottomright",
       legend = c("S0", "S10", "S3"),
       horiz = FALSE,
       bty = "n",
       pch = 15,
       col = c("#74C9BC", "#F79C7B", "#8FB4DC"),
       text.col = "black",
       cex = 1, pt.cex = 1.5)
dev.off()

plotdat <- radar_source("Fig6h")
pdf(file.path(output_dir, "Figure6h_EPB41L1_up_radar.pdf"), width = 8, height = 8)
radarchart(plotdat,
           pty = c(16, 16, 32),
           axistype = 1,
           pcol = c("#008A89", "#D92927"),
           plwd = c(3, 3, 3),
           plty = 1,
           cglcol = "grey60",
           cglty = 1,
           cglwd = 1,
           axislabcol = "grey60",
           vlcex = 0.8,
           vlabels = colnames(plotdat),
           caxislabels = c(0, 10, 20, 30, 40),
           calcex = 0.8)
legend(x = "bottomright",
       legend = c("S0", "S10"),
       horiz = FALSE,
       bty = "n",
       pch = 15,
       col = c("#008A89", "#D92927"),
       text.col = "black",
       cex = 1, pt.cex = 1.5)
dev.off()

plotdat <- radar_source("Fig6i")
pdf(file.path(output_dir, "Figure6i_EPB41L1_down_radar.pdf"), width = 8, height = 8)
radarchart(plotdat,
           pty = c(16, 16, 32),
           axistype = 1,
           pcol = c("#74C9BC", "#F79C7B"),
           plwd = c(3, 3, 3),
           plty = 1,
           cglcol = "grey60",
           cglty = 1,
           cglwd = 1,
           axislabcol = "grey60",
           vlcex = 0.8,
           vlabels = colnames(plotdat),
           caxislabels = c(0, 10, 20, 30, 40),
           calcex = 0.8)
legend(x = "bottomright",
       legend = c("S0", "S10"),
       horiz = FALSE,
       bty = "n",
       pch = 15,
       col = c("#74C9BC", "#F79C7B"),
       text.col = "black",
       cex = 1, pt.cex = 1.5)
dev.off()

plotdat <- radar_source("Fig6j")
pdf(file.path(output_dir, "Figure6j_GPR1_up_radar.pdf"), width = 8, height = 8)
radarchart(plotdat,
           pty = c(16, 16, 32),
           axistype = 1,
           pcol = c("#008A89", "#764FA0"),
           plwd = c(3, 3, 3),
           plty = 1,
           cglcol = "grey60",
           cglty = 1,
           cglwd = 1,
           axislabcol = "grey60",
           vlcex = 0.8,
           vlabels = colnames(plotdat),
           caxislabels = c(0, 10, 20, 30, 40),
           calcex = 0.8)
legend(x = "bottomright",
       legend = c("S0", "S1"),
       horiz = FALSE,
       bty = "n",
       pch = 15,
       col = c("#008A89", "#764FA0"),
       text.col = "black",
       cex = 1, pt.cex = 1.5)
dev.off()

plotdat <- radar_source("Fig6k")
pdf(file.path(output_dir, "Figure6k_GPR1_down_radar.pdf"), width = 8, height = 8)
radarchart(plotdat,
           pty = c(16, 16, 32),
           axistype = 1,
           pcol = c("#74C9BC", "#AB99C8"),
           plwd = c(3, 3, 3),
           plty = 1,
           cglcol = "grey60",
           cglty = 1,
           cglwd = 1,
           axislabcol = "grey60",
           vlcex = 0.8,
           vlabels = colnames(plotdat),
           caxislabels = c(0, 10, 20, 30, 40),
           calcex = 0.8)
legend(x = "bottomright",
       legend = c("S0", "S1"),
       horiz = FALSE,
       bty = "n",
       pch = 15,
       col = c("#74C9BC", "#AB99C8"),
       text.col = "black",
       cex = 1, pt.cex = 1.5)
dev.off()
