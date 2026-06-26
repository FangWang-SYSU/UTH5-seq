library(tidyverse)
library(ggraph)
library(igraph)
library(scales)
library(ggplot2)
library(ggpubr)
library(Seurat)
library(ComplexHeatmap)
library(colorRamps)
library(circlize) 
library(dplyr)
library(reshape2)
library(data.table)
library(ggplot2)
library(dplyr)
library(readxl)
library(linkET)

Sys.setlocale("LC_ALL", "en_US.UTF-8")

source_data_file <- "../source_data/Figure_4_Source_Data.xlsx"
if (!file.exists(source_data_file)) {
  stop("Cannot find source data file: ", source_data_file)
}

output_dir <- "Figure4_output"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

read_xlsx_sheet <- function(sheet_name) {
  read_excel(source_data_file, sheet = sheet_name) %>%
    as.data.frame(check.names = FALSE)
}

as_numeric_column <- function(x) {
  as.numeric(as.character(x))
}

read_source <- function(file_name) {
  sheet_map <- c(
    "Figure4a_source_data.csv" = "Fig4a",
    "Figure4b_source_data.csv" = "Fig4b",
    "Figure4c_source_data.csv" = "Fig4c"
  )

  sheet_name <- unname(sheet_map[file_name])
  if (is.na(sheet_name)) {
    stop("No source-data sheet mapping for: ", file_name)
  }

  dat <- read_xlsx_sheet(sheet_name)
  if (sheet_name == "Fig4a") {
    dat$CytoTRACE2_Score <- as_numeric_column(dat$CytoTRACE2_Score)
  }
  if (sheet_name == "Fig4b") {
    dat$proliferation <- as_numeric_column(dat$proliferation)
  }
  if (sheet_name == "Fig4c") {
    dat$matrix_remodling1 <- as_numeric_column(dat$ECM)
  }
  dat
}

write_legacy_scatter_csv <- function(sheet_name, condition, output_file, rename_cols = NULL) {
  dat <- read_xlsx_sheet(sheet_name) %>%
    filter(Condition == condition) %>%
    select(-Condition)

  if (!is.null(rename_cols)) {
    names(dat) <- rename_cols
  }

  dat[] <- lapply(dat, as_numeric_column)
  write.csv(dat, output_file, quote = FALSE)
}

write_legacy_expression_csv <- function(sheet_name, gene_name, output_file) {
  dat <- read_xlsx_sheet(sheet_name) %>%
    filter(Gene == gene_name) %>%
    transmute(
      group = Condition,
      expression = as_numeric_column(Expression)
    )
  write.csv(dat, output_file, quote = FALSE)
}

materialize_legacy_inputs <- function() {
  write_legacy_scatter_csv("Fig4f", "NTC", "GLIS1_MECOM_NTC.csv", c("GLIS1", "MECOM"))
  write_legacy_scatter_csv("Fig4f", "POSTN_CRISPR_KO", "GLIS1_MECOM_PODTN_KO.csv", c("GLIS1", "MECOM"))
  write_legacy_scatter_csv("Fig4g", "NTC", "ZNF117_MECOM_NTC.csv", c("ZNF117", "MECOM"))
  write_legacy_scatter_csv("Fig4g", "POSTN_CRISPR_KO", "ZNF117_MECOM_POSTN_KO.csv", c("ZNF117", "MECOM"))
  write_legacy_scatter_csv("Fig4h", "NTC", "MYBL1_NR2C2_NTC.csv", c("MYBL1", "NR2C2"))
  write_legacy_scatter_csv("Fig4h", "POSTN_CRISPR_KO", "MYBL1_NR2C2_POSTN_KO.csv", c("MYBL1", "NR2C2"))
  write_legacy_scatter_csv("Fig4i", "NTC", "MYBL2_NR2C2_NTC.csv", c("MYBL2", "NR2C2"))
  write_legacy_scatter_csv("Fig4i", "POSTN_CRISPR_KO", "MYBL2_NR2C2_POSTN_KO.csv", c("MYBL2", "NR2C2"))

  write_legacy_scatter_csv("Fig4l", "NTC", "GATA6_ZNF608_NTC.csv", c("GATA6", "ZNF608"))
  write_legacy_scatter_csv("Fig4l", "ACTG2_KO", "GATA6_ZNF608_ACTG2_KO.csv", c("GATA6", "ZNF608"))
  write_legacy_expression_csv("Fig4m", "GATA6", "GATA6_expression.csv")
  write_legacy_expression_csv("Fig4m", "ZNF608", "ZNF608_expression.csv")

  write_legacy_scatter_csv("Fig4n", "NTC", "RARB_TWIST2_NTC.csv", c("RARB", "TWIST2"))
  write_legacy_scatter_csv("Fig4n", "KRTAP3-1_KO", "RARB_TWIST2_KRTAP3-1_KO.csv", c("RARB", "TWIST2"))
  write_legacy_expression_csv("Fig4o", "RARB", "RARB_expression.csv")
  write_legacy_expression_csv("Fig4o", "TWIST2", "TWIST2_expression.csv")

  fig4a <- read_source("Figure4a_source_data.csv")
  fig4b <- read_source("Figure4b_source_data.csv")
  fig4c <- read_source("Figure4c_source_data.csv")
  ntc <- data.frame(
    CytoTRACE2_Score = fig4a$CytoTRACE2_Score,
    proliferationScore = fig4b$proliferation,
    matrix_remodling1 = fig4c$matrix_remodling1,
    group = fig4a$Condition,
    check.names = FALSE
  )
  saveRDS(ntc, "postn_proliferation.rds")

  tf6 <- read_xlsx_sheet("Fig4e right")
  tf6_names <- tf6$TF
  ntc_mat <- tf6 %>%
    rename(RowName = TF) %>%
    mutate(across(-RowName, as_numeric_column)) %>%
    as.data.frame(check.names = FALSE)
  rownames(ntc_mat) <- ntc_mat$RowName
  ntc_mat$RowName <- NULL

  top_ntc <- read_xlsx_sheet("Fig4e_top") %>%
    mutate(Correlation = as_numeric_column(Correlation))
  left_ko <- read_xlsx_sheet("Fig4_left") %>%
    mutate(Correlation = as_numeric_column(Correlation))
  tf7_names <- setdiff(unique(c(top_ntc$TF1, top_ntc$TF2, left_ko$TF1, left_ko$TF2)), tf6_names)
  all_tfs <- unique(c(tf6_names, tf7_names))

  make_square <- function() {
    mat <- matrix(0, nrow = length(all_tfs), ncol = length(all_tfs), dimnames = list(all_tfs, all_tfs))
    common <- intersect(rownames(ntc_mat), colnames(ntc_mat))
    mat[common, common] <- as.matrix(ntc_mat[common, common, drop = FALSE])
    mat
  }

  ntc_square <- make_square()
  ko_square <- make_square()

  for (i in seq_len(nrow(top_ntc))) {
    ntc_square[top_ntc$TF1[i], top_ntc$TF2[i]] <- top_ntc$Correlation[i]
    ntc_square[top_ntc$TF2[i], top_ntc$TF1[i]] <- top_ntc$Correlation[i]
  }
  for (i in seq_len(nrow(left_ko))) {
    ko_square[left_ko$TF1[i], left_ko$TF2[i]] <- left_ko$Correlation[i]
    ko_square[left_ko$TF2[i], left_ko$TF1[i]] <- left_ko$Correlation[i]
  }

  tf_cluster <- data.frame(
    RowName = all_tfs,
    Group = ifelse(all_tfs %in% tf6_names, "6", "7"),
    check.names = FALSE
  )

  write.csv(tf_cluster, "TF_clusters_NTC.csv", quote = FALSE)
  write.csv(ntc_square, "TF_correlation_NTC.csv", quote = FALSE)
  write.csv(ko_square, "TF_correlation_POSTN_KO.csv", quote = FALSE)
}

materialize_legacy_inputs()

melt <- function(x, ...) {
  df <- as.data.frame(x, check.names = FALSE)
  id_values <- df$TF1
  value_cols <- setdiff(names(df), "TF1")
  data.frame(
    TF1 = rep(id_values, times = length(value_cols)),
    variable = rep(value_cols, each = nrow(df)),
    value = as.vector(as.matrix(df[, value_cols, drop = FALSE])),
    check.names = FALSE
  )
}

qcorrplot <- function(correlate, type = "full", diag = TRUE, ...) {
  mat <- as.matrix(correlate)
  storage.mode(mat) <- "numeric"
  row_levels <- rownames(mat)
  col_levels <- colnames(mat)
  plot_data <- expand.grid(
    y = row_levels,
    x = col_levels,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  plot_data$r <- as.vector(mat)
  plot_data$row_index <- match(plot_data$y, row_levels)
  plot_data$col_index <- match(plot_data$x, col_levels)

  if (!diag) {
    plot_data <- plot_data[plot_data$row_index != plot_data$col_index, , drop = FALSE]
  }
  if (type == "upper") {
    plot_data <- plot_data[plot_data$row_index < plot_data$col_index, , drop = FALSE]
  }
  if (type == "lower") {
    plot_data <- plot_data[plot_data$row_index > plot_data$col_index, , drop = FALSE]
  }

  ggplot(plot_data, aes(x = factor(x, levels = col_levels), y = factor(y, levels = rev(row_levels)))) +
    coord_fixed() +
    labs(x = NULL, y = NULL) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.ticks = element_blank()
    )
}

cols=c("#084594","#2c4ca0", "#4292c6", 
       "#9ecae1", "#deebf7","#f7fbff",
       "white", "white","#fee0d2","#fcbba1"
       ,"#fc9272","#ef3b2c", "#cb181d","#99000d")

## Figure 4a
fig4a <- read_source("Figure4a_source_data.csv")
plotdat <- data.frame(gene_score = fig4a$CytoTRACE2_Score, groups = fig4a$Condition)
plotdat <- plotdat %>% mutate(group = forcats::fct_reorder(groups, gene_score, .fun = 'median', .desc = FALSE))
p <- ggviolin(plotdat,
              x = "group",
              y = "gene_score",
              color = 'group',
              fill = "group",
              legend = "right",
              outlier = FALSE,
              add = "boxplot",
              palette = c("#eb7e60", "grey"),
              add.params = list(fill = "white", outlier.shape = NA),
              error.plot = "errorbar", border = "white") + stat_compare_means()
ggsave(p, file = file.path(output_dir, "Figure4a_cytotrace_score.pdf"))

## Figure 4b
fig4b <- read_source("Figure4b_source_data.csv")
plotdat <- data.frame(gene_score = fig4b$proliferation, groups = fig4b$Condition)
plotdat <- plotdat %>% mutate(group = forcats::fct_reorder(groups, gene_score, .fun = 'median', .desc = FALSE))
p <- ggviolin(plotdat,
              x = "group",
              y = "gene_score",
              color = 'group',
              fill = "group",
              legend = "right",
              outlier = FALSE,
              add = "boxplot",
              palette = c("#eb7e60", "grey"),
              add.params = list(fill = "white", outlier.shape = NA),
              error.plot = "errorbar", border = "white") + stat_compare_means()
ggsave(p, file = file.path(output_dir, "Figure4b_proliferation_score.pdf"))

## Figure 4c
fig4c <- read_source("Figure4c_source_data.csv")
plotdat <- data.frame(gene_score = fig4c$matrix_remodling1, groups = fig4c$Condition)
plotdat <- plotdat %>% mutate(group = forcats::fct_reorder(groups, gene_score, .fun = 'median', .desc = FALSE))
p <- ggviolin(plotdat,
              x = "group",
              y = "gene_score",
              color = 'group',
              fill = "group",
              legend = "right",
              outlier = FALSE,
              add = "boxplot",
              palette = c("#eb7e60", "grey"),
              add.params = list(fill = "white", outlier.shape = NA),
              error.plot = "errorbar", border = "white") + stat_compare_means()
ggsave(p, file = file.path(output_dir, "Figure4c_matrix_remodling_score.pdf"))

postn_tf_cors <- read_source("Figure4c_source_data.csv")

postn_tf_cors<-read.csv("TF_correlation_POSTN_KO.csv",
                        header = T,
                        row.names = 1,
                        check.names = F)

NTC_TF_cluster<-read.csv("TF_clusters_NTC.csv",
                         header = T,
                         row.names = 1,
                         check.names = F)

NTC_TF_correlation<-read.csv("TF_correlation_NTC.csv",
                             header = T,
                             row.names = 1,
                             check.names = F)

NTC_TF_correlation=NTC_TF_correlation[NTC_TF_cluster$RowName,NTC_TF_cluster$RowName]
postn_tf_cors=postn_tf_cors[NTC_TF_cluster$RowName,NTC_TF_cluster$RowName]


tf_cor_difference<-postn_tf_cors-NTC_TF_correlation

c6=subset(NTC_TF_cluster,Group=="6")
c7=subset(NTC_TF_cluster,Group=="7")

mat<-tf_cor_difference[unique(c(c6$RowName,c7$RowName)),unique(c(c6$RowName,c7$RowName))]
mat$TF1=rownames(mat)
mat=melt(mat)
colnames(mat)<-c("TF1","TF2","cors")
mat$TF_pair <- apply(mat[, c("TF1", "TF2")], 1, function(x) paste(sort(x), collapse = "_"))
setDT(mat)
mat[,c("TF1","TF2"):=tstrsplit(TF_pair,"_")]
mat$groups="TFc6-7"
mat$groups[mat$TF1%in%c(c6$RowName) &mat$TF2%in%c(c6$RowName) ]="TFc6-6"
mat$groups[mat$TF1%in%c(c7$RowName) &mat$TF2%in%c(c7$RowName) ]="TFc7-7"
mat <- unique(mat)

write.csv(mat,file="POSTN_KO_6-7-network.csv",quote = F)


## POSTN KO heatmap of TFc6 and TFc7
postn_ko_6=postn_tf_cors[c6$RowName,c6$RowName]
ntc_6=NTC_TF_correlation[c6$RowName,c6$RowName]

postn_ko_6[postn_ko_6>0.6]=0.6
p2=qcorrplot(postn_ko_6, 
             type = "upper",
             diag = FALSE) +
  geom_tile(aes(fill = r), 
            width = 0.98, 
            height = 0.98,
            color = "white",
            linewidth = 0.2) +
  scale_fill_gradientn(
    colours = cols,
    oob     = scales::squish
  )

ntc_6[ntc_6>0.6]=0.6
p3=qcorrplot(ntc_6, 
             type = "lower", 
             diag = FALSE) +
  geom_tile(aes(fill = r), 
            width = 0.98,
            height = 0.98,
            color = "white",
            linewidth = 0.2)+
  scale_fill_gradientn(
    colours = cols[8:14],
    oob     = scales::squish
  )

## PSOTN CRISPR_KO  

tfc_6_7_ntc=NTC_TF_correlation[c7$RowName,c6$RowName]
tfc_6_7_ntc[tfc_6_7_ntc>0.6]=0.6
tfc_6_7_ntc[tfc_6_7_ntc<(-0.6)]=-0.6
bk <- c(seq(-0.6,-0.2,by=0.01),seq(0.21,0.6,by=0.01))
pheatmap::pheatmap(tfc_6_7_ntc, 
                   cluster_rows = F,
                   cluster_cols = F,
                   treeheight_row = 0,
                   scale = 'none',
                   color = c(colorRampPalette(c("#084594","white"))(length(bk)/2),
                             "white",colorRampPalette(c("white","#99000d"))(length(bk)/2)),
                   breaks = bk,
                   cellwidth = 14, 
                   cellheight = 14, 
                   fontsize = 12,
                   border_color = '#ffffff',
                   filename = file.path(output_dir, "Figure4e_NTC_TFc6_TFc7_heatmap.pdf")
                  )

tfc_6_7_postn_ko=postn_tf_cors[c6$RowName,c7$RowName]
tfc_6_7_postn_ko[tfc_6_7_postn_ko>0.6]=0.6
tfc_6_7_postn_ko[tfc_6_7_postn_ko<(-0.6)]=-0.6

pheatmap::pheatmap(tfc_6_7_postn_ko, 
                   cluster_rows = F,
                   cluster_cols = F,
                   treeheight_row = 0,
                   scale = 'none',
                   color = c(colorRampPalette(c("#084594","white"))(length(bk)/2),
                             "white",colorRampPalette(c("white","#99000d"))(length(bk)/2)),
                   breaks = bk,
                   cellwidth = 14, 
                   cellheight = 14,
                   fontsize = 12,
                   border_color = '#ffffff',
                   filename = file.path(output_dir, "Figure4e_POSTN_KO_TFc6_TFc7_heatmap.pdf")
)





glis1_mecom_NTC<-read.csv("GLIS1_MECOM_NTC.csv",
                          header = T,
                          row.names = 1,
                          check.names = F)
glis1_mecom_KO<-read.csv("GLIS1_MECOM_PODTN_KO.csv",
                         header = T,
                         row.names = 1,
                         check.names = F)

p1<-ggscatter(glis1_mecom_NTC, 
              x = "GLIS1",
              y = "MECOM",
              alpha=0.5,
              size = 0.005,
              color = "#A1A1A1",
              add = "reg.line",
              conf.int = FALSE,
              add.params = list(color = "black", 
                                fill = "lightgray",
                                size=0.4),        
                legend='right')+stat_cor()+ggtitle("NTC")

p2<-ggscatter(glis1_mecom_KO,
              x = "GLIS1",
              y = "MECOM",
              alpha=0.5,
              size = 1.5,
              color = "#0A526D",
              add = "reg.line",
              conf.int = FALSE,
              add.params = list(color = "black",
                                fill = "lightgray",
                                size=0.4),      
                  legend='right')+stat_cor()+ggtitle("POSTN_CRISPR_KO")
p=p1+p2
ggsave(p, file = file.path(output_dir, "Figure4f_GLIS1_MECOM.pdf"))




znf117_mecom_NTC<-read.csv("ZNF117_MECOM_NTC.csv",
                          header = T,
                          row.names = 1,
                          check.names = F)
znf117_mecom_KO<-read.csv("ZNF117_MECOM_POSTN_KO.csv",
                         header = T,
                         row.names = 1,
                         check.names = F)

p1<-ggscatter(znf117_mecom_NTC, 
              x = "ZNF117",
              y = "MECOM",
              alpha=0.5,
              size = 0.005,
              color = "#A1A1A1",
              add = "reg.line",
              conf.int = FALSE,
              add.params = list(color = "black", 
                                fill = "lightgray",
                                size=0.4),        
              legend='right')+stat_cor()+ggtitle("NTC")

p2<-ggscatter(znf117_mecom_KO,
              x = "ZNF117",
              y = "MECOM",
              alpha=0.5,
              size = 1.5,
              color = "#0A526D",
              add = "reg.line",
              conf.int = FALSE,
              add.params = list(color = "black",
                                fill = "lightgray",
                                size=0.4),      
              legend='right')+stat_cor()+ggtitle("POSTN_CRISPR_KO")
p=p1+p2
ggsave(p, file = file.path(output_dir, "Figure4g_ZNF117_MECOM.pdf"))





mybl1_nr2c2_NTC<-read.csv("MYBL1_NR2C2_NTC.csv",
                           header = T,
                           row.names = 1,
                           check.names = F)
mybl1_nr2c2_KO<-read.csv("MYBL1_NR2C2_POSTN_KO.csv",
                          header = T,
                          row.names = 1,
                          check.names = F)

p1<-ggscatter(mybl1_nr2c2_NTC, 
              x = "MYBL1",
              y = "NR2C2",
              alpha=0.5,
              size = 0.005,
              color = "#A1A1A1",
              add = "reg.line",
              conf.int = FALSE,
              add.params = list(color = "black", 
                                fill = "lightgray",
                                size=0.4),        
              legend='right')+stat_cor()+ggtitle("NTC")

p2<-ggscatter(mybl1_nr2c2_KO,
              x = "MYBL1",
              y = "NR2C2",
              alpha=0.5,
              size = 1.5,
              color = "#0A526D",
              add = "reg.line",
              conf.int = FALSE,
              add.params = list(color = "black",
                                fill = "lightgray",
                                size=0.4),      
              legend='right')+stat_cor()+ggtitle("POSTN_CRISPR_KO")
p=p1+p2
ggsave(p, file = file.path(output_dir, "Figure4h_MYBL1_NR2C2.pdf"))


mybl2_nr2c2_NTC<-read.csv("MYBL2_NR2C2_NTC.csv",
                          header = T,
                          row.names = 1,
                          check.names = F)
mybl2_nr2c2_KO<-read.csv("MYBL2_NR2C2_POSTN_KO.csv",
                         header = T,
                         row.names = 1,
                         check.names = F)

p1<-ggscatter(mybl2_nr2c2_NTC, 
              x = "MYBL2",
              y = "NR2C2",
              alpha=0.5,
              size = 0.005,
              color = "#A1A1A1",
              add = "reg.line",
              conf.int = FALSE,
              add.params = list(color = "black", 
                                fill = "lightgray",
                                size=0.4),        
              legend='right')+stat_cor()+ggtitle("NTC")

p2<-ggscatter(mybl2_nr2c2_KO,
              x = "MYBL2",
              y = "NR2C2",
              alpha=0.5,
              size = 1.5,
              color = "#0A526D",
              add = "reg.line",
              conf.int = FALSE,
              add.params = list(color = "black",
                                fill = "lightgray",
                                size=0.4),      
              legend='right')+stat_cor()+ggtitle("POSTN_CRISPR_KO")
p=p1+p2
ggsave(p, file = file.path(output_dir, "Figure4i_MYBL2_NR2C2.pdf"))

# Cytotrace scores
ntc<-readRDS("postn_proliferation.rds")
plotdat<-data.frame(gene_score=ntc$CytoTRACE2_Score,groups=ntc$group)
plotdat<- plotdat %>% mutate(group=forcats::fct_reorder(groups,gene_score, .fun = 'median',.desc = F) )
p=ggviolin(plotdat,
           x = "group", 
           y = "gene_score",
           color='group',
           fill = "group",
           legend="right",
           outlier=FALSE,
           add = "boxplot",
           palette = c("#eb7e60","grey"),
           add.params = list(fill = "white",outlier.shape = NA), 
           error.plot = "errorbar",border="white")+stat_compare_means()
p

# Proliferation scores
plotdat<-data.frame(gene_score=ntc$proliferationScore,groups=ntc$group)
plotdat<- plotdat %>% mutate(group=forcats::fct_reorder(groups,gene_score, .fun = 'median',.desc = F) )
p=ggviolin(plotdat,
           x = "group", 
           y = "gene_score",
           color='group',
           fill = "group",
           legend="right",
           outlier=FALSE,
           add = "boxplot",
           palette = c("#eb7e60","grey"),
           add.params = list(fill = "white",outlier.shape = NA), 
           error.plot = "errorbar",border="white")+stat_compare_means()
p

# Matrix remodling scores
plotdat<-data.frame(gene_score=ntc$matrix_remodling1,groups=ntc$group)
plotdat<- plotdat %>% mutate(group=forcats::fct_reorder(groups,gene_score, .fun = 'median',.desc = F) )
p=ggviolin(plotdat,
           x = "group", 
           y = "gene_score",
           color='group',
           fill = "group",
           legend="right",
           outlier=FALSE,
           add = "boxplot",
           palette = c("#eb7e60","grey"),
           add.params = list(fill = "white",outlier.shape = NA), 
           error.plot = "errorbar",border="white")+stat_compare_means()
p

## GATA6 - ZNF608
gata6_expr<-read.csv("GATA6_expression.csv",header = T,row.names = 1,check.names = F)
znf608_expr<-read.csv("ZNF608_expression.csv",header = T,row.names = 1,check.names = F)
ntc_expr<-read.csv("GATA6_ZNF608_NTC.csv",header = T,row.names = 1,check.names = F)
actg2_ko<-read.csv("GATA6_ZNF608_ACTG2_KO.csv",header = T,row.names = 1,check.names = F)

p1=ggviolin(gata6_expr,
              x="group",
              y="expression",
              fill ="group",
              color = "group",
              palette = c("grey","#114D72"),
              col="TFs",width = 0.5
              )+stat_compare_means()+ggtitle("GATA6")
p2=ggviolin(znf608_expr,
              x="group",
              y="expression",
              fill ="group",
              color = "group",
              palette = c("grey","#114D72"),
              col="TFs",width = 0.5
              )+stat_compare_means()+ggtitle("ZNF608")

p3<-ggscatter(ntc_expr,
              x = "GATA6",
              y = "ZNF608",
              alpha=0.8,
              size = 0.5,
              color = "#A1A1A1",
              add = "reg.line",
              conf.int = TRUE,
              add.params = list(color = "black", 
                                fill = "lightgray"), 
              legend='right')+stat_cor(method = "pearson")+ggtitle("NTC")

p4<-ggscatter(actg2_ko,
              x = "GATA6",
              y = "ZNF608",
              alpha=0.8,
              size = 0.5,
              color = "#114D72",
              add = "reg.line",
              conf.int = TRUE,
              add.params = list(color = "black", 
                                fill = "lightgray"),        
              legend='right')+stat_cor()+ggtitle("ACTG2_KO")
pp=ggarrange(p1,p2,p3, p4,ncol = 4, nrow = 1)
ggsave(pp, file = file.path(output_dir, "Figure4l_m_GATA6_ZNF608_ACTG2_KO.pdf"))


## GATA6 - ZNF608
rarb_expr<-read.csv("RARB_expression.csv",header = T,row.names = 1,check.names = F)
twist2_expr<-read.csv("TWIST2_expression.csv",header = T,row.names = 1,check.names = F)
ntc_expr<-read.csv("RARB_TWIST2_NTC.csv",header = T,row.names = 1,check.names = F)
krtap3_1_ko<-read.csv("RARB_TWIST2_KRTAP3-1_KO.csv",header = T,row.names = 1,check.names = F)

p1=ggviolin(rarb_expr,
            x="group",
            y="expression",
            fill ="group",
            color = "group",
            palette = c("grey","#114D72"),
            col="TFs",
            width = 0.5
)+stat_compare_means()+ggtitle("RARB")
p2=ggviolin(twist2_expr,
            x="group",
            y="expression",
            fill ="group",
            color = "group",
            palette = c("grey","#114D72"),
            col="TFs",
            width = 0.5
)+stat_compare_means()+ggtitle("TWIST2")

p3<-ggscatter(ntc_expr,
              x = "RARB",
              y = "TWIST2",
              alpha=0.8,
              size = 0.5,
              color = "#A1A1A1",
              add = "reg.line",
              conf.int = TRUE,
              add.params = list(color = "black", 
                                fill = "lightgray"), 
              legend='right')+stat_cor(method = "pearson")+ggtitle("NTC")

p4<-ggscatter(krtap3_1_ko,
              x = "RARB",
              y = "TWIST2",
              alpha=0.8,
              size = 0.5,
              color = "#114D72",
              add = "reg.line",
              conf.int = TRUE,
              add.params = list(color = "black", 
                                fill = "lightgray"),        
              legend='right')+stat_cor()+ggtitle("ACTG2_KO")
pp=ggarrange(p1,p2,p3, p4,ncol = 4, nrow = 1)
ggsave(pp, file = file.path(output_dir, "Figure4n_o_RARB_TWIST2_KRTAP3-1_KO.pdf"))




