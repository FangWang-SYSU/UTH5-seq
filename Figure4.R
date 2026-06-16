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
ggsave(p2,file="Fig2b.pdf")
ggsave(p3,file="Fig2b.pdf")

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
                   border_color = '#ffffff'
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
                   border_color = '#ffffff'
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
ggsave(p,file="Fig2c.pdf")




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
ggsave(p,file="Fig2d.pdf")





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
ggsave(p,file="Fig2e.pdf")


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
ggsave(p,file="Fig2f.pdf")



e2f7_mecom_NTC<-read.csv("E2F7_MECOM_NTC.csv",
                          header = T,
                          row.names = 1,
                          check.names = F)
e2f7_mecom_KO<-read.csv("E2F7_MECOM_PODTN_KO.csv",
                         header = T,
                         row.names = 1,
                         check.names = F)

p1<-ggscatter(e2f7_mecom_NTC, 
              x = "E2F7",
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

p2<-ggscatter(e2f7_mecom_KO,
              x = "E2F7", 
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
ggsave(p,file="Fig2f.pdf")

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
ggsave(p,file="Cytotrace_score.pdf")

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
ggsave(p,file="Proliferation_score.pdf")

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
ggsave(p,file="Proliferation_score.pdf")

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
ggsave(pp,file="GATA6_ZNF608_ACTG2_KO.pdf")


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
ggsave(pp,file="RARB_TWIST2_KRTAP3-1_KO.pdf")












