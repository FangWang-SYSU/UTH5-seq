library(data.table)
library(reshape2)
library(ggplot2)
library(parallel)
library(dplyr)
library(ggpubr)
library(fmsb)

seu<-readRDS("firbo.rds")
modules<-paste0("Module",1:8)
for(ms in modules){
  Idents(seu)<- seu@meta.data[,ms]
  CellstoHighlight <- WhichCells(seu, idents = ms)
  plotpl1 <- UMAPPlot(seu,combine = FALSE, 
                      cells.highlight=CellstoHighlight,
                      cols.highlight='#FC7216')
  scale.col <- cet_pal(18,name='fire')[3:18]
  highlight_data <- seu@reductions$umap@cell.embeddings[CellstoHighlight, , drop = FALSE]
  highlight_df <- as.data.frame(highlight_data)
  colnames(highlight_df) <- c("UMAP_1", "UMAP_2")
  
  pp <- plotpl1[[1]]
  pp <- pp + stat_density_2d(
    data = highlight_df,  # only show highlight cell density
    aes(x = UMAP_1, 
        y = UMAP_2,
        fill = after_stat(level)),
    geom = "density_2d_filled",
    colour = "ivory",
    n = 150,
    h = c(1, 1)
  ) + scale_fill_gradientn(colours = scale.col) +
    ggtitle(mm)+theme_few()
  
  
  ggsave(pp,file=paste0(ms,"_density.pdf"))

}

# Heatmap of module enrichment 

moduleEnrich<-read.csv("Module_cluster_odds_ratio.csv",
                       header = T,
                       row.names = 1,
                       check.names = F)
moduleEnrich=t(moduleEnrich)
moduleEnrich=moduleEnrich[c("S5","S6","S10","S12","S11","S0","S1","S3","S7","S8","S9","S2","S4"),
                          c("Module1","Module2","Module4","Module3","Module5","Module6","Module7","Module8")]

moduleEnrich[moduleEnrich>1.5]=1.5
moduleEnrich[moduleEnrich<(-1.5)]=-1.5
pheatmap::pheatmap(moduleEnrich,
                   scale = "none",
                   cluster_rows = F,
                   cluster_cols = F,
                   border_color = "white",
                   color = colorRampPalette(c("#33ABC1","white", "#B11927"))(20)
                   )

## Top10 perturbation bubble plot 
plotdat<-read.csv("Top_10_cluster_prop_and ration_in pertub_module_fisher_plot_buble.csv",
                  header = T,
                  row.names = 1,
                  check.names = F)
plotdat$Module=gsub("M_","Module",plotdat$Module)

plotdat$Module=factor(plotdat$Module,levels=rev(c("Module1","Module2","Module4","Module3","Module5","Module6","Module7","Module8")))
plotdat=df_long[order(plotdat$Module),]
plotdat$Perturbation=factor(plotdat$Perturbation,levels = unique(plotdat$Perturbation))
plotdat$cluster=paste0("S",plotdat$cluster)
plotdat$cluster=factor(plotdat$cluster,levels = c("S5","S6","S10","S12","S11","S0","S1","S3","S7","S8","S9","S2","S4"))

cols=c("#313772","#2c4ca0","#326db6","#478ecc","#75b5dc","white","#c44438","#b7282e")
pp=ggscatter(plotdat,
             x="cluster",
             y="Perturbation",
             color = 'log2ratio',
             size='log10pval',
             legend='right')+ 
  scale_color_gradientn(
  colours = cols,     
  limits  = c(-6, 2), 
  oob     = squish  
)

ggsave(pp,file="top10_perturbation_plot_buble.pdf",width = 6,height = 17)


### ACTG2 and MEGF8 function
enrichDat<-read.table("cluster_0_ACTG2_MFGE8.txt",header = T,check.names = F,sep="\t")

enrichDat$Description=factor(enrichDat$Description,levels = unique(enrichDat$Description))
p=ggscatter(enrichDat,
            x="groupby",
            y="Description",
            color = "groups",
            size ="log10qvalue",
            legend="right",
            palette = c("#45174E","#D86B33"))+
  rotate_x_text(45)
ggsave(p,file="Fig5d")

## fig5e
plotdat<-read.csv("MFGE8_compare_up_GOBP_plot_compare.csv",header = T,row.names = 1,check.names = F)

pdf("Fig5e.pdf",width =8, height =8)
radarchart(plotdat,
           pty = c(16,16,32),
           axistype = 1,
           pcol = c("#008A89", "#D92927","#3272AF"),
           plwd = c(3,3,3),
           plty = 1, 
           cglcol = "grey60", 
           cglty = 1, 
           cglwd = 1, 
           axislabcol = "grey60", 
           vlcex = 0.8,
           vlabels = colnames(colnames(plotdat)),
           caxislabels = c(0, 10, 20, 30, 40),
           calcex=0.8 
)

legend(x = "bottomright", 
       legend = c("S0","S10","S3"),
       horiz = F,
       bty = "n",
       pch = 15 ,
       col = c("#008A89", "#D92927","#3272AF"),
       text.col = "black",
       cex = 1, pt.cex = 1.5)

dev.off()

## Fig5f
plotdat<-read.csv("MFGE8_compare_Down_GOBP_plot_compare.csv",header = T,row.names = 1,check.names = F)
pdf("Fig5f.pdf",width =8, height =8)
radarchart(plotdat,
           pty = c(16,16,32),
           axistype = 1,
           pcol = c("#74C9BC", "#F79C7B","#8FB4DC"),
           plwd = c(3,3,3),
           plty = 1, 
           cglcol = "grey60", 
           cglty = 1, 
           cglwd = 1, 
           axislabcol = "grey60", 
           vlcex = 0.8,
           vlabels = colnames(colnames(plotdat)),
           caxislabels = c(0, 10, 20, 30, 40),
           calcex=0.8 
)

legend(x = "bottomright", 
       legend = c("S0","S10","S3"),
       horiz = F,
       bty = "n",
       pch = 15 ,
       col = c("#74C9BC", "#F79C7B","#8FB4DC"),
       text.col = "black",
       cex = 1, pt.cex = 1.5)

dev.off()


## Fig5g
plotdat<-read.csv("EPB41L1_compare_Up_GOBP_plot_compare.csv",header = T,row.names = 1,check.names = F)
pdf("Fig5g.pdf",width =8, height =8)
radarchart(plotdat,
           pty = c(16,16,32),
           axistype = 1,
           pcol = c("#008A89", "#D92927"),
           plwd = c(3,3,3),
           plty = 1, 
           cglcol = "grey60", 
           cglty = 1, 
           cglwd = 1, 
           axislabcol = "grey60", 
           vlcex = 0.8,
           vlabels = colnames(colnames(plotdat)),
           caxislabels = c(0, 10, 20, 30, 40),
           calcex=0.8 
)

legend(x = "bottomright", 
       legend = c("S0","S10"),
       horiz = F,
       bty = "n",
       pch = 15 ,
       col = c("#008A89", "#D92927"),
       text.col = "black",
       cex = 1, pt.cex = 1.5)

dev.off()


## Fig5h
plotdat<-read.csv("EPB41L1_compare_Down_GOBP_plot_compare.csv",header = T,row.names = 1,check.names = F)
pdf("Fig5h.pdf",width =8, height =8)
radarchart(plotdat,
           pty = c(16,16,32),
           axistype = 1,
           pcol = c("#74C9BC", "#F79C7B"),
           plwd = c(3,3,3),
           plty = 1, 
           cglcol = "grey60", 
           cglty = 1, 
           cglwd = 1, 
           axislabcol = "grey60", 
           vlcex = 0.8,
           vlabels = colnames(colnames(plotdat)),
           caxislabels = c(0, 10, 20, 30, 40),
           calcex=0.8 
)

legend(x = "bottomright", 
       legend = c("S0","S10"),
       horiz = F,
       bty = "n",
       pch = 15 ,
       col = c("#74C9BC", "#F79C7B"),
       text.col = "black",
       cex = 1, pt.cex = 1.5)

dev.off()

## Fig5i
plotdat<-read.csv("GPR1_compare_Up_GOBP_plot_compare.csv",header = T,row.names = 1,check.names = F)
pdf("Fig5i.pdf",width =8, height =8)
radarchart(plotdat,
           pty = c(16,16,32),
           axistype = 1,
           pcol = c("#008A89", "#764FA0"),
           plwd = c(3,3,3),
           plty = 1, 
           cglcol = "grey60", 
           cglty = 1, 
           cglwd = 1, 
           axislabcol = "grey60", 
           vlcex = 0.8,
           vlabels = colnames(colnames(plotdat)),
           caxislabels = c(0, 10, 20, 30, 40),
           calcex=0.8 
)

legend(x = "bottomright", 
       legend = c("S0","S1"),
       horiz = F,
       bty = "n",
       pch = 15 ,
       col = c("#008A89", "#764FA0"),
       text.col = "black",
       cex = 1, pt.cex = 1.5)

dev.off()




## Fig5j
plotdat<-read.csv("GPR1_compare_Down_GOBP_plot_compare.csv",header = T,row.names = 1,check.names = F)
pdf("Fig5j.pdf",width =8, height =8)
radarchart(plotdat,
           pty = c(16,16,32),
           axistype = 1,
           pcol = c("#74C9BC", "#AB99C8"),
           plwd = c(3,3,3),
           plty = 1, 
           cglcol = "grey60", 
           cglty = 1, 
           cglwd = 1, 
           axislabcol = "grey60", 
           vlcex = 0.8,
           vlabels = colnames(colnames(plotdat)),
           caxislabels = c(0, 10, 20, 30, 40),
           calcex=0.8 
)

legend(x = "bottomright", 
       legend = c("S0","S1"),
       horiz = F,
       bty = "n",
       pch = 15 ,
       col = c("#74C9BC", "#AB99C8"),
       text.col = "black",
       cex = 1, pt.cex = 1.5)

dev.off()
