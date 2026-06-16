library(Seurat)
library(Seurat)
library(ggpubr)
library(Nebulosa)

replot_data <- read.csv("Figure5a_source_data.csv", check.names = FALSE)
cluster_values <- unique(replot_data[["Seurat cluster"]])
cluster_numbers <- suppressWarnings(as.numeric(cluster_values))
if (all(!is.na(cluster_numbers))) {
  cluster_levels <- as.character(sort(cluster_numbers))
} else {
  cluster_levels <- sort(as.character(cluster_values))
}
replot_data[["Seurat cluster"]] <- factor(
  replot_data[["Seurat cluster"]],
  levels = cluster_levels
)

figure5a_plot <- ggplot(
  replot_data,
  aes(x = `UMAP 1`, y = `UMAP 2`, color = `Seurat cluster`)
) +
  geom_point(size = 0.25, stroke = 0, alpha = 1) +
  coord_equal() +
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
  plot_file,
  figure5a_plot,
  width = 4.2,
  height = 3.8,
  units = "in",
  useDingbats = FALSE
)




## marker gene expression
sc_dat.markers<-read.csv("Figure5b_scouce_data.csv",header = T,row.names = 1,check.names = F)

#write.csv(df_unique,file="df_unique_merge5611.csv")

subdat<-subset(ntcs,features = sc_dat.markers$gene)
Idents(subdat)<-subdat$cellstates
#subdat=NormalizeData(subdat)

subdat[["RNA"]]$data <- as(object = subdat[["RNA"]]$data, Class = "dgCMatrix")
avemat<-AverageExpression(subdat,layer = 'data')
avemat=data.frame(avemat)
colnames(avemat)<-gsub("RNA.","S",colnames(avemat),fixed = T)

rowcol=as.list(umapColor1)
names(rowcol)=paste0("S",0:12)
rowanno<-rowAnnotation(Clusters=sc_dat.markers$cluster,col=rowcol)

avemat=t(avemat)
avemat1=scale(avemat,center = T)
markerexpression=data.frame(t(avemat1))
markerexpression=markerexpression[sc_dat.markers$gene,]
color =colorRamp2(c(-2, -0.1,0.1, 2), c("#5ab4ac", "white","white", "#b2182b"))
sc_dat.markers$cluster=as.character(sc_dat.markers$cluster)
sc_dat.markers$cluster=factor(sc_dat.markers$cluster,levels = c("S0","S1","S2","S3","S8","S5","S6","S7","S10","S11","S4","S9","S12"))
splitrow=sc_dat.markers$cluster
pdf("marker_gene_expression.pdf",width =40, height =20)
ht=ComplexHeatmap:: Heatmap(
  markerexpression,
  cluster_rows = F,
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
  column_names_rot=0
)
draw(ht,
     heatmap_legend_side="right", 
     annotation_legend_side="right",
     legend_grouping = "original")
dev.off() 



plotdat<-read.csv("Figure5c_source_data.csv", check.names = FALSE)
plotdat<- plotdat %>% mutate(group=forcats::fct_reorder(groups,cytotrace_score, .fun = 'median',.desc = T) )

p=ggviolin(plotdat,
           x = "group",
           y = "cytotrace_score",
           color='group', fill = "group",
           legend="right",
           add = "boxplot",
           palette = umapColor1,
           add.params = list(fill = "white"),
           error.plot = "errorbar",border="white")
p

ggsave(p,file="cluster_cytotrace_score.pdf")


igure5d_source_data <- read.csv("Figure5d_source_data.csv", header = TRUE, check.names = FALSE)
p1 <- ggplot(Figure5d_source_data, aes(`UMAP 1`, `UMAP 2`, color = `PRRX1 density`)) +
  geom_point(size = 0.1) +
  scale_color_viridis_c(option = "plasma", name = "Density") +
  ggtitle("PRRX1") +
  theme_void()
p2 <- ggplot(Figure5d_source_data, aes(`UMAP 1`, `UMAP 2`, color = `ACTA2 density`)) +
  geom_point(size = 0.1) +
  scale_color_viridis_c(option = "plasma", name = "Density") +
  ggtitle("ACTA2") +
  theme_void()
p <- p1 + p2
ggsave(p,file="Figure5d_figuree.pdf")

## marker genes bubble plot
featureplot<-read.csv("Figure5f_source_data.csv",header = T,row.names = 1,check.names = F)
featureplot$anno=factor(featureplot$anno,levels = c(unique(featureplot$anno)))

featureplot$features.plot=factor(featureplot$features.plot,
                                 levels=rev(c(unique(featureplot$features.plot))))
featureplot$features=featureplot$`features.plot`
featureplot$cluster=paste0("S",featureplot$cluster)
featureplot$cluster=factor(featureplot$cluster,
                           levels=c("S0","S1","S2","S3","S8","S5",
                                    "S6","S7","S10","S11","S4","S9","S12"))

cols=rev(c('#b2182b','#ef8a62','#fddbc7','#e0e0e0','#999999','#4d4d4d'))

pp=ggscatter(featureplot,
             x="cluster",
             y="features",
             color = 'avg.exp.scaled',
             size='pct_scaled',
             legend='right')+ 
  scale_color_gradientn(
  colours = cols,    
  limits  = c(-2, 2), 
  oob     = squish   
)
pp
ggsave(pp,file="feature_marker_plot_buble.pdf")

prop<-read.csv("Figure5k.csv",header=T)
p=ggplot(prop, aes( x = Group,y=Proportion,fill = clusters, stratum = clusters, 
                          alluvium = clusters))+ geom_col(width = 0.5,  color = 'white', size = 0.5) + 
  geom_flow(width = 0.5, alpha = 0.3, knot.pos = 0.2, color = 'white', size = 0.5)+ 
  coord_flip()+ 
  scale_fill_manual(values = umapColor2)+ 
  scale_y_continuous(expand = c(0,0),name="", 
                     label=c("0%","25%","50%","75%","100%"))+
  scale_x_discrete(expand = c(0,0),name="")+ 
  theme(panel.background = element_blank(),
        panel.grid = element_blank(), 
        axis.line = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text = element_text(color="black",size=10),
        axis.ticks.length.x = unit(0.1,"cm"),
        plot.margin = margin(10, 10, 10, 10)) 

ggsave(p,file="Figure5k.pdf",width=6,height = 4)

cols=c('#ffffcc','#ffeda0','#fed976','#feb24c','#fd8d3c','#fc4e2a','#e31a1c','#bd0026','#800026')
vln.df<-read.csv("Figure5g.csv",header=T,row.names=1,check.names=F)
vln.df$group=paste0("S",vln.df$group)
 p<-vln.df%>%ggplot(aes(group,exp))+geom_violin(aes(fill=mean),scale = "width",linewidth=0)+geom_boxplot(width=0.3,outlier.shape = NA,color="white")+
  # facet_wrap(vln.df$gene~.,scales = "free_y",nrow=10)+
   #  labs(title = paste0(genes," expression"))+
   scale_fill_gradientn(colors =cols)+
   #scale_fill_viridis_c(option = "B",
   scale_x_discrete("")+
   theme_bw()+theme(
     
     panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
     axis.text.x.bottom = element_text(angle = 0,hjust = 1,vjust = 1,size=8)
   )+theme(strip.text.y = element_text(angle = 0)) +  theme(
     strip.text = element_text(face = "bold", size = rel(0.5)),
     strip.background = element_rect(fill = "white", colour = "black", size = 1),strip.placement = "outside")+
   force_panelsizes(rows = unit(1.2, "in"),cols = unit(4, "in"))
 
 




