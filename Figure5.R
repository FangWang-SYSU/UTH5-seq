library(Seurat)
library(Seurat)
library(ggpubr)
library(Nebulosa)


ntcs<-readRDS("fibroblast.rds")

umapColor1=c("#7F3C8D","#11A579","#3969AC","#E73F74","#80BA5A",
             "#E68310", "#008695","#CF1C90", "#f97b72", "#E7ABFD"
             ,"#9a221c" ,"#FED8A3" ,"#fcb93e")

p1=DimPlot(ntcs,
           group.by = "seurat_clusters",
           cols = umapColor1,
           alpha=0.8,
           pt.size = 0.4)
ggsave(p1,file="NTC_umap.pdf")


plotdat<-data.frame(cytotrace_score=ntcs$CytoTRACE2_Score,
                    groups=ntcs$cellstates)



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

ggsave("cluster_cytotrace_score.pdf",p,)

## marker gene expression
sc_dat.markers<-read.csv("/Volumes/zz_20t/Project/CRISPR_CAF/For_public/Ref/Fig2/marker_gene_functions_df_unique.csv",header = T,row.names = 1,check.names = F)

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

### PRRX1 and ACTA2 density
p <- plot_density(controlcell, 
                  c("PRRX1","ACTA2"),
                  slot = "data",
                  adjust = 0.01,
                  pal = "plasma",
                  size = 0.1)

ggsave(p,file="PRRX1_ACTA2_density_umap.pdf")

## marker genes bubble plot
featureplot<-read.csv("feature_marker_plot_buble.csv",header = T,row.names = 1,check.names = F)
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
ggsave(pp,file="feature_marker_plot_buble1111.pdf")







