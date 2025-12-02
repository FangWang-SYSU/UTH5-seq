library(ggplot2)
library(data.table)
library(reshape2)
library(ggplot2)
library(parallel)
library(Seurat)
library(reticulate)
library(anndata)
library(SeuratDisk) 
library(reticulate)

seu=readRDS("fibro.rds")

## Response clusters
umapColor1=c("#7F3C8D","#11A579","#3969AC","#E73F74","#80BA5A",
             "#E68310", "#008695","#CF1C90", "#f97b72", "#E7ABFD"
             ,"#9a221c" ,"#FED8A3" ,"#fcb93e" ,"#F9BEBB")

umapColor2 <- c('#53A85F','#E5D2DD','#E59CC4','#D6E7A3','#F1BB72','#F3B1A0')
perturbtions =  c("ACTG2","MFGE8","GPR1","PTX4","EGFL6","IL1RAPL1","OLFM3")
for (pt in perturbtions){
  data <- read_h5ad(paste0(,pt,"_leiden_clusters.h5ad"))
  subdat<-subset(seu,cellnewname%in%data$obs[,"cellnewname"])
  umapdata<-data.frame(data$obsm[["X_umap"]])
  rownames(umapdata)<-data$obs[,"cellnewname"]
  colnames(umapdata)<-c("UMAP_1","UMAP_2")
  umapdata=umapdata[subdat$cellnewname,]
  umapdata=as.matrix(umapdata[,c("UMAP_1","UMAP_2")])
  umap_reduction <- CreateDimReducObject(embeddings = umapdata, key = "UMAP_", assay = DefaultAssay(subdat))
  subdat@reductions$umap <- umap_reduction
  subdat$response_cluster=data$obs[,"leiden"]
  p1=DimPlot(subdat,group.by = "seurat_clusters",cols = umapColor2,alpha=0.8,pt.size = 0.4)
p2=DimPlot(subdat,group.by = "response_cluster",cols =cols4,alpha=0.8,pt.size = 0.4)
pp=p1/p2
ggsave(pp,file=paste0(pt,"_scdat_mRNA.pdf"))

}


synergyScore=read.csv("synergy_score.csv",header = T,check.names = F)
synergyScore$clusters=factor(synergyScore$clusters,levels = paste0("S",0:12))
p<-synergyScore%>%ggplot(aes(clusters,synergyScore))+
  geom_violin(aes(fill=clusters),
              scale = "width",
              linewidth=0)+
  facet_wrap(~groups,
             scales = "fixed",
             nrow=3)+
  scale_fill_manual(values = umapColor1) +
  theme_classic2()

ggsave(p,file="synergy_score.pdf")



seu=readRDS("synergy_score.rds")
genes<-c("IGFBP5","COL8A1","USP53","CHSY3",
         "CITED2","CTSB","PRR16","SOX5",
         "NRXN3","PLCB4","CHRM2")
mat=as.matrix(GetAssayData(seu,layer = "data"))
mat=data.frame(t(mat))
mat=mat[,genes]
mat$groups=seu$groups
df<-reshape2::melt(mat)
colnames(df)<-c("group","gene","exp")
df$group_gene=paste(df$group,df$gene,sep = "_")
#stat.df=as.data.frame(vln.df%>%group_by(celltype,gene)%>%summarise(mean=mean(expm1(exp))))
expr.df=as.data.frame(df%>%group_by(group,gene)%>%summarise(mean=mean((exp))))
expr.df$group_gene=paste(expr.df$group,expr.df$gene,sep = "_")
expr.df=expr.df[,c("mean","group_gene")]
df=inner_join(df,expr.df,by="group_gene")
df$group<-reorder(df$group,-df$mean)
df$group=factor(df$group,levels = 
                  rev(c("NTC","PTX4","IL1RAPL1",
                        "OLFM3","GPR1","OLFM3_PTX4",
                        "IL1RAPL1_PTX4","GPR1_PTX4")))

cols=c('#fcc5c0','#f768a1','#dd3497','#7a0177')
df$gene=factor(df$gene,levels = genes)

p<-df%>%ggplot(aes(gene,exp))+
  geom_violin(aes(fill=mean),
              scale = "width",linewidth=0)+
  facet_wrap(~group,scales = "free_y",nrow=10)+
  scale_fill_gradientn(colors =cols)+
  theme_classic2()
p
ggsave(p,file="shared_synergy_genes.pdf")




genes<-c("GPR176","KCNQ1OT1","TBC1D19","C5orf30",
         "LIMCH1","CALM2","MYH9","MSN","OXR1",
         "HEG1","EIF4A2","SDC2","SKA2","ENO1")
mat=as.matrix(GetAssayData(seu,layer = "data"))
mat=data.frame(t(mat))
mat=mat[,genes]
mat$groups=seu$groups
df<-reshape2::melt(mat)
colnames(df)<-c("group","gene","exp")
df$group_gene=paste(df$group,df$gene,sep = "_")
#stat.df=as.data.frame(vln.df%>%group_by(celltype,gene)%>%summarise(mean=mean(expm1(exp))))
expr.df=as.data.frame(df%>%group_by(group,gene)%>%summarise(mean=mean((exp))))
expr.df$group_gene=paste(expr.df$group,expr.df$gene,sep = "_")
expr.df=expr.df[,c("mean","group_gene")]
df=inner_join(df,expr.df,by="group_gene")
df$group<-reorder(df$group,-df$mean)
df$group=factor(df$group,levels = 
                  rev(c("NTC","PTX4","IL1RAPL1",
                        "OLFM3","GPR1","OLFM3_PTX4",
                        "IL1RAPL1_PTX4","GPR1_PTX4")))

cols=c('#fcc5c0','#f768a1','#dd3497','#7a0177')
df$gene=factor(df$gene,levels = genes)

p<-df%>%ggplot(aes(gene,exp))+
  geom_violin(aes(fill=mean),
              scale = "width",linewidth=0)+
  facet_wrap(~group,scales = "free_y",nrow=10)+
  scale_fill_gradientn(colors =cols)+
  theme_classic2()
p

ggsave(p,file="Specific_synergy_genes.pdf")



plotdat<-read.table("GO_BP_top10.txt",header = T,sep="\t",check.names = F)
plotdat$Description=factor(plotdat$Description,levels = unique(plotdat$Description))
pp=ggscatter(plotdat,
             x="groups",
             y="Description",
             color = 'qvalue',
             size='Count',
             legend='right')+
  scale_color_gradientn(colors =cols) +
  scale_size_continuous(name = "Count",  breaks = c(3, 5),
                        range = c(4, 6))+rotate_x_text(45)
ggsave(pp,file="GO_BP_top10.pdf")


library(fmsb)
df<-read.csv("synergy_genes_all_enrich_pathway_ra_kegg.csv",header = T,check.names = F,row.names = 1)
pdf("Eneichment_KEGG_Reactome.pdf",width =8, height =8)
radarchart(df[1:5,],
           pty = c(16,16,16),
           axistype = 1,
           pcol = c("#4B114A", "#55998C","#E4AD3D"),
           #pfcol=c(scales::alpha(c("#4B114A", "#55998C","#E4AD3D"),c(0.5,0.5,0.5))),
           plwd = c(3,3,3),
           plty = 1, 
           cglcol = "grey60", 
           cglty = 1, 
           cglwd = 1, 
           axislabcol = "grey60", 
           vlcex = 0.8,
           vlabels = colnames(df),
           caxislabels = c(0, 10, 20, 30, 40),
           calcex=0.8 
)

legend(x = "topleft", legend = c("GPR1_PTX4","IL1RAPL1_PTX4","OLFM3_PTX4"), horiz = F,  bty = "n", pch = 15 , col = c("#4B114A", "#55998C","#E4AD3D"),  text.col = "black", cex = 0.8, pt.cex = 1)
dev.off()

###pathway activity
pathwayActive<-readRDS("synergy_cel_pathways.rds")
df <- t(as.matrix(pathwayActive@assays$pathwaysmlm@data)) %>%
  as.data.frame() %>%
  dplyr::mutate(cluster = Seurat::Idents(pathwayActive)) %>%
  tidyr::pivot_longer(cols = -cluster, 
                      names_to = "source", 
                      values_to = "score") %>%
  dplyr::group_by(cluster, source) %>%
  dplyr::summarise(mean = mean(score))

# Transform to wide matrix
top_acts_mat <- df %>%
  tidyr::pivot_wider(id_cols = 'cluster', 
                     names_from = 'source',
                     values_from = 'mean') %>%
  tibble::column_to_rownames(var = 'cluster') %>%
  as.matrix()
colors <- rev(c('#8c510a','white','#35978f'))
colors.use <- grDevices::colorRampPalette(colors = colors)(100)
my_breaks <- c(seq(-0.3, 0, length.out = ceiling(100 / 2) + 1),
               seq(0.01, 0.3, length.out = floor(100 / 2)))

top_acts_mat=data.frame(top_acts_mat)
top_acts_mat$Androgen=NULL
top_acts_mat$Estrogen=NULL
pheatmap::pheatmap(mat = top_acts_mat,
                      color = colors.use,
                      cluster_rows = F,
                      border_color = "white",
                      breaks = my_breaks,
                      cellwidth = 15,
                      cellheight = 15,
                      treeheight_row = 20,
                      treeheight_col = 20) 




