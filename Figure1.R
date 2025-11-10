

library(ggplot2)
library(ggpubr)
library(Seurat)
library(ComplexHeatmap)
library(colorRamps)
library(circlize) 
library(dplyr)
## Fig1b
techidat<-read.table("gRNA_cells.txt",header = T,sep="\t")

techidat$cellnum=techidat$cells/100000

pp=ggscatter(techidat,x="cellnum",
             y="percent",
             fill = "Expriment",color = "Expriment",
             palette = "npg",
             repel = T,
             legend=NULL,
             label = "Expriment",
             label.select = techidat$Expriment,
             xlab="Number of cells(x10^5)",
             ylab="Cell with sgRNA(%)",
             size = 3)

ggsave(pp,file="Fig1b.pdf",width = 6,height = 4)

##  Fig1c
dat<-readRDS("dat.rds")
per_cell_df <- dat@meta.data[, c("cell_new_name", "Pertubation")]
cell_pert_df <- separate_rows(per_cell_df, Pertubation, sep = ":")
cell_pert_df <- data.frame(cell_pert_df)
cell_pert_df <- unique(cell_pert_df)

pt_cell_num <- data.frame(table(cell_pert_df$Pertubation))
colnames(pt_cell_num) <- c("perturbation", "cells")
pt_cell_num<-read.csv("cell_number_for_each_perturbation.csv",header = T,row.names = 1,check.names = F)

p_hist <- gghistogram(
  pt_cell_num,
  x = "cells",
  bins = 20000,
  xlim = c(0, 1200),
  xlab = "# Number of cells for each perturbation",
  ylab = "# gRNA counts",
  col = NA,
  fill = "#1E4768",
  alpha = 0.8
) +geom_vline(xintercept = c(50), 
              linetype = "dashed", 
              color = "darkgrey")



ggsave(p_hist, file = "Fig1c.pdf", width = 3.5, height = 3)

### Heatmap of effect size 

highlight_pts <- c(
  "DSCAML1","FOXS1","FOXH1","IL1RL1","LPS1",
  "PEG10","TRAF3","GATA4","ASB11","EPB41L1",
  "MFGE8","POSTN","INPPL1","MAP2K3","C1QL1",
  "C1QL4","PTX3"
)

highlight_genes<-c(
  "CAV1","CCL2","MYC","INPPL1","STAT3","COL6A1",
  "PTX3","EGF","NID2","COL25A1","IL6","LOX","FGF7",
  "COL4A5","IL1RL1","SMOC2", "POSTN","ETV1","TNC",
  "MMP3","FOXO1","PRRX1","RUNX2","GATA4",
  "FOXB1","F3","COL5A2","AREG"
)

effectMat<-read.csv('Effect_size_filtered.csv',
                    header = T,
                    check.names = F,
                    row.names = 1)
geneProgram<-read.csv("geneProgram.csv",
                      header = T,
                      check.names = F)
pertubation<-read.csv("perturbationModule.csv",
                      header = T,
                      check.names = F)

pt=rownames(effectMat)
ptindex<-which(pt%in%highlight_pts)
ptlabels<-pt[ptindex]


gp=colnames(effectMat)
gpindex<-which(gp%in%highlight_genes)
gplabs<-gp[gpindex]

colanno=data.frame(Program=geneProgram$Program,row.names =rownames(geneProgram))
rowanno=data.frame(Module=pertubation$Module,row.names =rownames(pertubation))

rowcol <- list(
  Module = c(
    Module1 = "#003C30",
    Module2 = "#43AC5E",
    Module3 = "#2D4999",
    Module4 = "#C9B2D6",
    Module5 = "#E6780C",
    Module6 = "#7E9A29",
    Module7 = "#C51E7E",
    Module8 = "#723F91"
  )
)


colcol <- list(
  Program = c(
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
)

rowanno<-rowAnnotation(Module=rowanno$Module,
                       col=rowcol,
                       Perturbation=anno_mark(at=ptindex,
                                              labels=ptlabels,
                                              labels_gp=gpar(fontsize=10),
                                              lines_gp=gpar()
                                              )
                       
)

colanno = HeatmapAnnotation(Genes=anno_mark(at=gpindex,
                                            labels=gplabs,
                                            labels_gp=gpar(fontsize=10),
                                            side = "top",
                                            lines_gp=gpar()),
                            Program = colanno$Program,
                            col=colcol)

color =colorRamp2(c(-0.2,-0, 0.2), c("#00768B","white", "#C36600"))

pdf("Pertubation_effect_size_heatmap.pdf",width =25, height =20)

effectMat[effectMat>0.2]=0.2
effectMat[effectMat<(-0.2)]=-0.2

p1=ComplexHeatmap::Heatmap(effectMat,
                           show_row_names = FALSE, 
                           name="Effect size",
                           show_column_names = FALSE, 
                           show_column_dend = FALSE,
                           cluster_columns = F,
                           cluster_row_slices =T, 
                           cluster_rows = F,
                           col=color,
                           row_split = pertubation$Module,
                           column_split=geneProgram$Program ,
                           top_annotation =colanno,
                           right_annotation = rowanno,
                           heatmap_height = unit(0.05, "cm")*nrow(mat),
                           heatmap_width =unit(0.01, "cm")*ncol(effectMat))

draw(p1, 
     show_heatmap_legend =T,
     heatmap_legend_side="left", 
     annotation_legend_side="left",
     legend_grouping = "original") # 图例位置

dev.off()

##Density plot for correlation in geneprogram
interProgramCors<-read.csv("GeneProgram_expresion_cors_a1.csv",
                           header = T,
                           row.names = 1,
                           check.names = F)
intraProgramCors<-read.csv("GeneProgram_expresion_cors_a2.csv",
                           header = T,
                           row.names = 1,
                           check.names = F)

dd<-rbind(interProgramCors,intraProgramCors)
write.csv(dd,file="GeneProgram_expression_correlation.csv",quote = F,row.names = F)

geneeProgramCors<-read.csv("GeneProgram_expression_correlation.csv",
                           header = T,
                           check.names = F)


geneeProgramCors<-readRDS("/Volumes/zzzz/CAF_crispr/final_version/articles/Figs/Fig1/Gene_program_expression_cors/gene_cors_all.rds")
geneeProgramCors=data.frame(geneeProgramCors$r)
gghistogram(geneeProgramCors,x="cors",color = "groups",y="density",binwidth=0.01)

p=ggdensity(geneeProgramCors,
            x="cors",
            color = "groups",
            palette = c("#2078B5","#D72829")
)
p
ggsave(p,file="GeneProgram_expression_correlation_dencity.pdf",height = 3,width = 4.5)




## Density plot for effect size in backgaroud and modules
perturbationEffctsize<-read.csv("Perturbation_effect_density.csv",
                                header = T,
                                check.names = F,
                                row.names = 1)

p=ggdensity(perturbationEffctsize,
            x="effectSize",
            color = "groups",
            palette = c("#91B1DE","#F18516")
            )
p
ggsave(p,file="Perturbation_expression_dencity.pdf",height = 3,width = 4.5)


## geneProgram UMAP
geneProgram<-read.csv("geneProgram.csv",header = T,check.names = F)

showGenes<- c(
  "DHX9","YWHAQ","MCM6","SF3B1","NCL","MLH1","POLR2B",
  "PPP2CA","NPM1","HDAC2","CUL1","YWHAZ","TFDP1",
  "HNRNPC","SRSF1","CAPZB","LMNA","ARF1","DCTN1",
  "PSMB1","CYCS","YWHAG","GSN","TUBA1A","CTSA","ATP2B4","LPIN1",
  "HADHA","HADHB","KLKB1","HSD17B4","TFPI2","MYOF",
  "CD9","SMAD3","WNT2B","DDR2","TGFB2","COL3A1","WNT5A","PLOD2",
  "FGF2","LOX","SOX4","DAAM2","COL12A1","COL10A1",
  "CDC42","LYST","MAP4K3","KIF13A","PRKN","NTRK2","FGF7",
  "VPS13C","VPS35","ITGA9","CCK","SUCNR1","CXCL13",
  "SCIN","CNTNAP2","VIPR2","TACR2","BDKRB1","MAP2K1","CLDN6",
  "ITGAM","PVR","PDGFB","KIT","PDGFC","VEGFC","GRIA1",
  "EGFR","SHANK2","SYT1","IGF1","NTRK3","DLGAP1","NEGR1",
  "ERBB4","ADGRL3","EREG","PTPN13","GRID2","RASA1","ELMO1",
  "DOCK1","HIF1A","PRKACB","CDC42EP3","ERLEC1","GCC2",
  "GOLGA4","AREG","DKK2","PEBP1","PPP2R5E","CAV1","CCL2",
  "MYC","INPPL1","STAT3","COL6A1","PTX3","EGF","NID2",
  "COL25A1","IL6","COL4A5","IL1RL1","SMOC2","POSTN",
  "ETV1","TNC","MMP3","FOXO1","PRRX1","RUNX2","GATA4","FOXB1",
  "F3","COL5A2"
)

geneProgramColor <-c( "#4EBCD5", "#E74C35","#07A087", "#A75728",
                      "#F193A8","#BDBD21", "#FABE70", "#3D5589", "#B19D86")

p=ggscatter(geneProgram,
            x="UMAP_1",
            y="UMAP_2",
            color = "Program",
            size=0.3,
            label = "gene",
            label.select = showGenes,
            palette=geneProgramColor 
            )
ggsave(p,file="GeneProgram_UMAP.pdf",height = 3,width = 4.5)



## Perturbation module UMAP plot
modulecol <- c("#003C30", "#43AC5E", "#2D4999", "#C9B2D6", 
               "#E6780C", "#7E9A29","#C51E7E", "#723F91")
p=ggscatter(pertubation,
            x="UMAP_1",
            y="UMAP_2",
            color = "Module",
            size=2,
            label = "Pertubation",
           label.select =highlight_pts ,
            palette=modulecol )

ggsave(p,file="Module_UMAP.pdf",height = 3,width = 4.5)

dat<-readRDS("CAFs.rds")


controlcell<-subset(dat,Pertubation=="Control")

geneProgram$gene1=geneProgram$gene
geneProgram$gene2=geneProgram$gene
geneProgram$Program1=geneProgram$Program
geneProgram$Program2=geneProgram$Program

controlcell[["RNA"]]$counts <- as(object = controlcell[["RNA"]]$counts, Class = "dgCMatrix")
controlcell<-NormalizeData(controlcell)
Idents(controlcell)=controlcell$seurat_clusters

matexpr=AverageExpression(controlcell)

matexpr=PseudobulkExpression(controlcell,method = "average",layer = "data")
matexpr=data.frame(matexpr,check.names = F)

#matexpr<-data.frame(GetAssayData(controlcell,slot = "data"),check.names = F)
matexpr=matexpr[unique(geneeProgramCors1$gene1),]
matexpr1=t(matexpr)

genecors=cor(matexpr1,method = "spearman")
genecors=cor(matexpr1,method = "pearson")
genecors=genecors[geneProgram$gene,geneProgram$gene]
#genecors[genecors>0.5]=0.5
#genecors[genecors<(-0.5)]=-0.5

pheatmap(genecors,cluster_rows = F,cluster_cols = F,show_rownames = F,show_colnames = F)

genecors=data.frame(genecors,check.names = F)
genecors$gene1=rownames(genecors)

genecorslong=reshape2::melt(genecors)
colnames(genecorslong)=c("gene1","gene2","cors")

genecorslong=left_join(genecorslong,geneProgram[,c("gene1","Program1")],by="gene1")

genecorslong=left_join(genecorslong,geneProgram[,c("gene2","Program2")],by="gene2")
genecorslong$groups="inter-module"
genecorslong$groups[genecorslong$Program1==genecorslong$Program2]="intra-module"

ggdensity(genecorslong$value)

pp=gghistogram(genecorslong,x="cors",y="density",
               color = "groups",fill = "groups",
               legend="right",
               palette = c(c("#2078B5","#D72829")))

+ggthemes::theme_few()


ggsave(pp,file="gene_program_expression_correlation.pdf",width = 6,height = 4)
write.csv(genecorslong,file="gene_program_expression_correlation.csv")









