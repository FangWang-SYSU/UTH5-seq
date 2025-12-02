library(CytoTRACE2)
library(Seurat)
library(ggpubr)


outdir<-"~/work"
setwd(outdir)

dat<-readRDS("CAFs.rds")
controlcell<-subset(dat,Pertubation=="Control")
controlcell[["RNA"]]$counts <- as(object = controlcell[["RNA"]]$counts, Class = "dgCMatrix")
controlcell<-NormalizeData(controlcell)

phenodat<-controlcell$seurat_clusters
phenodat=data.frame(phenodat)
names(phenodat)<-gsub("-",".",names(phenodat),fixed = T)
int.embed <- data.frame(Embeddings(controlcell, reduction = "umap"))
rownames(int.embed)<-gsub("-",".",rownames(int.embed),fixed = T)
annotation<-controlcell$cell_states
results <- cytotrace2(controlcell, ncores = 5, is_seurat = TRUE, 
                      slot_type = "counts",
                      batch_size = 50000,
                      smooth_batch_size=3000,
                      seed = 20,
                      species = 'human')
plotdat<-data.frame(cytotrace_score=results$CytoTRACE2_Score,groups=results$cell_states)

medianRes<-plotdat%>% group_by(groups) %>%
  dplyr::summarize(median_value = median(cytotrace_score)) %>%
  filter(median_value == max(median_value))

#plotdat$group<-phenodat

plotdat<- plotdat %>% mutate(group=forcats::fct_reorder(groups,cytotrace_score, .fun = 'median',.desc = T) )

p=ggviolin(plotdat, x = "group", y = "cytotrace_score",color='group', fill = "group", #palette = zmcols,
           legend="right",
           add = "boxplot", add.params = list(fill = "white"), error.plot = "errorbar",border="white")

p