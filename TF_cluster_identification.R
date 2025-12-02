
library(NbClust)
library(ggoubr)
library(ggplot2)

res1<-NbClust(NTC_TF_correlation, distance = "euclidean", min.nc=6, max.nc=15, 
             method = "complete", index = "cindex")
res2<-NbClust(NTC_TF_correlation, distance = "euclidean", min.nc=6, max.nc=15, 
              method = "complete", index = "sdindex")

cindexDf<-data.frame(K=6:15,cindex=res1$All.index,silhouette=res2$All.index)
cindexDf$K=paste0("TFc",cindexDf$K)
cindexDf=reshape2::melt(cindexDf)
colnames(cindexDf)<-c("cluster","parameter","value")
pp=ggline(cindexDf,x="cluster",y="value",color = "parameter",palette = c("#0F6056",  "#F25AA6"))+geom_vline(xintercept="TFc9")+ggthemes::theme_few()

ggsave(pp,file="TFcluster_determination.pdf",width = 7,height = 5)
write.csv(cindexDf,file="TFcluster_determination.csv")




