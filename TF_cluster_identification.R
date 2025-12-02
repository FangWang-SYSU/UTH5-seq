###



library(cluster)
D <- as.dist(1 - abs(NTC_TF_correlation))   # 或 1 - abs(cor_mat)
D <- as.dist(1-NTC_TF_correlation) # 或 1 - abs(cor_mat)
D <- dist(NTC_TF_correlation) 
hc <- hclust(D, method="complete")
ks <- 2:15
sil <- sapply(ks, function(k){
  cl <- cutree(hc, k)
  mean(silhouette(cl, D)[,3],na.rm=T)
})
plot(ks, sil, type="b")
which.max(sil)  # 给出最优K


sil_pac <- Silhouette(dist_matrix, method = "pac", sort = TRUE)
head(sil_pac)



library(cluster)

# 1) 相似度 -> 距离
C <- as.matrix(NTC_TF_correlation)
C[is.na(C)] <- 0
C <- pmax(pmin(C, 1), -1)       # 夹到 [-1, 1]
diag(C) <- 1
D <- as.dist(1 - abs(C))        # 若要保留方向，用 as.dist(1 - C)

# 2) 层次聚类
hc <- hclust(D, method = "complete")  # 也可试 "average"

# 3) 扫描 K 计算平均轮廓系数
ks  <- 2:20
sil <- vapply(ks, function(k){
  cl <- cutree(hc, k)
  mean(silhouette(cl, D)[, "sil_width"])
}, numeric(1))

# 4) 可视化 & 最优 K
plot(ks, sil, type = "b")
ks[which.max(sil)]






library(ComplexHeatmap)
library(cluster)

# cor_mat: 你的 TF×TF 相关矩阵
C <- as.matrix(NTC_TF_correlation)
diag(C) <- 1
C[is.na(C)] <- 0

# 与热图一致：忽略方向的距离（常用）或保留方向(看你需求)
dist_fun <- function(m) as.dist(1 - abs(m))   # 若想保留方向改成: as.dist(1 - m)
dist_fun <- function(m) as.dist(1 - m)   # 若想保留方向改成: as.dist(1 - m)

# —— 画热图时这样写，保证和下面评估用的是同一度量/方法 ——
ht <- Heatmap(
  C,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  # clustering_distance_rows = dist_fun,
  # clustering_distance_columns = dist_fun,
  clustering_method_rows = "complete",         # 或你实际用的 "complete"/"average"
  clustering_method_columns = "complete"
)
ht_draw <- draw(ht)
row_dnd <- row_dend(ht_draw)    # 和热图一致的树
D <- dist_fun(C)                # 与热图一致的距离

ks <- 2:20
# 把树裁成不同K的标签
sil <- sapply(ks, function(k){
  cl <- cutree(as.hclust(row_dnd), k = k)
  mean(silhouette(cl, D)[, "sil_width"])
})
plot(ks, sil, type="b", xlab="K", ylab="Average Silhouette")
ks[which.max(sil)]   # 候选最优K














library(cluster)
D <- as.dist(1 - abs(NTC_TF_correlation))   # 或 1 - abs(cor_mat)
D <- dist(NTC_TF_correlation) # 或 1 - abs(cor_mat)
hc <- hclust(D, method="complete")
ks <- 2:20
sil <- sapply(ks, function(k){
  cl <- cutree(hc, k)
  mean(silhouette(cl, D)[,3])
})
plot(ks, sil, type="b")
which.max(sil)  # 给出最优K


sil_pac <- Silhouette(dist_matrix, method = "pac", sort = TRUE)
head(sil_pac)



BiocManager::install("ConsensusClusterPlus")
library(ConsensusClusterPlus)
results = ConsensusClusterPlus(dist(NTC_TF_correlation),maxK=20,reps=10,pItem=0.8,pFeature=1,   
                               title="./ConsensusCluster",clusterAlg="hc",
                               seed=123456,plot="png")

maxK = 20
Kvec = 2:maxK
x1 = 0.1; x2 = 0.9 # threshold defining the intermediate sub-interval
PAC = rep(NA,length(Kvec))
names(PAC) = paste("K=",Kvec,sep="") # from 2 to maxK

for(i in Kvec){
  M = results[[i]]$consensusMatrix
  Fn = ecdf(M[lower.tri(M)])
  PAC[i-1] = Fn(x2) - Fn(x1)
}#end for i

# The optimal K
optK = Kvec[which.min(PAC)]
optK## [1] 6









library(cluster)

hc <- hclust(D, method = "average")  # 或 "complete"；与生信里常用的平均链接更稳
ks <- 2:20
asw_curve <- sapply(ks, function(k){
  cl_k <- cutree(hc, k)
  mean(silhouette(cl_k, D)[, "sil_width"])
})
plot(ks, asw_curve, type="b", xlab="K", ylab="Average silhouette")
abline(v = 9, lty = 2)
ks[which.max(asw_curve)]           # 统计学上的“最优 K”
asw_curve[ks == 9]                 # K=9 的 ASW（便于与最优K比较）




###
###
library(ComplexHeatmap)
library(cluster)

# 你的相似性矩阵
S <- as.matrix(NTC_TF_correlation)
S[is.na(S)] <- 0

## 1) 画热图（保持你的原设置，不指定 clustering_distance_*）
ht <- Heatmap(
  S,
  show_row_names = FALSE,
  show_column_names = FALSE,
  show_column_dend = FALSE,
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  #col = color,
  name = "Correlation",
  #row_split = 9,
  #column_split = 9,
  top_annotation = colanno,
  right_annotation = rowanno
)
htd <- draw(ht)

## 2) 复现 Heatmap 用的“二阶”距离：把每一行当特征向量，计算行间 dist()
#    这与未指定 clustering_distance_rows 时 Heatmap 内部的做法一致（欧氏距离）
D_row <- dist(S, method = "euclidean")   # 行向量的欧氏距离
row_tree <- as.hclust(row_dend(htd))     # 与热图一致的行树

## 3) 扫描 K 的平均轮廓系数（在同一棵树、同一距离上评估）
ks  <- 2:20
asw <- sapply(ks, function(k){
  cl_k <- cutree(row_tree, k = k)
  mean(silhouette(cl_k, D_row)[, "sil_width"])
})

# 结果与可视化
k_best  <- ks[which.max(asw)]
asw_best <- max(asw)
asw_k9   <- asw[ks == 9]

plot(ks, asw, type = "b", xlab = "K", ylab = "Average Silhouette")
abline(v = 9, lty = 2, col = "gray40")

cat("Best K =", k_best, "  ASW_best =", round(asw_best, 4),
    "\nASW at K=9 =", round(asw_k9, 4), "\n")



library(NbClust)


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




