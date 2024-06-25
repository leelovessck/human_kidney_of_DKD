rm(list = ls())
gc()



#打开必要的package
{
  if(!require(Seurat))install.packages("Seurat")
  if(!require(clustree))install.packages("clustree")
  if(!require(devtools))install.packages("devtools")
  if(!require(presto))devtools::install_github("immunogenomics/presto")
  if(!require(tidyverse))install.packages("tidyverse")
  if(!require(Matrix))install.packages("Matrix")
  if(!require(SeuratDisk))remotes::install_github("mojaveazure/seurat-disk")
  if(!require(Rcpp))install.packages("Rcpp")
  if(!require(RcppEigen))install.packages("RcppEigen")
  if(!require(ggsci))install.packages("ggsci")
  if(!require(viridis))install.packages("viridis")
  if(!require(tidyverse))install.packages("tidyverse")
  if(!require(scibetR))devtools::install_github("zwj-tina/scibetR")
  if(!require(ggplot2))install.packages("ggplot2")
  if(!require(SeuratDisk))remotes::install_github("mojaveazure/seurat-disk")
}



setwd("D:/2023.10/调试/糖尿病肾病/1.人DKD分析")
getwd()
gc()
human_all <- readRDS("./data source/merge_rds/全部样本（分群）.rds")
Idents(human_all) <- human_all$subclass.l1
inter <- subset(human_all, subclass.l1 == "Interstitial")
table(inter$subclass.l1)
table(inter$subclass.l2)



####降维###########################
#降维
inter <- inter %>% 
  NormalizeData(normalization.method = "LogNormalize",
                scale.factor = 10000,
                margin = 1,
                verbose = TRUE) %>%
  FindVariableFeatures(selintertion.method = "vst",
                       nfeatures = 2000)
inter <- ScaleData(inter)
inter <- RunPCA(inter,
             features = VariableFeatures(object = inter))

inter<- JackStraw(inter,
               num.replicate = 100,
               dims = 40)
inter <- ScoreJackStraw(inter, dims = 1:40)
p1 <- JackStrawPlot(inter, dims = 1:40)
p2 <- ElbowPlot(inter, ndims = 50) 
p1 + p2
  ##得到1.间质细胞降维结果图



####分群############################
#分群
inter <- FindNeighbors(inter, dims = 1:30)


##设置不同resolutions查看效果
object_for_clustree <- inter
set.resolutions = seq(0, 1.2, by = 0.1)
object_for_clustree <- FindClusters(object = object_for_clustree ,
                                    resolution = set.resolutions,
                                    verbose = TRUE) 
clustree(object_for_clustree)
  ##得到2.间质细胞分群树状图


##分群
inter <- FindClusters(inter, resolution = 0.7)
inter <- RunUMAP(inter, dims = 1:30)


##间质细胞细胞marker（cluster）
markers.all <- FindAllMarkers(inter,
                              only.pos = TRUE,
                              min.pct = 0.25, logfc.threshold = 0.25,
                              verbose = TRUE)
markers.all <- markers.all %>% group_by(cluster)
write.csv(markers.all,file = "间质细胞各cluster的markers.csv")
markers.top20 <- markers.all %>% group_by(cluster) %>% 
  top_n(n = 20, wt = avg_log2FC)
write.csv(markers.top20,file = "间质细胞各cluster的markers（top20）.csv")



####间质细胞cluster的marker####################################
#间质细胞cluster的marker
inter.marker <- c("COL1A1", "COL1A2", "COL6A3", "COL5A1", "SLC24A3", "DCN",
                  "MYH11", "RGS6", "MCAM", "ROBO1", "PIEZO2", "POSTN", "REN", 
                  "GRID2", "TRPC", "PDGFRA", "PDGFRB", "ITGA8", "CCN1", 
                  "SLIT3", "ACTA2", "RGS5", "ADGRB3")

p1 <- FeaturePlot(inter, features = inter.marker,raster=FALSE)
p1
  ##得到3.UMAP图-间质marker
p2 <- DotPlot(inter, features = inter.marker)+coord_flip()
p2
  ##得到4.间质marker气泡图
p3 <- VlnPlot(inter, features = inter.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p3
  ##得到5.间质marker小提琴图

table(Idents(inter))



####写入注释##########################
#写入注释
new.cluster.ids <- c("0" = "MC",
                     "1" = "VSMC",
                     "2" = "VSMC",
                     "3" = "FIB",
                     "4" = "FIB",
                     "5" = "MC")
names(new.cluster.ids) <- levels(inter)
inter <- RenameIdents(inter, new.cluster.ids)
inter$subclass.l2 <- Idents(inter)
Idents(inter) <- inter$subclass.l2
saveRDS(inter, "./data source/merge_rds/间质细胞亚群.rds")




####间质细胞可视化############################
#间质细胞可视化
Idents(inter) <- inter$seurat_clusters
DimPlot(inter,
        reduction = "umap", 
        label = TRUE, pt.size = 0.5,
        raster = FALSE)
  ##得到6.间质细胞UMAP-细胞分群

Idents(inter) <- inter$orig.ident
DimPlot(inter,
        reduction = "umap", 
        label = FALSE, pt.size = 0.5,
        raster = FALSE)
  ##得到7.间质细胞UMAP-样本

Idents(inter) <- inter$dataSource
DimPlot(inter,
        reduction = "umap", 
        label = FALSE, pt.size = 0.5,
        raster = FALSE)
  ##得到8.间质细胞UMAP-数据库

Idents(inter) <- inter$sampletype
DimPlot(inter,
        reduction = "umap", 
        label = FALSE, pt.size = 0.5,
        raster = FALSE)
  ##得到9.间质细胞UMAP-分组

Idents(inter) <- inter$sex
DimPlot(inter,
        reduction = "umap", 
        label = FALSE, pt.size = 0.5,
        raster = FALSE)
  ##得到10.间质细胞UMAP-性别

Idents(inter) <- inter$subclass.l2
DimPlot(inter,
        reduction = "umap", 
        label = TRUE, pt.size = 0.5,
        raster = FALSE,
        cols = c("#DC050C", "#1965B0", "#882E72", "#33A02C", "#ffff00"))
  ##得到11.间质细胞UMAP-注释



####间质细胞皮尔逊检验##########################
#间质细胞皮尔逊检验
table(inter$subclass.l2)
human_pearson <- AverageExpression(inter,
                                   group.by = "subclass.l2",
                                   assays = "RNA")
human_pearson <- human_pearson[[1]]
human_pearson <- as.matrix(human_pearson)
gene_pearson <- names(tail(sort(apply(human_pearson, 1, sd)),1000))
cell_pearson <- cor(human_pearson[gene_pearson,],method = 'spearman')
pheatmap::pheatmap(cell_pearson,
                   method = "spearman")
  ##得到12.间质细胞注释的spearman相关性



####可视化#########################
#可视化
##系膜细胞
MC.marker <- c("ROBO1", "PIEZO2", "POSTN", "REN", 
               "GRID2", "TRPC", "PDGFRA", "PDGFRB", "ITGA8", "CCN1", 
               "SLIT3", "ACTA2", "RGS5", "ADGRB3")
Idents(inter) <- inter$subclass.l2
p4 <- FeaturePlot(inter, features = MC.marker,raster=FALSE)
p5 <- DotPlot(inter, features = MC.marker)+coord_flip()
p6 <- VlnPlot(inter, features = MC.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p4 / p5 / p6
  ##得到13.系膜细胞marker（celltype）


##血管平滑肌细胞
VSMC.marker <- c("MYH11", "RGS6", "MCAM")
Idents(inter) <- inter$subclass.l2
p4 <- FeaturePlot(inter, features = VSMC.marker,raster=FALSE)
p5 <- DotPlot(inter, features = VSMC.marker)+coord_flip()
p6 <- VlnPlot(inter, features = VSMC.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p4 / p5 / p6
  ##得到14.血管平滑肌细胞marker（celltype）


##成纤维细胞
FIB.marker <- c("COL1A1", "COL1A2", "COL6A3", "COL5A1", "SLC24A3", "DCN")
Idents(inter) <- inter$subclass.l2
p4 <- FeaturePlot(inter, features = FIB.marker,raster=FALSE)
p5 <- DotPlot(inter, features = FIB.marker)+coord_flip()
p6 <- VlnPlot(inter, features = FIB.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
p4 / p5 / p6
  ##得到15.成纤维细胞marker（celltype）
