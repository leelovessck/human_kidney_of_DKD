rm(list = ls())
gc()



#打开必要的package
{
if(!require(harmony))devtools::install_github("immunogenomics/harmony")
if(!require(ggplot2))install.packages("ggplot2")
if(!require(multtest))BiocManager::install("multtest")
if(!require(Seurat))install.packages("Seurat")
if(!require(dplyr))install.packages("dplyr")
if(!require(tidyverse))install.packages("tidyverse")
if(!require(clustree))install.packages("clustree")
if(!require(patchwork))install.packages("patchwork")
if(!require(presto))devtools::install_github("immunogenomics/presto")
}



setwd("D:/2023.10/调试/糖尿病肾病/1.人DKD分析")
getwd()
gc()
human_all <- readRDS("./data source/merge_rds/全部样本（原始）.rds")
Idents(human_all) <- human_all$orig.ident



####质量控制###########################
#质量控制
ribo_genes <- rownames(human_all)[grep("^Rp[sl]", 
                                       rownames(human_all),
                                       ignore.case = T)]
human_all <- PercentageFeatureSet(human_all, 
                                  "^RP[SL]", 
                                  col.name = "percent.ribo")

VlnPlot(human_all, 
        features = c("nFeature_RNA", "nCount_RNA", 
                     "percent.mt", "percent.ribo"), 
        ncol = 2,
        pt.size = 0) 
  ##得到1.QC小提琴图（按样本）

p1 <- FeatureScatter(human_all, 
                     feature1 = "nCount_RNA", 
                     feature2 = "percent.mt") +
  NoLegend()
p2 <- FeatureScatter(human_all,
                     feature1 = "nCount_RNA",
                     feature2 = "nFeature_RNA") +
  NoLegend()
p3 <- FeatureScatter(human_all,
                     feature1 = "nCount_RNA",
                     feature2 = "percent.ribo")
p1 + p2 + p3 
  ##得到2.QC散点图（按样本）



####降维###########################
#降维
human_all <- human_all %>% 
  NormalizeData(normalization.method = "LogNormalize",
                scale.factor = 10000,
                margin = 1,
                verbose = TRUE) %>%
  FindVariableFeatures(selection.method = "vst",
                       nfeatures = 2000)
human_all <- ScaleData(human_all)
human_all <- RunPCA(human_all,
                    features = VariableFeatures(object = human_all))

human_all<- JackStraw(human_all,
                      num.replicate = 100,
                      dims = 40)
human_all <- ScoreJackStraw(human_all, dims = 1:40)
p1 <- JackStrawPlot(human_all, dims = 1:40)
p2 <- ElbowPlot(human_all, ndims = 50) 
p1 + p2
  ##得到3.降维结果图



####harmony去批次#######################
#harmony去批次
human_all <- RunHarmony(human_all, "orig.ident")
harmony.embeddings <- Embeddings(human_all, reduction = "harmony")
human_all <- human_all %>% 
  ScaleData(verbose = TRUE) %>% 
  RunPCA(features = VariableFeatures(human_all), npcs = 50, verbose = TRUE)
human_all<- JackStraw(human_all,
                      num.replicate = 100,
                      dims = 40)
human_all <- ScoreJackStraw(human_all, dims = 1:40)
p1 <- JackStrawPlot(human_all, dims = 1:40)
p2 <- ElbowPlot(human_all, ndims = 50) 
p1 + p2
  ##得到4.降维结果图（去批次后）
human_all <- FindNeighbors(human_all, dims = 1:20)



####分群############################
#分群
##设置不同resolutions查看效果
object_for_clustree <- human_all
set.resolutions = seq(0, 1.2, by = 0.1)
object_for_clustree <- FindClusters(object = object_for_clustree ,
                                    resolution = set.resolutions,
                                    verbose = TRUE) 
clustree(object_for_clustree)
  ##得到5.分群树状图


##分群
human_all <- FindClusters(human_all, resolution = 0.6)
human_all <- RunUMAP(human_all, dims = 1:20)
saveRDS(human_all, "./data source/merge_rds/全部样本（分群）.rds")


##细胞marker（cluster）
markers.all <- FindAllMarkers(human_all,
                              only.pos = TRUE,
                              min.pct = 0.25, logfc.threshold = 0.25,
                              verbose = TRUE)
markers.all <- markers.all %>% group_by(cluster)
write.csv(markers.all,file = "各cluster的markers.csv")
markers.top20 <- markers.all %>% group_by(cluster) %>% 
  top_n(n = 20, wt = avg_log2FC)
write.csv(markers.top20,file = "各cluster的markers（top20）.csv")



####分群可视化############################
#分群可视化
DimPlot(human_all,
        reduction = "umap", 
        label = TRUE, pt.size = 0.5,
        raster = FALSE)
  ##得到6.UMAP-细胞分群

Idents(human_all) <- human_all$orig.ident
DimPlot(human_all,
        reduction = "umap", 
        label = FALSE, pt.size = 0.5,
        raster = FALSE)
  ##得到7.UMAP-样本

Idents(human_all) <- human_all$dataSource
DimPlot(human_all,
        reduction = "umap", 
        label = FALSE, pt.size = 0.5,
        raster = FALSE)
  ##得到8.UMAP-数据库

Idents(human_all) <- human_all$sampletype
DimPlot(human_all,
        reduction = "umap", 
        label = FALSE, pt.size = 0.5,
        raster = FALSE)
  ##得到9.UMAP-分组

Idents(human_all) <- human_all$sex
DimPlot(human_all,
        reduction = "umap", 
        label = FALSE, pt.size = 0.5,
        raster = FALSE)
  ##得到10.UMAP-性别

Idents(human_all) <- human_all$subclass.l1
DimPlot(human_all,
        reduction = "umap", 
        label = TRUE, pt.size = 0.5,
        raster = FALSE,
        cols = c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
                 "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#666666", 
                 "#8600bf", "#56ff0d", "#ffff00"))
  ##得到11.UMAP-原注释（粗）

Idents(human_all) <- human_all$subclass.l2
DimPlot(human_all,
        reduction = "umap", 
        label = TRUE, pt.size = 0.5,
        raster = FALSE)
  ##得到12.UMAP-原注释（细）
