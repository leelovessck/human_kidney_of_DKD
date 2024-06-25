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
EC <- subset(human_all, subclass.l1 == "EC")
table(EC$subclass.l1)
table(EC$subclass.l2)



####降维###########################
#降维
EC <- EC %>% 
  NormalizeData(normalization.method = "LogNormalize",
                scale.factor = 10000,
                margin = 1,
                verbose = TRUE) %>%
  FindVariableFeatures(selection.method = "vst",
                       nfeatures = 2000)
EC <- ScaleData(EC)
EC <- RunPCA(EC,
             features = VariableFeatures(object = EC))

EC<- JackStraw(EC,
               num.replicate = 100,
               dims = 40)
EC <- ScoreJackStraw(EC, dims = 1:40)
p1 <- JackStrawPlot(EC, dims = 1:40)
p2 <- ElbowPlot(EC, ndims = 50) 
p1 + p2
  ##得到1.EC部位降维结果图



####分群############################
#分群
EC <- FindNeighbors(EC, dims = 1:15)


##设置不同resolutions查看效果
object_for_clustree <- EC
set.resolutions = seq(0, 1.2, by = 0.1)
object_for_clustree <- FindClusters(object = object_for_clustree ,
                                    resolution = set.resolutions,
                                    verbose = TRUE) 
clustree(object_for_clustree)
  ##得到2.EC部位分群树状图


##分群
EC <- FindClusters(EC, resolution = 0.6)
EC <- RunUMAP(EC, dims = 1:15)


##EC部位细胞marker（cluster）
markers.all <- FindAllMarkers(EC,
                              only.pos = TRUE,
                              min.pct = 0.25, logfc.threshold = 0.25,
                              verbose = TRUE)
markers.all <- markers.all %>% group_by(cluster)
write.csv(markers.all,file = "EC部位各cluster的markers.csv")
markers.top20 <- markers.all %>% group_by(cluster) %>% 
  top_n(n = 20, wt = avg_log2FC)
write.csv(markers.top20,file = "EC部位各cluster的markers（top20）.csv")



####EC部位亚群注释#########################
#EC部位亚群注释
subclass <- EC@meta.data
table(subclass$subclass.l2)
subclass <- subclass %>%  
  mutate(subclass.l2 = ifelse(subclass.l2 == "dEC-PTC",
                              "EC-PTC",
                              subclass.l2))  
table(subclass$subclass.l2)
EC@meta.data <- subclass
table(EC$subclass.l2)
saveRDS(EC, "./data source/merge_rds/EC部位亚群.rds")
EC <- readRDS("./data source/merge_rds/EC部位亚群.rds")


####EC部位可视化############################
#EC部位可视化
Idents(EC) <- EC$seurat_clusters
DimPlot(EC,
        reduction = "umap", 
        label = TRUE, pt.size = 0.5,
        raster = FALSE)
  ##得到3.EC部位UMAP-细胞分群

Idents(EC) <- EC$orig.ident
DimPlot(EC,
        reduction = "umap", 
        label = FALSE, pt.size = 0.5,
        raster = FALSE)
  ##得到4.EC部位UMAP-样本

Idents(EC) <- EC$dataSource
DimPlot(EC,
        reduction = "umap", 
        label = FALSE, pt.size = 0.5,
        raster = FALSE)
  ##得到5.EC部位UMAP-数据库

Idents(EC) <- EC$sampletype
DimPlot(EC,
        reduction = "umap", 
        label = FALSE, pt.size = 0.5,
        raster = FALSE)
  ##得到6.EC部位UMAP-分组

Idents(EC) <- EC$sex
DimPlot(EC,
        reduction = "umap", 
        label = FALSE, pt.size = 0.5,
        raster = FALSE)
  ##得到7.EC部位UMAP-性别

Idents(EC) <- EC$subclass.l2
DimPlot(EC,
        reduction = "umap", 
        label = TRUE, pt.size = 0.5,
        raster = FALSE,
        cols = c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
                 "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#666666", 
                 "#8600bf", "#56ff0d", "#ffff00"))
  ##得到8.EC部位UMAP-原注释



####内皮部位（亚群）皮尔逊检验##########################
#内皮部位（亚群）皮尔逊检验
table(EC$subclass.l2)
human_pearson <- AverageExpression(EC,
                                   group.by = "subclass.l2",
                                   assays = "RNA")
human_pearson <- human_pearson[[1]]
human_pearson <- as.matrix(human_pearson)
gene_pearson <- names(tail(sort(apply(human_pearson, 1, sd)),1000))
cell_pearson <- cor(human_pearson[gene_pearson,],method = 'spearman')
pheatmap::pheatmap(cell_pearson,
                   method = "spearman")
  ##得到9.EC部位注释的spearman相关性



####注释置信分数##############################
#注释置信分数
##获得参考矩阵
human_kpmp <- LoadH5Seurat("D:/2023.10/调试/糖尿病肾病/使用的原始数据/人/KPMP/KPMP.h5Seurat")
Idents(human_kpmp) <- human_kpmp$subclass.l1
table(human_kpmp$subclass.l1)
human_EC <- subset(human_kpmp, subclass.l1 == "EC")
 

table(human_EC$subclass.l1)
human_EC$subclass.l1 <- droplevels(human_EC$subclass.l1,
                                   exclude = setdiff(
                                     levels(human_EC$subclass.l1),
                                     unique(human_EC$subclass.l1)))
table(human_EC$subclass.l2)
human_EC$subclass.l2 <- droplevels(human_EC$subclass.l2,
                                   exclude = setdiff(
                                     levels(human_EC$subclass.l2),
                                     unique(human_EC$subclass.l2)))
Idents(human_EC) <- human_EC$subclass.l2

reference_matrix <- GetAssayData(object = human_EC,
                                 layer = "counts") 
reference_matrix <- as.data.frame(t(reference_matrix))

cell_labels <- human_EC$subclass.l2 
cell_labels_df <- as.data.frame(cell_labels)  
colnames(cell_labels_df) <- "label"  
reference <- cbind(reference_matrix, cell_labels_df) 


##获得注解矩阵
query_matrix <- GetAssayData(object = EC,
                             layer = "counts") 
query_matrix <- as.data.frame(t(query_matrix))

cell_labels <- EC$subclass.l2  
cell_labels_df <- as.data.frame(cell_labels)  
colnames(cell_labels_df) <- "label"  
query <- cbind(query_matrix, cell_labels_df)  
ori_label <- query$label


##计算置信度分数
query <- query[,-ncol(query)]
null <- readRDS("./data source/merge_rds/null.rds.gz")
c_score <- conf_score_R(ref = reference,
                        query = query,
                        null_expr = null,
                        gene_num = 500)
tibble(
  ori = ori_label,
  prd = SciBet_R(reference, query),
  c_score = c_score
) -> res

Confusion_heatmap_negctrl(res, cutoff = 0.05)
  ##得到10.EC部位注释与参考细胞集的c_score相关性



####可视化#########################
#可视化
##肾小球毛细血管内皮细胞
Idents(EC) <- EC$seurat_clusters
GC.marker <- c("PECAM1", "CDH5", "EMCN", "HECW2", "PLAT", "ITGA8", "BTNL9", 
               "CEACAM1", "PITPNC1", "SLCO2A1")
p1 <- FeaturePlot(EC, features = GC.marker,raster=FALSE)
p2 <- DotPlot(EC, features = GC.marker)+coord_flip()
p3 <- VlnPlot(EC, features = GC.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
Idents(EC) <- EC$subclass.l2
p4 <- FeaturePlot(EC, features = GC.marker,raster=FALSE)
p5 <- DotPlot(EC, features = GC.marker)+coord_flip()
p6 <- VlnPlot(EC, features = GC.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)

p1 / p2 / p3
  ##得到11.肾小球毛细血管内皮细胞marker（cluster）
p4 / p5 / p6
  ##得到12.肾小球毛细血管内皮细胞marker（celltype）


##入、出球小动脉内皮细胞
Idents(EC) <- EC$seurat_clusters
AEA.marker <- c("PECAM1", "CDH5", "EMCN", "BTNL9", "PALMD", "TM4SF1", 
                "SERPINE2", "SLC14A1", "ENPP2", "SLCO2A1")
p1 <- FeaturePlot(EC, features = AEA.marker,raster=FALSE)
p2 <- DotPlot(EC, features = AEA.marker)+coord_flip()
p3 <- VlnPlot(EC, features = AEA.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
Idents(EC) <- EC$subclass.l2
p4 <- FeaturePlot(EC, features = AEA.marker,raster=FALSE)
p5 <- DotPlot(EC, features = AEA.marker)+coord_flip()
p6 <- VlnPlot(EC, features = AEA.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)

p1 / p2 / p3
  ##得到13.入、出球小动脉内皮细胞marker（cluster）
p4 / p5 / p6
  ##得到14.入、出球小动脉内皮细胞marker（celltype）


##管周毛细血管内皮细胞
Idents(EC) <- EC$seurat_clusters
PTC.marker <- c("PECAM1", "CDH5", "EMCN", "DNASE1L3", "CEACAM1", "PITPNC1", 
                "SLCO2A1", "PLVAP", "TLL1")
p1 <- FeaturePlot(EC, features = PTC.marker,raster=FALSE)
p2 <- DotPlot(EC, features = PTC.marker)+coord_flip()
p3 <- VlnPlot(EC, features = PTC.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
Idents(EC) <- EC$subclass.l2
p4 <- FeaturePlot(EC, features = PTC.marker,raster=FALSE)
p5 <- DotPlot(EC, features = PTC.marker)+coord_flip()
p6 <- VlnPlot(EC, features = PTC.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)

p1 / p2 / p3
  ##得到15.管周毛细血管内皮细胞marker（cluster）
p4 / p5 / p6
  ##得到16.管周毛细血管内皮细胞marker（celltype）


##淋巴内皮细胞
Idents(EC) <- EC$seurat_clusters
LYM.marker <- c("PECAM1", "CDH5", "EMCN", "MMRN1")
p1 <- FeaturePlot(EC, features = LYM.marker,raster=FALSE)
p2 <- DotPlot(EC, features = LYM.marker)+coord_flip()
p3 <- VlnPlot(EC, features = LYM.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
Idents(EC) <- EC$subclass.l2
p4 <- FeaturePlot(EC, features = LYM.marker,raster=FALSE)
p5 <- DotPlot(EC, features = LYM.marker)+coord_flip()
p6 <- VlnPlot(EC, features = LYM.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)

p1 / p2 / p3
  ##得到17.淋巴内皮细胞marker（cluster）
p4 / p5 / p6
  ##得到18.淋巴内皮细胞marker（celltype）


##循环内皮细胞
Idents(EC) <- EC$seurat_clusters
CYC.marker <- c("PECAM1", "CDH5", "EMCN", "PLAT", "PITPNC1", "SLCO2A1", "PLVAP")
p1 <- FeaturePlot(EC, features = CYC.marker,raster=FALSE)
p2 <- DotPlot(EC, features = CYC.marker)+coord_flip()
p3 <- VlnPlot(EC, features = CYC.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
Idents(EC) <- EC$subclass.l2
p4 <- FeaturePlot(EC, features = CYC.marker,raster=FALSE)
p5 <- DotPlot(EC, features = CYC.marker)+coord_flip()
p6 <- VlnPlot(EC, features = CYC.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)

p1 / p2 / p3
  ##得到19.循环内皮细胞marker（cluster）
p4 / p5 / p6
  ##得到20.循环内皮细胞marker（celltype）
