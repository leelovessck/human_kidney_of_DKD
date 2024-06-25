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
}



setwd("D:/2023.10/调试/糖尿病肾病/1.人DKD分析")
getwd()
gc()
human_all <- readRDS("./data source/merge_rds/全部样本（分群）.rds")
Idents(human_all) <- human_all$subclass.l1




####细胞marker（celltype）####################
#细胞marker（celltype）
markers.all <- FindAllMarkers(human_all,
                              only.pos = TRUE,
                              min.pct = 0.25, logfc.threshold = 0.25,
                              verbose = TRUE)
markers.all <- markers.all %>% group_by(cluster)
write.csv(markers.all,file = "各celltype的markers.csv")
markers.top20 <- markers.all %>% group_by(cluster) %>% 
  top_n(n = 20, wt = avg_log2FC)
write.csv(markers.top20,file = "各celltype的markers（top20）.csv")



####皮尔逊检验##########################
#皮尔逊检验
table(human_all$subclass.l1)
human_pearson <- AverageExpression(human_all,
                                   group.by = "subclass.l1",
                                   assays = "RNA")
human_pearson <- human_pearson[[1]]
human_pearson <- as.matrix(human_pearson)
gene_pearson <- names(tail(sort(apply(human_pearson, 1, sd)),20000))
cell_pearson <- cor(human_pearson[gene_pearson,],method = 'spearman')
pheatmap::pheatmap(cell_pearson,
                   method = "spearman")
  ##得到1.细胞注释的spearman相关性



####注释置信分数##############################
#注释置信分数
##获得参考矩阵
human_kpmp <- readRDS("./data source/merge_rds/标准注释.rds")
Idents(human_kpmp) <- human_kpmp$subclass.l1
subclass.l1 <- human_kpmp@meta.data
subclass.l1 <- subclass.l1 %>%  
  mutate(subclass.l1 = ifelse(subclass.l1 == "ATL/TAL", "ATL", subclass.l1))
human_kpmp@meta.data <- subclass.l1

reference_matrix <- GetAssayData(object = human_kpmp,
                                 layer = "counts") 
reference_seurat <- CreateSeuratObject(counts = reference_matrix,
                                       min.cells = 200)
reference_matrix <- GetAssayData(object = reference_seurat,
                                 layer = "counts")
reference_matrix <- as.data.frame(t(reference_matrix))

cell_labels <- human_kpmp$subclass.l1  
cell_labels_df <- as.data.frame(cell_labels)  
colnames(cell_labels_df) <- "label"  
reference <- cbind(reference_matrix, cell_labels_df) 


##获得注解矩阵
query_matrix <- GetAssayData(object = human_all,
                             layer = "counts") 
query_seurat <- CreateSeuratObject(counts = query_matrix,
                                   min.cells = 200)
query_matrix <- GetAssayData(object = query_seurat,
                             layer = "counts")
query_matrix <- as.data.frame(t(query_matrix))

cell_labels <- human_all$subclass.l1  
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

Confusion_heatmap_negctrl(res, cutoff = 0.01)
  ##得到2.与参考细胞集的c_score相关性



####可视化##############################
#可视化
##足细胞
Idents(human_all) <- human_all$seurat_clusters
POD.marker <- c("PTPRQ", "WT1", "NPHS2", "CDKN1C", "SPOCK2", "IGFBP2")
p1 <- FeaturePlot(human_all, features = POD.marker,raster=FALSE)
p2 <- DotPlot(human_all, features = POD.marker)+coord_flip()
p3 <- VlnPlot(human_all, features = POD.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
Idents(human_all) <- human_all$subclass.l1
p4 <- FeaturePlot(human_all, features = POD.marker,raster=FALSE)
p5 <- DotPlot(human_all, features = POD.marker)+coord_flip()
p6 <- VlnPlot(human_all, features = POD.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)

p1 / p2 / p3
  ##得到3.足细胞marker（cluster）
p4 / p5 / p6
  ##得到4.足细胞marker（celltype）


##肾小球壁层上皮
Idents(human_all) <- human_all$seurat_clusters
PEC.marker <- c("CLDN1", "CFH", "ALDH1A2")
p1 <- FeaturePlot(human_all, features = PEC.marker,raster=FALSE)
p2 <- DotPlot(human_all, features = PEC.marker)+coord_flip()
p3 <- VlnPlot(human_all, features = PEC.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
Idents(human_all) <- human_all$subclass.l1
p4 <- FeaturePlot(human_all, features = PEC.marker,raster=FALSE)
p5 <- DotPlot(human_all, features = PEC.marker)+coord_flip()
p6 <- VlnPlot(human_all, features = PEC.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)

p1 / p2 / p3
  ##得到5.肾小球壁层上皮marker（cluster）
p4 / p5 / p6
  ##得到6.肾小球壁层上皮marker（celltype）


##近端小管
Idents(human_all) <- human_all$seurat_clusters
PT.marker <- c("LRP2", "SLC5A12", "SLC22A6", "PRODH2", "SLC5A2", "SLC22A8", 
               "SLC7A8", "SLC34A1", "SLC5A10", "SLC5A11", "SLC7A13")
p1 <- FeaturePlot(human_all, features = PT.marker,raster=FALSE)
p2 <- DotPlot(human_all, features = PT.marker)+coord_flip()
p3 <- VlnPlot(human_all, features = PT.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
Idents(human_all) <- human_all$subclass.l1
p4 <- FeaturePlot(human_all, features = PT.marker,raster=FALSE)
p5 <- DotPlot(human_all, features = PT.marker)+coord_flip()
p6 <- VlnPlot(human_all, features = PT.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)

p1 / p2 / p3
  ##得到7.近端小管marker（cluster）
p4 / p5 / p6
  ##得到8.近端小管marker（celltype）


##髓袢升支细段
Idents(human_all) <- human_all$seurat_clusters
ATL.marker <- c("CRYAB", "TACSTD2", "SLC44A5", "CLDN10", "SOD3", "PAPPA2")
p1 <- FeaturePlot(human_all, features = ATL.marker,raster=FALSE)
p2 <- DotPlot(human_all, features = ATL.marker)+coord_flip()
p3 <- VlnPlot(human_all, features = ATL.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
Idents(human_all) <- human_all$subclass.l1
p4 <- FeaturePlot(human_all, features = ATL.marker,raster=FALSE)
p5 <- DotPlot(human_all, features = ATL.marker)+coord_flip()
p6 <- VlnPlot(human_all, features = ATL.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)

p1 / p2 / p3
  ##得到9.髓袢升支细段marker（cluster）
p4 / p5 / p6
  ##得到10.髓袢升支细段marker（celltype）


##髓袢降支细段
Idents(human_all) <- human_all$seurat_clusters
DTL.marker <- c("CRYAB", "TACSTD2", "ID1", "AKR1B1")
p1 <- FeaturePlot(human_all, features = DTL.marker,raster=FALSE)
p2 <- DotPlot(human_all, features = DTL.marker)+coord_flip()
p3 <- VlnPlot(human_all, features = DTL.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
Idents(human_all) <- human_all$subclass.l1
p4 <- FeaturePlot(human_all, features = DTL.marker,raster=FALSE)
p5 <- DotPlot(human_all, features = DTL.marker)+coord_flip()
p6 <- VlnPlot(human_all, features = DTL.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)

p1 / p2 / p3
  ##得到11.髓袢降支细段marker（cluster）
p4 / p5 / p6
  ##得到12.髓袢降支细段marker（celltype）


##髓袢升支粗段
Idents(human_all) <- human_all$seurat_clusters
TAL.marker <- c("CASR", "SLC12A1", "UMOD", "PROM1", "ITGB1", "FGF13", "CLDN14",
                "KCTD16", "ANK2", "ESRRB", "EGF", "ENOX1")
p1 <- FeaturePlot(human_all, features = TAL.marker,raster=FALSE)
p2 <- DotPlot(human_all, features = TAL.marker)+coord_flip()
p3 <- VlnPlot(human_all, features = TAL.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
Idents(human_all) <- human_all$subclass.l1
p4 <- FeaturePlot(human_all, features = TAL.marker,raster=FALSE)
p5 <- DotPlot(human_all, features = TAL.marker)+coord_flip()
p6 <- VlnPlot(human_all, features = TAL.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)

p1 / p2 / p3
  ##得到13.髓袢升支粗段marker（cluster）
p4 / p5 / p6
  ##得到14.髓袢升支粗段marker（celltype）


##集合管
Idents(human_all) <- human_all$seurat_clusters
DCT.marker <- c("SLC12A3", "TRPM6")
p1 <- FeaturePlot(human_all, features = DCT.marker,raster=FALSE)
p2 <- DotPlot(human_all, features = DCT.marker)+coord_flip()
p3 <- VlnPlot(human_all, features = DCT.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
Idents(human_all) <- human_all$subclass.l1
p4 <- FeaturePlot(human_all, features = DCT.marker,raster=FALSE)
p5 <- DotPlot(human_all, features = DCT.marker)+coord_flip()
p6 <- VlnPlot(human_all, features = DCT.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)

p1 / p2 / p3
  ##得到15.集合管marker（cluster）
p4 / p5 / p6
  ##得到16.集合管marker（celltype）


##连接管
Idents(human_all) <- human_all$seurat_clusters
CNT.marker <- c("SLC8A1", "SCN2A", "HSD11B2", "CALB1", "TRPV5")
p1 <- FeaturePlot(human_all, features = CNT.marker,raster=FALSE)
p2 <- DotPlot(human_all, features = CNT.marker)+coord_flip()
p3 <- VlnPlot(human_all, features = CNT.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
Idents(human_all) <- human_all$subclass.l1
p4 <- FeaturePlot(human_all, features = CNT.marker,raster=FALSE)
p5 <- DotPlot(human_all, features = CNT.marker)+coord_flip()
p6 <- VlnPlot(human_all, features = CNT.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)

p1 / p2 / p3
  ##得到17.连接管marker（cluster）
p4 / p5 / p6
  ##得到18.连接管marker（celltype）


##主细胞
Idents(human_all) <- human_all$seurat_clusters
PC.marker <- c("SCNN1G", "SCNN1B", "GATA3", "AQP2", "AQP3", "PDE10A", 
               "KCNK13", "FXYD4", "PHACTR1", "SLC14A2")
p1 <- FeaturePlot(human_all, features = PC.marker,raster=FALSE)
p2 <- DotPlot(human_all, features = PC.marker)+coord_flip()
p3 <- VlnPlot(human_all, features = PC.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
Idents(human_all) <- human_all$subclass.l1
p4 <- FeaturePlot(human_all, features = PC.marker,raster=FALSE)
p5 <- DotPlot(human_all, features = PC.marker)+coord_flip()
p6 <- VlnPlot(human_all, features = PC.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)

p1 / p2 / p3
  ##得到19.主细胞marker（cluster）
p4 / p5 / p6
  ##得到20.主细胞marker（celltype）


##闰细胞
Idents(human_all) <- human_all$seurat_clusters
IC.marker <- c("ATP6V0D2", "SLC4A1", "SLC26A7", "KIT", "AQP6", 
               "CALCA", "SLC4A9", "INSRR")
p1 <- FeaturePlot(human_all, features = IC.marker,raster=FALSE)
p2 <- DotPlot(human_all, features = IC.marker)+coord_flip()
p3 <- VlnPlot(human_all, features = IC.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
Idents(human_all) <- human_all$subclass.l1
p4 <- FeaturePlot(human_all, features = IC.marker,raster=FALSE)
p5 <- DotPlot(human_all, features = IC.marker)+coord_flip()
p6 <- VlnPlot(human_all, features = IC.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)

p1 / p2 / p3
  ##得到21.闰细胞marker（cluster）
p4 / p5 / p6
  ##得到22.闰细胞marker（celltype）


##内皮细胞
Idents(human_all) <- human_all$seurat_clusters
EC.marker <- c("PECAM1", "EMCN", "CDH5", "RAMP2", "RGCC", "AQP1", 
               "KDR", "PLVAP", "SLC14A1", "PTPRB")
p1 <- FeaturePlot(human_all, features = EC.marker,raster=FALSE)
p2 <- DotPlot(human_all, features = EC.marker)+coord_flip()
p3 <- VlnPlot(human_all, features = EC.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
Idents(human_all) <- human_all$subclass.l1
p4 <- FeaturePlot(human_all, features = EC.marker,raster=FALSE)
p5 <- DotPlot(human_all, features = EC.marker)+coord_flip()
p6 <- VlnPlot(human_all, features = EC.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)

p1 / p2 / p3
  ##得到23.内皮细胞marker（cluster）
p4 / p5 / p6
  ##得到24.内皮细胞marker（celltype）


##血管平滑肌细胞和周细胞
Idents(human_all) <- human_all$seurat_clusters
VSMP.marker <- c("PDGFRB", "ROBO1", "PIEZO2", "POSTN", "REN", "GRID2", 
                 "MYH11", "RGS6", "MCAM", "RGS5", "ADGRB3")
p1 <- FeaturePlot(human_all, features = VSMP.marker,raster=FALSE)
p2 <- DotPlot(human_all, features = VSMP.marker)+coord_flip()
p3 <- VlnPlot(human_all, features = VSMP.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
Idents(human_all) <- human_all$subclass.l1
p4 <- FeaturePlot(human_all, features = VSMP.marker,raster=FALSE)
p5 <- DotPlot(human_all, features = VSMP.marker)+coord_flip()
p6 <- VlnPlot(human_all, features = VSMP.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)

p1 / p2 / p3
  ##得到25.血管平滑肌细胞和周细胞marker（cluster）
p4 / p5 / p6
  ##得到26.血管平滑肌细胞和周细胞marker（celltype）


##成纤维细胞
Idents(human_all) <- human_all$seurat_clusters
FIB.marker <- c("COL1A1", "COL1A2", "C7", "DCN", "SYNPO2", "PCDH7")
p1 <- FeaturePlot(human_all, features = FIB.marker,raster=FALSE)
p2 <- DotPlot(human_all, features = FIB.marker)+coord_flip()
p3 <- VlnPlot(human_all, features = FIB.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
Idents(human_all) <- human_all$subclass.l1
p4 <- FeaturePlot(human_all, features = FIB.marker,raster=FALSE)
p5 <- DotPlot(human_all, features = FIB.marker)+coord_flip()
p6 <- VlnPlot(human_all, features = FIB.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)

p1 / p2 / p3
  ##得到27.成纤维细胞marker（cluster）
p4 / p5 / p6
  ##得到28.成纤维细胞marker（celltype）


##免疫细胞
Idents(human_all) <- human_all$seurat_clusters
IMM.marker <- c("PTPRC", "IGKC", "IL7R", "CD96", "CD14")
p1 <- FeaturePlot(human_all, features = IMM.marker,raster=FALSE)
p2 <- DotPlot(human_all, features = IMM.marker)+coord_flip()
p3 <- VlnPlot(human_all, features = IMM.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)
Idents(human_all) <- human_all$subclass.l1
p4 <- FeaturePlot(human_all, features = IMM.marker,raster=FALSE)
p5 <- DotPlot(human_all, features = IMM.marker)+coord_flip()
p6 <- VlnPlot(human_all, features = IMM.marker,
              pt.size = 0, layer = "counts",
              log = TRUE ,raster=FALSE)

p1 / p2 / p3
  ##得到29.免疫细胞marker（cluster）
p4 / p5 / p6
  ##得到30.免疫细胞marker（celltype）
