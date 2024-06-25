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
table(human_all$subclass.l1)
table(human_all$subclass.l2)



####处理免疫细胞###################
#处理免疫细胞
IMM <- subset(human_all, subclass.l1 == "Immune")
table(IMM$subclass.l1)
table(IMM$subclass.l2)

subclass.l2 <- IMM@meta.data
subclass.l2 <- subclass.l2 %>%  
  mutate(subclass.l2 = case_when(  
    subclass.l2 %in% c("cDC", "pDC") ~ "DC",  
    subclass.l2 == "MAC-M2" ~ "MAC",  
    subclass.l2 == "MDC" ~ "MON",  
    subclass.l2 == "ncMON" ~ "MON",  
    subclass.l2 %in% c("NK1", "NK2") ~ "NK",  
    subclass.l2 %in% c("T-CYT", "T-REG", "cycT") ~ "T",  
    subclass.l2 %in% "NKT" ~ "NK",
    TRUE ~ subclass.l2
  ))
table(subclass.l2$subclass.l2)

IMM@meta.data <- subclass.l2
table(IMM$subclass.l1)
table(IMM$subclass.l2)
Idents(IMM) <- IMM$subclass.l2

unique(Idents(human_all))
Idents(human_all, cells = colnames(IMM)) <- Idents(IMM)
DimPlot(human_all, 
        reduction = "umap", 
        label = TRUE, 
        pt.size = 0.5,raster=FALSE)



####处理内皮细胞###################
#处理内皮细胞
EC <- readRDS("./data source/merge_rds/EC部位亚群.rds")
Idents(EC) <- EC$subclass.l2
unique(Idents(EC))
Idents(human_all, cells = colnames(EC)) <- Idents(EC)
DimPlot(human_all, 
        reduction = "umap", 
        label = TRUE, 
        pt.size = 0.5,raster=FALSE)



####处理间质细胞###################
#处理间质细胞
inter <- readRDS("./data source/merge_rds/间质细胞亚群.rds")
Idents(inter) <- inter$subclass.l2
unique(Idents(inter))
Idents(human_all, cells = colnames(inter)) <- Idents(inter)
DimPlot(human_all, 
        reduction = "umap", 
        label = TRUE, 
        pt.size = 0.5,raster=FALSE)



####处理足细胞####################
#处理足细胞
POD <- subset(human_all, subclass.l1 == "PEC")
table(POD$subclass.l1)
table(POD$subclass.l2)

POD$subclass.l1 <- "POD"
POD$subclass.l2 <- "POD"

table(POD$subclass.l1)
table(POD$subclass.l2)
Idents(POD) <- POD$subclass.l2
unique(Idents(human_all))
Idents(human_all, cells = colnames(POD)) <- Idents(POD)
DimPlot(human_all, 
        reduction = "umap", 
        label = TRUE, 
        pt.size = 0.5,raster=FALSE)



####保存对象###################
#保存对象
human_all$usetype <- Idents(human_all)
unique(human_all$usetype)
saveRDS(human_all, "./data source/merge_rds/全部样本（自定义注释）.rds")



####可视化#########################
#可视化
DimPlot(human_all, 
        reduction = "umap", 
        label = TRUE, 
        pt.size = 0.5,raster=FALSE)
  ##得到1.UAMP图-自定义注释

DimPlot(human_all, 
        reduction = "umap", 
        label = FALSE, 
        pt.size = 0.5,raster=FALSE)
  ##得到2.UAMP图-自定义注释（无标签）



####细胞marker（自注释celltype）####################
#细胞marker（自注释celltype）
markers.all <- FindAllMarkers(human_all,
                              only.pos = TRUE,
                              min.pct = 0.25, logfc.threshold = 0.25,
                              verbose = TRUE)
markers.all <- markers.all %>% group_by(cluster)
write.csv(markers.all,file = "自注释celltype的markers.csv")
markers.top20 <- markers.all %>% group_by(cluster) %>% 
  top_n(n = 20, wt = avg_log2FC)
write.csv(markers.top20,file = "自注释celltype的markers（top20）.csv")
