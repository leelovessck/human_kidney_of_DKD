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
  if(!require(MySeuratWrappers))install.packages("MySeuratWrappers")
}



setwd("D:/2023.10/调试/糖尿病肾病/1.人DKD分析")
getwd()
gc()
human_all <- readRDS("./data source/merge_rds/全部样本（自定义注释）.rds")
Idents(human_all) <- human_all$usetype
cellnum <- as.data.frame(table(human_all$usetype))



####设定细胞顺序###############
#设定细胞顺序
human_all$usetype <- factor(x = human_all$usetype,
                            levels = c("POD", "PT", "ATL", "DTL", "TAL",
                                       "DCT", "CNT", "PC", "IC", "EC-PTC",
                                       "EC-GC", "EC-AEA", "EC-LYM", "cycEC",
                                       "MC", "VSMC", "FIB", "NK", "B", "T", 
                                       "MON", "MAC", "DC", "PL", "cycMNP",
                                       "MAST"))
levels(human_all) <- c("MAST", "cycMNP", "PL", "DC", "MAC", "MON", "T", "B",
                       "NK", "FIB", "VSMC", "MC", "cycEC", "EC-LYM", "EC-AEA",
                       "EC-GC", "EC-PTC", "IC", "PC", "CNT", "DCT", "TAL", 
                       "DTL", "ATL", "PT", "POD")
saveRDS(human_all, "./data source/merge_rds/全部样本（自定义注释）.rds")



####设定marker######################
#设定marker
POD.marker <- c("PTPRQ", "WT1", "NPHS2", "CDKN1C", "SPOCK2", "IGFBP2")
PT.marker <- c("LRP2", "SLC5A12", "SLC22A6", "PRODH2", "SLC5A2", "SLC22A8", 
               "SLC7A8")
ATL.marker <- c("CRYAB", "TACSTD2", "SLC44A5", "SOD3", "PAPPA2")
DTL.marker <- c("CRYAB", "TACSTD2", "ID1", "AKR1B1")
TAL.marker <- c("CASR", "SLC12A1", "UMOD", "PROM1", "ITGB1", "FGF13", "CLDN14",
                "KCTD16", "ANK2", "ESRRB", "EGF", "ENOX1")
DCT.marker <- c("SLC12A3", "TRPM6")
CNT.marker <- c("SLC8A1", "SCN2A", "HSD11B2", "CALB1", "TRPV5")
PC.marker <- c("SCNN1G", "SCNN1B", "GATA3", "AQP2", "AQP3", "PDE10A", 
               "KCNK13", "FXYD4", "PHACTR1", "SLC14A2")
IC.marker <- c("ATP6V0D2", "SLC4A1", "SLC26A7", "KIT", "AQP6", 
               "CALCA", "SLC4A9", "INSRR")
EC.marker <- c("PECAM1", "EMCN", "CDH5", "RAMP2", "RGCC", "AQP1", "KDR", 
               "PLVAP", "SLC14A1", "PTPRB", "HECW2", "PLAT", "ITGA8", "BTNL9", 
               "PALMD", "TM4SF1", "SERPINE2", "ENPP2", "DNASE1L3", 
               "CEACAM1", "PITPNC1", "SLCO2A1", "TLL1", "MMRN1")
MC.marker <- c("ROBO1", "PIEZO2", "POSTN", "REN", 
               "GRID2", "TRPC", "PDGFRA", "PDGFRB", "ITGA8", "CCN1", 
               "SLIT3", "ACTA2", "RGS5", "ADGRB3")
VSMC.marker <- c("MYH11", "RGS6", "MCAM")
FIB.marker <- c("COL1A1", "COL1A2", "COL6A3", "COL5A1", "SLC24A3", "DCN")
IMM.marker <- c("PTPRC", "BANK1", "MS4A1", "IGKC", "MZB1", "THEMIS", "IL7R",
                "CD96", "CD247", "GNLY", "NKG7", "GZMA", "MS4A2", "MRC1", 
                "CD163", "CD14", "MSR1", "DIAPH3", "CENPF", "MKI67", "ITGAX", 
                "HLA-DQA1", "CSF2RA", "FLT3", "CLEC9A", "IL3RA", "CLEC4C", 
                "CTSS", "FCN1", "FCGR3A", "S100A9", "S100A8", "FCGR3B")
all.marker <- c("PTPRQ", "WT1", "NPHS2", "CDKN1C", "SPOCK2", "IGFBP2", "LRP2", 
                "SLC5A12", "SLC22A6", "PRODH2", "SLC5A2", "SLC22A8", "SLC7A8", 
                "CRYAB", "TACSTD2", "SLC44A5", "SOD3", "PAPPA2", "ID1", 
                "AKR1B1", "CASR", "SLC12A1", "UMOD", "PROM1", "ITGB1", "FGF13", 
                "CLDN14", "KCTD16", "ANK2", "ESRRB", "EGF", "ENOX1", "SLC12A3", 
                "TRPM6", "SLC8A1", "SCN2A", "HSD11B2", "CALB1", "TRPV5", 
                "SCNN1G", "SCNN1B", "GATA3", "AQP2", "AQP3", "PDE10A", 
                "KCNK13", "FXYD4", "PHACTR1", "SLC14A2", "ATP6V0D2", "SLC4A1", 
                "SLC26A7", "KIT", "AQP6", "CALCA", "SLC4A9", "INSRR", "PECAM1", 
                "EMCN", "CDH5", "RAMP2", "RGCC", "AQP1", "KDR", "PLVAP", 
                "SLC14A1", "PTPRB", "HECW2", "PLAT", "BTNL9", "PALMD", 
                "TM4SF1", "SERPINE2", "ENPP2", "DNASE1L3", "CEACAM1", 
                "PITPNC1", "SLCO2A1", "TLL1", "MMRN1", "ROBO1", "PIEZO2", 
                "POSTN", "REN", "GRID2", "TRPC", "PDGFRA", "PDGFRB", "CCN1", 
                "SLIT3", "ACTA2", "RGS5", "ADGRB3", "MYH11", "RGS6", "MCAM", 
                "COL1A1", "COL1A2", "COL6A3", "COL5A1", "SLC24A3", "DCN", 
                "PTPRC", "BANK1", "MS4A1", "IGKC", "MZB1", "THEMIS", "IL7R", 
                "CD96", "CD247", "GNLY", "NKG7", "GZMA", "MS4A2", "MRC1", 
                "CD163", "CD14", "MSR1", "DIAPH3", "CENPF", "MKI67", "ITGAX", 
                "HLA-DQA1", "CSF2RA", "FLT3", "CLEC9A", "IL3RA", "CLEC4C", 
                "CTSS", "FCN1", "FCGR3A", "S100A9", "S100A8", "FCGR3B")



####UMAP图###############################
#UMAP图
FeaturePlot(human_all, features = POD.marker, raster=FALSE, ncol = 4)
  ##得到1.UMAP图-POD.marker

FeaturePlot(human_all, features = PT.marker, raster=FALSE, ncol = 4)  
  ##得到2.UMAP图-PT.marker  

FeaturePlot(human_all, features = ATL.marker, raster=FALSE, ncol = 4)  
  ##得到3.UMAP图-ATL.marker  

FeaturePlot(human_all, features = DTL.marker, raster=FALSE, ncol = 4)  
  ##得到4.UMAP图-DTL.marker  

FeaturePlot(human_all, features = TAL.marker, raster=FALSE, ncol = 4)  
  ##得到5.UMAP图-TAL.marker  

FeaturePlot(human_all, features = DCT.marker, raster=FALSE, ncol = 4)  
  ##得到6.UMAP图-DCT.marker  

FeaturePlot(human_all, features = CNT.marker, raster=FALSE, ncol = 4)  
  ##得到7.UMAP图-CNT.marker  

FeaturePlot(human_all, features = PC.marker, raster=FALSE, ncol = 4)  
  ##得到8.UMAP图-PC.marker  

FeaturePlot(human_all, features = IC.marker, raster=FALSE, ncol = 4)  
  ##得到9.UMAP图-IC.marker  

FeaturePlot(human_all, features = EC.marker, raster=FALSE, ncol = 4)  
  ##得到10.UMAP图-EC.marker  

FeaturePlot(human_all, features = MC.marker, raster=FALSE, ncol = 4)  
  ##得到11.UMAP图-MC.marker  

FeaturePlot(human_all, features = VSMC.marker, raster=FALSE, ncol = 4)  
  ##得到12.UMAP图-VSMC.marker

FeaturePlot(human_all, features = FIB.marker, raster=FALSE, ncol = 4)  
  ##得到13.UMAP图-FIB.marker

FeaturePlot(human_all, features = IMM.marker, raster=FALSE, ncol = 4)  
  ##得到14.UMAP图-IMM.marker



####点图##################################
#点图
DotPlot(human_all, features = all.marker) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ##得到15.点图-全部marker



###小提琴图###########################
#小提琴图
levels(human_all) <- c("POD", "PT", "ATL", "DTL", "TAL", "DCT", "CNT", "PC",
                       "IC", "EC-PTC", "EC-GC", "EC-AEA", "EC-LYM", "cycEC",
                       "MC", "VSMC", "FIB", "NK", "B", "T", "MON", "MAC", 
                       "DC", "PL", "cycMNP", "MAST")
VlnPlot(human_all, 
        features = all.marker,
        stacked=T,
        pt.size=0,
        direction = "horizontal", #水平作图
        x.lab = '', y.lab = '')+#横纵轴不标记任何东西
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())#不显示坐标刻度
  ##得到16.小提琴图-全部marker



####皮尔逊检验##########################
#皮尔逊检验
table(human_all$usetype)
human_pearson <- AverageExpression(human_all,
                                   group.by = "usetype",
                                   assays = "RNA")
human_pearson <- human_pearson[[1]]
human_pearson <- as.matrix(human_pearson)
gene_pearson <- names(tail(sort(apply(human_pearson, 1, sd)),1000))
cell_pearson <- cor(human_pearson[gene_pearson,],method = 'spearman')
pheatmap::pheatmap(cell_pearson,
                   method = "spearman")
  ##得到17.全部细胞注释的spearman相关性
