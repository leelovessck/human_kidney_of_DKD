rm(list = ls())
gc()


#打开必要的package
{
if(!require(NMF))install.packages('NMF')
if(!require(circlize))devtools::install_github("jokergoo/circlize")
if(!require(ComplexHeatmap))devtools::install_github("jokergoo/ComplexHeatmap")
if(!require(CellChat))devtools::install_github("jinworks/CellChat")
if(!require(Seurat))install.packages("Seurat")
if(!require(patchwork))install.packages("patchwork")
if(!require(cowplot))install.packages("cowplot")
if(!require(dplyr))install.packages("dplyr")
}



####数据准备############################
#数据准备
setwd("D:/2023.10/调试/糖尿病肾病/1.人DKD分析")
getwd()
human_all <- readRDS("./data source/merge_rds/全部样本（自定义注释）.rds")
Idents(human_all) <- human_all$usetype
table(Idents(human_all))


usetype <- human_all@meta.data

usetype <- usetype %>%  
  mutate(usetype = case_when(  
    usetype %in% c("B", "PL") ~ "B",  
    TRUE ~ as.character(usetype)))

usetype <- usetype %>%  
  mutate(usetype = case_when(  
    usetype %in% c("MON", "MAC") ~ "MAC",  
    TRUE ~ as.character(usetype)))

usetype <- usetype %>%  
  mutate(usetype = case_when(  
    usetype %in% c("DC", "DC2") ~ "DC",  
    TRUE ~ as.character(usetype)))

usetype <- usetype %>%  
  mutate(usetype = case_when(  
    usetype %in% c("cycMNP", "MAST") ~ "Myeloid",  
    TRUE ~ as.character(usetype)))

usetype <- usetype %>%  
  mutate(usetype = case_when(  
    usetype %in% c("DTL", "ATL") ~ "LOH",  
    TRUE ~ as.character(usetype)))

usetype <- usetype %>%  
  mutate(usetype = case_when(  
    usetype %in% c("DCT", "TAL") ~ "DT",  
    TRUE ~ as.character(usetype)))

usetype <- usetype %>%  
  mutate(usetype = case_when(  
    usetype %in% c("CNT", "PC", "IC") ~ "CT",  
    TRUE ~ as.character(usetype)))

table(usetype$usetype)
human_all@meta.data <- usetype

human_all <- subset(human_all, subset = usetype != "EC-LYM")
human_all <- subset(human_all, subset = usetype != "cycEC")
human_all <- subset(human_all, subset = usetype != "NK")
human_all <- subset(human_all, subset = usetype != "Myeloid")
table(human_all$usetype)
unique(human_all$usetype)

saveRDS(human_all, 
        "D:/2023.10/调试/糖尿病肾病/1.人DKD分析/cellchat详细/cellchat的rds/全部细胞（seurat对象）.rds")
