rm(list = ls())
gc()
setwd("D:/2023.10/调试/糖尿病肾病/1.人DKD分析")
getwd()

#打开必要的package
{
if(!require(org.Hs.eg.db))BiocManager::install("org.Hs.eg.db")
if(!require(dplyr))install.packages("dplyr")
if(!require(clusterProfiler))BiocManager::install("clusterProfiler")
if(!require(ggplot2))install.packages("ggplot2")
if(!require(R.utils))install.packages("R.utils")
if(!require(tidyverse))install.packages("tidyverse")
if(!require(stringr))install.packages("stringr")
if(!require(enrichplot))install.packages("enrichplot")
if(!require(msigdbr))install.packages("msigdbr")
if(!require(GSVA))BiocManager::install("GSVA")
if(!require(pheatmap))install.packages("pheatmap")
if(!require(limma))BiocManager::install("limma")
if(!require(BiocParallel))install.packages("BiocParallel")
  if(!require(ReactomePA))BiocManager::install("ReactomePA")
}



####修改下载协议#####################
#修改下载协议
R.utils::setOption("clusterProfiler.download.method","auto")



####载入数据###############
#载入数据
ALL.DEG <- read.csv("./result/4.2.1 EC部位亚群DEG/全部EC的DEG.csv")
GC.DEG <- read.csv("./result/4.2.1 EC部位亚群DEG/EC-GC的DEG.csv")
AEA.DEG <- read.csv("./result/4.2.1 EC部位亚群DEG/EC-AEA的DEG.csv")
PTC.DEG <- read.csv("./result/4.2.1 EC部位亚群DEG/EC-PTC的DEG.csv")



####全部EC（all）####################
#全部EC（all）
rownames(ALL.DEG) <- ALL.DEG[,1]
colnames(ALL.DEG)[1] <- "SYMBOL"
ALL.DEG_rich <- ALL.DEG %>%     
  filter(p_val_adj < 0.05) %>%  
  mutate(abs_log2FC = abs(avg_log2FC)) %>%  
  arrange(desc(abs_log2FC)) %>% 
  slice(1:750) %>% 
  rownames_to_column("gene_id") %>%
  select(gene_id) %>%
  unlist() 

ALL.DEG_df <- bitr(ALL.DEG_rich,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL"),
               OrgDb = org.Hs.eg.db)
ALL.DEG_df <- ALL.DEG_df %>% distinct(SYMBOL, .keep_all = T)


##KEGG富集分析
kegg <- enrichKEGG(unique(ALL.DEG_df$ENTREZID), organism='hsa',
                          pvalueCutoff=0.05,pAdjustMethod='BH',qvalueCutoff=0.2,
                          minGSSize=10,maxGSSize=500,use_internal_data=F)
kegg <- setReadable(kegg,"org.Hs.eg.db","ENTREZID")  
  ##把ENTREZID转换为SYMBOL
write.csv(kegg@result, "1.1 全部EC的KEGG（all）.csv")


##GO（molecular function）
goMF <- enrichGO(ALL.DEG_df$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", ont = "MF",
                 pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,
                 maxGSSize = 500, readable = FALSE, pool = FALSE)
write.csv(goMF@result, "1.2 全部EC的GO_MF（all）.csv")


#GO（cell component）
goCC <- enrichGO(ALL.DEG_df$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", ont = "CC",
                 pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,
                 maxGSSize = 500, readable = FALSE, pool = FALSE)
write.csv(goCC@result, "1.3 全部EC的GO_CC（all）.csv")


#GO（biological process）
goBP <- enrichGO(ALL.DEG_df$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                 pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,
                 maxGSSize = 500, readable = FALSE, pool = FALSE)
write.csv(goBP@result, "1.4 全部EC的GO_BP（all）.csv")


#Reactome
Reactome <- enrichPathway(unique(ALL.DEG_df$ENTREZID),
                          organism = "human",
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.2,
                          minGSSize = 10,
                          maxGSSize = 500,
                          readable = FALSE)
Reactome <- setReadable(Reactome,"org.Hs.eg.db","ENTREZID")  
##把ENTREZID转换为SYMBOL
write.csv(Reactome@result, "1.5 全部EC的Reactome（all）.csv")



####EC-AEA（all）####################
#EC-AEA（all）
rownames(AEA.DEG) <- AEA.DEG[,1]
colnames(AEA.DEG)[1] <- "SYMBOL"
AEA.DEG_rich <- AEA.DEG %>%   
  filter(p_val_adj < 0.05) %>%  
  mutate(abs_log2FC = abs(avg_log2FC)) %>%  
  arrange(desc(abs_log2FC)) %>% 
  slice(1:2000) %>% 
  rownames_to_column("gene_id") %>%
  select(gene_id) %>%
  unlist() 

AEA.DEG_df <- bitr(AEA.DEG_rich,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL"),
                   OrgDb = org.Hs.eg.db)
AEA.DEG_df <- AEA.DEG_df %>% distinct(SYMBOL, .keep_all = T)


##KEGG富集分析
kegg <- enrichKEGG(unique(AEA.DEG_df$ENTREZID), organism='hsa',
                   pvalueCutoff=0.05,pAdjustMethod='BH',qvalueCutoff=0.2,
                   minGSSize=10,maxGSSize=500,use_internal_data=F)
kegg <- setReadable(kegg,"org.Hs.eg.db","ENTREZID")  
##把ENTREZID转换为SYMBOL
write.csv(kegg@result, "2.1 EC-AEA的KEGG（all）.csv")


##GO（molecular function）
goMF <- enrichGO(AEA.DEG_df$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", ont = "MF",
                 pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,
                 maxGSSize = 500, readable = FALSE, pool = FALSE)
write.csv(goMF@result, "2.2 EC-AEA的GO_MF（all）.csv")


#GO（cell component）
goCC <- enrichGO(AEA.DEG_df$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", ont = "CC",
                 pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,
                 maxGSSize = 500, readable = FALSE, pool = FALSE)
write.csv(goCC@result, "2.3 EC-AEA的GO_CC（all）.csv")


#GO（biological process）
goBP <- enrichGO(AEA.DEG_df$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                 pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,
                 maxGSSize = 500, readable = FALSE, pool = FALSE)
write.csv(goBP@result, "2.4 EC-AEA的GO_BP（all）.csv")

#Reactome
Reactome <- enrichPathway(unique(AEA.DEG_df$ENTREZID),
                          organism = "human",
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.2,
                          minGSSize = 10,
                          maxGSSize = 500,
                          readable = FALSE)
Reactome <- setReadable(Reactome,"org.Hs.eg.db","ENTREZID")  
##把ENTREZID转换为SYMBOL
write.csv(Reactome@result, "2.5 EC-AEA的Reactome（all）.csv")



####EC-GC（all）####################
#EC-GC（all）
rownames(GC.DEG) <- GC.DEG[,1]
colnames(GC.DEG)[1] <- "SYMBOL"
GC.DEG_rich <- GC.DEG %>%   
  filter(p_val_adj < 0.05) %>%  
  mutate(abs_log2FC = abs(avg_log2FC)) %>%  
  arrange(desc(abs_log2FC)) %>% 
  slice(1:1000) %>% 
  rownames_to_column("gene_id") %>%
  select(gene_id) %>%
  unlist() 

GC.DEG_df <- bitr(GC.DEG_rich,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL"),
                   OrgDb = org.Hs.eg.db)
GC.DEG_df <- GC.DEG_df %>% distinct(SYMBOL, .keep_all = T)


##KEGG富集分析
kegg <- enrichKEGG(unique(GC.DEG_df$ENTREZID), organism='hsa',
                   pvalueCutoff=0.05,pAdjustMethod='BH',qvalueCutoff=0.2,
                   minGSSize=10,maxGSSize=500,use_internal_data=F)
kegg <- setReadable(kegg,"org.Hs.eg.db","ENTREZID")  
##把ENTREZID转换为SYMBOL
write.csv(kegg@result, "3.1 EC-GC的KEGG（all）.csv")


##GO（molecular function）
goMF <- enrichGO(GC.DEG_df$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", ont = "MF",
                 pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,
                 maxGSSize = 500, readable = FALSE, pool = FALSE)
write.csv(goMF@result, "3.2 EC-GC的GO_MF（all）.csv")


#GO（cell component）
goCC <- enrichGO(GC.DEG_df$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", ont = "CC",
                 pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,
                 maxGSSize = 500, readable = FALSE, pool = FALSE)
write.csv(goCC@result, "3.3 EC-GC的GO_CC（all）.csv")


#GO（biological process）
goBP <- enrichGO(GC.DEG_df$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                 pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,
                 maxGSSize = 500, readable = FALSE, pool = FALSE)
write.csv(goBP@result, "3.4 EC-GC的GO_BP（all）.csv")


#Reactome
Reactome <- enrichPathway(unique(GC.DEG_df$ENTREZID),
                          organism = "human",
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.2,
                          minGSSize = 10,
                          maxGSSize = 500,
                          readable = FALSE)
Reactome <- setReadable(Reactome,"org.Hs.eg.db","ENTREZID")  
##把ENTREZID转换为SYMBOL
write.csv(Reactome@result, "3.5 EC-GC的Reactome（all）.csv")



####EC-PTC（all）####################
#EC-PTC（all）
rownames(PTC.DEG) <- PTC.DEG[,1]
colnames(PTC.DEG)[1] <- "SYMBOL"
PTC.DEG_rich <- PTC.DEG %>%   
  filter(p_val_adj < 0.05) %>%  
  mutate(abs_log2FC = abs(avg_log2FC)) %>%  
  arrange(desc(abs_log2FC)) %>% 
  slice(1:1000) %>% 
  rownames_to_column("gene_id") %>%
  select(gene_id) %>%
  unlist() 

PTC.DEG_df <- bitr(PTC.DEG_rich,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL"),
                  OrgDb = org.Hs.eg.db)
PTC.DEG_df <- PTC.DEG_df %>% distinct(SYMBOL, .keep_all = T)


##KEGG富集分析
kegg <- enrichKEGG(unique(PTC.DEG_df$ENTREZID), organism='hsa',
                   pvalueCutoff=0.05,pAdjustMethod='BH',qvalueCutoff=0.2,
                   minGSSize=10,maxGSSize=500,use_internal_data=F)
kegg <- setReadable(kegg,"org.Hs.eg.db","ENTREZID")  
##把ENTREZID转换为SYMBOL
write.csv(kegg@result, "4.1 EC-PTC的KEGG（all）.csv")


##GO（molecular function）
goMF <- enrichGO(PTC.DEG_df$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", ont = "MF",
                 pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,
                 maxGSSize = 500, readable = FALSE, pool = FALSE)
write.csv(goMF@result, "4.2 EC-PTC的GO_MF（all）.csv")


#GO（cell component）
goCC <- enrichGO(PTC.DEG_df$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", ont = "CC",
                 pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,
                 maxGSSize = 500, readable = FALSE, pool = FALSE)
write.csv(goCC@result, "4.3 EC-PTC的GO_CC（all）.csv")


#GO（biological process）
goBP <- enrichGO(PTC.DEG_df$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                 pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,
                 maxGSSize = 500, readable = FALSE, pool = FALSE)
write.csv(goBP@result, "4.4 EC-PTC的GO_BP（all）.csv")


#Reactome
Reactome <- enrichPathway(unique(PTC.DEG_df$ENTREZID),
                          organism = "human",
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.2,
                          minGSSize = 10,
                          maxGSSize = 500,
                          readable = FALSE)
Reactome <- setReadable(Reactome,"org.Hs.eg.db","ENTREZID")  
##把ENTREZID转换为SYMBOL
write.csv(Reactome@result, "4.5 EC-PTC的Reactome（all）.csv")



####全部EC（up）####################
#全部EC（up）
rownames(ALL.DEG) <- ALL.DEG[,1]
colnames(ALL.DEG)[1] <- "SYMBOL"
ALL.DEG_rich <- ALL.DEG %>%   
  filter(p_val_adj < 0.05) %>%  
  filter(avg_log2FC > 0) %>%
  mutate(abs_log2FC = abs(avg_log2FC)) %>%  
  arrange(desc(abs_log2FC)) %>% 
  slice(1:1000) %>% 
  rownames_to_column("gene_id") %>%
  select(gene_id) %>%
  unlist() 

ALL.DEG_df <- bitr(ALL.DEG_rich,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL"),
                   OrgDb = org.Hs.eg.db)
ALL.DEG_df <- ALL.DEG_df %>% distinct(SYMBOL, .keep_all = T)


##KEGG富集分析
kegg <- enrichKEGG(unique(ALL.DEG_df$ENTREZID), organism='hsa',
                   pvalueCutoff=0.05,pAdjustMethod='BH',qvalueCutoff=0.2,
                   minGSSize=10,maxGSSize=500,use_internal_data=F)
kegg <- setReadable(kegg,"org.Hs.eg.db","ENTREZID")  
##把ENTREZID转换为SYMBOL
write.csv(kegg@result, "1.6 全部EC的KEGG（up）.csv")


##GO（molecular function）
goMF <- enrichGO(ALL.DEG_df$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", ont = "MF",
                 pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,
                 maxGSSize = 500, readable = FALSE, pool = FALSE)
write.csv(goMF@result, "1.7 全部EC的GO_MF（up）.csv")


#GO（cell component）
goCC <- enrichGO(ALL.DEG_df$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", ont = "CC",
                 pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,
                 maxGSSize = 500, readable = FALSE, pool = FALSE)
write.csv(goCC@result, "1.8 全部EC的GO_CC（up）.csv")


#GO（biological process）
goBP <- enrichGO(ALL.DEG_df$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                 pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,
                 maxGSSize = 500, readable = FALSE, pool = FALSE)
write.csv(goBP@result, "1.9 全部EC的GO_BP（up）.csv")


#Reactome
Reactome <- enrichPathway(unique(ALL.DEG_df$ENTREZID),
                          organism = "human",
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.2,
                          minGSSize = 10,
                          maxGSSize = 500,
                          readable = FALSE)
Reactome <- setReadable(Reactome,"org.Hs.eg.db","ENTREZID")  
##把ENTREZID转换为SYMBOL
write.csv(Reactome@result, "1.10 全部EC的Reactome（up）.csv")



####EC-AEA（up）####################
#EC-AEA（up）
rownames(AEA.DEG) <- AEA.DEG[,1]
colnames(AEA.DEG)[1] <- "SYMBOL"
AEA.DEG_rich <- AEA.DEG %>%   
  filter(p_val_adj < 0.05) %>%  
  filter(avg_log2FC > 0) %>%
  mutate(abs_log2FC = abs(avg_log2FC)) %>%  
  arrange(desc(abs_log2FC)) %>% 
  slice(1:1000) %>% 
  rownames_to_column("gene_id") %>%
  select(gene_id) %>%
  unlist() 

AEA.DEG_df <- bitr(AEA.DEG_rich,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL"),
                   OrgDb = org.Hs.eg.db)
AEA.DEG_df <- AEA.DEG_df %>% distinct(SYMBOL, .keep_all = T)


##KEGG富集分析
kegg <- enrichKEGG(unique(AEA.DEG_df$ENTREZID), organism='hsa',
                   pvalueCutoff=0.05,pAdjustMethod='BH',qvalueCutoff=0.2,
                   minGSSize=10,maxGSSize=500,use_internal_data=F)
kegg <- setReadable(kegg,"org.Hs.eg.db","ENTREZID")  
##把ENTREZID转换为SYMBOL
write.csv(kegg@result, "2.6 EC-AEA的KEGG（up）.csv")


##GO（molecular function）
goMF <- enrichGO(AEA.DEG_df$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", ont = "MF",
                 pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,
                 maxGSSize = 500, readable = FALSE, pool = FALSE)
write.csv(goMF@result, "2.7 EC-AEA的GO_MF（up）.csv")


#GO（cell component）
goCC <- enrichGO(AEA.DEG_df$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", ont = "CC",
                 pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,
                 maxGSSize = 500, readable = FALSE, pool = FALSE)
write.csv(goCC@result, "2.8 EC-AEA的GO_CC（up）.csv")


#GO（biological process）
goBP <- enrichGO(AEA.DEG_df$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                 pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,
                 maxGSSize = 500, readable = FALSE, pool = FALSE)
write.csv(goBP@result, "2.9 EC-AEA的GO_BP（up）.csv")


#Reactome
Reactome <- enrichPathway(unique(AEA.DEG_df$ENTREZID),
                          organism = "human",
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.2,
                          minGSSize = 10,
                          maxGSSize = 500,
                          readable = FALSE)
Reactome <- setReadable(Reactome,"org.Hs.eg.db","ENTREZID")  
##把ENTREZID转换为SYMBOL
write.csv(Reactome@result, "2.10 EC-AEA的Reactome（up）.csv")



####EC-GC（up）####################
#EC-GC（up）
rownames(GC.DEG) <- GC.DEG[,1]
colnames(GC.DEG)[1] <- "SYMBOL"
GC.DEG_rich <- GC.DEG %>%   
  filter(p_val_adj < 0.05) %>%  
  filter(avg_log2FC > 0) %>%
  mutate(abs_log2FC = abs(avg_log2FC)) %>%  
  arrange(desc(abs_log2FC)) %>% 
  slice(1:1000) %>% 
  rownames_to_column("gene_id") %>%
  select(gene_id) %>%
  unlist() 

GC.DEG_df <- bitr(GC.DEG_rich,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL"),
                  OrgDb = org.Hs.eg.db)
GC.DEG_df <- GC.DEG_df %>% distinct(SYMBOL, .keep_all = T)


##KEGG富集分析
kegg <- enrichKEGG(unique(GC.DEG_df$ENTREZID), organism='hsa',
                   pvalueCutoff=0.05,pAdjustMethod='BH',qvalueCutoff=0.2,
                   minGSSize=10,maxGSSize=500,use_internal_data=F)
kegg <- setReadable(kegg,"org.Hs.eg.db","ENTREZID")  
##把ENTREZID转换为SYMBOL
write.csv(kegg@result, "3.6 EC-GC的KEGG（up）.csv")


##GO（molecular function）
goMF <- enrichGO(GC.DEG_df$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", ont = "MF",
                 pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,
                 maxGSSize = 500, readable = FALSE, pool = FALSE)
write.csv(goMF@result, "3.7 EC-GC的GO_MF（up）.csv")


#GO（cell component）
goCC <- enrichGO(GC.DEG_df$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", ont = "CC",
                 pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,
                 maxGSSize = 500, readable = FALSE, pool = FALSE)
write.csv(goCC@result, "3.8 EC-GC的GO_CC（up）.csv")


#GO（biological process）
goBP <- enrichGO(GC.DEG_df$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                 pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,
                 maxGSSize = 500, readable = FALSE, pool = FALSE)
write.csv(goBP@result, "3.9 EC-GC的GO_BP（up）.csv")


#Reactome
Reactome <- enrichPathway(unique(GC.DEG_df$ENTREZID),
                          organism = "human",
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.2,
                          minGSSize = 10,
                          maxGSSize = 500,
                          readable = FALSE)
Reactome <- setReadable(Reactome,"org.Hs.eg.db","ENTREZID")  
##把ENTREZID转换为SYMBOL
write.csv(Reactome@result, "3.10 EC-GC的Reactome（up）.csv")



####EC-PTC（up）####################
#EC-PTC（up）
rownames(PTC.DEG) <- PTC.DEG[,1]
colnames(PTC.DEG)[1] <- "SYMBOL"
PTC.DEG_rich <- PTC.DEG %>%   
  filter(p_val_adj < 0.05) %>%  
  filter(avg_log2FC > 0) %>%
  mutate(abs_log2FC = abs(avg_log2FC)) %>%  
  arrange(desc(abs_log2FC)) %>% 
  slice(1:1000) %>% 
  rownames_to_column("gene_id") %>%
  select(gene_id) %>%
  unlist() 

PTC.DEG_df <- bitr(PTC.DEG_rich,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL"),
                   OrgDb = org.Hs.eg.db)
PTC.DEG_df <- PTC.DEG_df %>% distinct(SYMBOL, .keep_all = T)


##KEGG富集分析
kegg <- enrichKEGG(unique(PTC.DEG_df$ENTREZID), organism='hsa',
                   pvalueCutoff=0.05,pAdjustMethod='BH',qvalueCutoff=0.2,
                   minGSSize=10,maxGSSize=500,use_internal_data=F)
kegg <- setReadable(kegg,"org.Hs.eg.db","ENTREZID")  
##把ENTREZID转换为SYMBOL
write.csv(kegg@result, "4.6 EC-PTC的KEGG（up）.csv")


##GO（molecular function）
goMF <- enrichGO(PTC.DEG_df$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", ont = "MF",
                 pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,
                 maxGSSize = 500, readable = FALSE, pool = FALSE)
write.csv(goMF@result, "4.7 EC-PTC的GO_MF（up）.csv")


#GO（cell component）
goCC <- enrichGO(PTC.DEG_df$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", ont = "CC",
                 pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,
                 maxGSSize = 500, readable = FALSE, pool = FALSE)
write.csv(goCC@result, "4.8 EC-PTC的GO_CC（up）.csv")


#GO（biological process）
goBP <- enrichGO(PTC.DEG_df$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                 pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,
                 maxGSSize = 500, readable = FALSE, pool = FALSE)
write.csv(goBP@result, "4.9 EC-PTC的GO_BP（up）.csv")


#Reactome
Reactome <- enrichPathway(unique(PTC.DEG_df$ENTREZID),
                          organism = "human",
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.2,
                          minGSSize = 10,
                          maxGSSize = 500,
                          readable = FALSE)
Reactome <- setReadable(Reactome,"org.Hs.eg.db","ENTREZID")  
##把ENTREZID转换为SYMBOL
write.csv(Reactome@result, "4.10 EC-PTC的Reactome（up）.csv")



####全部EC（down）####################
#全部EC（down）
rownames(ALL.DEG) <- ALL.DEG[,1]
colnames(ALL.DEG)[1] <- "SYMBOL"
ALL.DEG_rich <- ALL.DEG %>%   
  filter(p_val_adj < 0.05) %>%  
  filter(avg_log2FC < 0) %>%
  mutate(abs_log2FC = abs(avg_log2FC)) %>%  
  arrange(desc(abs_log2FC)) %>% 
  slice(1:1000) %>% 
  rownames_to_column("gene_id") %>%
  select(gene_id) %>%
  unlist() 

ALL.DEG_df <- bitr(ALL.DEG_rich,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL"),
                   OrgDb = org.Hs.eg.db)
ALL.DEG_df <- ALL.DEG_df %>% distinct(SYMBOL, .keep_all = T)


##KEGG富集分析
kegg <- enrichKEGG(unique(ALL.DEG_df$ENTREZID), organism='hsa',
                   pvalueCutoff=0.05,pAdjustMethod='BH',qvalueCutoff=0.2,
                   minGSSize=10,maxGSSize=500,use_internal_data=F)
kegg <- setReadable(kegg,"org.Hs.eg.db","ENTREZID")  
##把ENTREZID转换为SYMBOL
write.csv(kegg@result, "1.11 全部EC的KEGG（down）.csv")


##GO（molecular function）
goMF <- enrichGO(ALL.DEG_df$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", ont = "MF",
                 pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,
                 maxGSSize = 500, readable = FALSE, pool = FALSE)
write.csv(goMF@result, "1.12 全部EC的GO_MF（down）.csv")


#GO（cell component）
goCC <- enrichGO(ALL.DEG_df$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", ont = "CC",
                 pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,
                 maxGSSize = 500, readable = FALSE, pool = FALSE)
write.csv(goCC@result, "1.13 全部EC的GO_CC（down）.csv")


#GO（biological process）
goBP <- enrichGO(ALL.DEG_df$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                 pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,
                 maxGSSize = 500, readable = FALSE, pool = FALSE)
write.csv(goBP@result, "1.14 全部EC的GO_BP（down）.csv")


#Reactome
Reactome <- enrichPathway(unique(ALL.DEG_df$ENTREZID),
                          organism = "human",
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.2,
                          minGSSize = 10,
                          maxGSSize = 500,
                          readable = FALSE)
Reactome <- setReadable(Reactome,"org.Hs.eg.db","ENTREZID")  
##把ENTREZID转换为SYMBOL
write.csv(Reactome@result, "1.15 全部EC的Reactome（down）.csv")



####EC-AEA（down）####################
#EC-AEA（down）
rownames(AEA.DEG) <- AEA.DEG[,1]
colnames(AEA.DEG)[1] <- "SYMBOL"
AEA.DEG_rich <- AEA.DEG %>%   
  filter(p_val_adj < 0.05) %>%  
  filter(avg_log2FC < 0) %>%
  mutate(abs_log2FC = abs(avg_log2FC)) %>%  
  arrange(desc(abs_log2FC)) %>% 
  slice(1:1000) %>% 
  rownames_to_column("gene_id") %>%
  select(gene_id) %>%
  unlist() 

AEA.DEG_df <- bitr(AEA.DEG_rich,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL"),
                   OrgDb = org.Hs.eg.db)
AEA.DEG_df <- AEA.DEG_df %>% distinct(SYMBOL, .keep_all = T)


##KEGG富集分析
kegg <- enrichKEGG(unique(AEA.DEG_df$ENTREZID), organism='hsa',
                   pvalueCutoff=0.05,pAdjustMethod='BH',qvalueCutoff=0.2,
                   minGSSize=10,maxGSSize=500,use_internal_data=F)
kegg <- setReadable(kegg,"org.Hs.eg.db","ENTREZID")  
##把ENTREZID转换为SYMBOL
write.csv(kegg@result, "2.11 EC-AEA的KEGG（down）.csv")


##GO（molecular function）
goMF <- enrichGO(AEA.DEG_df$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", ont = "MF",
                 pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,
                 maxGSSize = 500, readable = FALSE, pool = FALSE)
write.csv(goMF@result, "2.12 EC-AEA的GO_MF（down）.csv")


#GO（cell component）
goCC <- enrichGO(AEA.DEG_df$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", ont = "CC",
                 pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,
                 maxGSSize = 500, readable = FALSE, pool = FALSE)
write.csv(goCC@result, "2.13 EC-AEA的GO_CC（down）.csv")


#GO（biological process）
goBP <- enrichGO(AEA.DEG_df$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                 pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,
                 maxGSSize = 500, readable = FALSE, pool = FALSE)
write.csv(goBP@result, "2.14 EC-AEA的GO_BP（down）.csv")


#Reactome
Reactome <- enrichPathway(unique(AEA.DEG_df$ENTREZID),
                          organism = "human",
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.2,
                          minGSSize = 10,
                          maxGSSize = 500,
                          readable = FALSE)
Reactome <- setReadable(Reactome,"org.Hs.eg.db","ENTREZID")  
##把ENTREZID转换为SYMBOL
write.csv(Reactome@result, "2.15 EC-AEA的Reactome（down）.csv")



####EC-GC（down）####################
#EC-GC（down）
rownames(GC.DEG) <- GC.DEG[,1]
colnames(GC.DEG)[1] <- "SYMBOL"
GC.DEG_rich <- GC.DEG %>%   
  filter(p_val_adj < 0.05) %>%  
  filter(avg_log2FC < 0) %>%
  mutate(abs_log2FC = abs(avg_log2FC)) %>%  
  arrange(desc(abs_log2FC)) %>% 
  slice(1:1000) %>% 
  rownames_to_column("gene_id") %>%
  select(gene_id) %>%
  unlist() 

GC.DEG_df <- bitr(GC.DEG_rich,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL"),
                  OrgDb = org.Hs.eg.db)
GC.DEG_df <- GC.DEG_df %>% distinct(SYMBOL, .keep_all = T)


##KEGG富集分析
kegg <- enrichKEGG(unique(GC.DEG_df$ENTREZID), organism='hsa',
                   pvalueCutoff=0.05,pAdjustMethod='BH',qvalueCutoff=0.2,
                   minGSSize=10,maxGSSize=500,use_internal_data=F)
kegg <- setReadable(kegg,"org.Hs.eg.db","ENTREZID")  
##把ENTREZID转换为SYMBOL
write.csv(kegg@result, "3.11 EC-GC的KEGG（down）.csv")


##GO（molecular function）
goMF <- enrichGO(GC.DEG_df$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", ont = "MF",
                 pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,
                 maxGSSize = 500, readable = FALSE, pool = FALSE)
write.csv(goMF@result, "3.12 EC-GC的GO_MF（down）.csv")


#GO（cell component）
goCC <- enrichGO(GC.DEG_df$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", ont = "CC",
                 pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,
                 maxGSSize = 500, readable = FALSE, pool = FALSE)
write.csv(goCC@result, "3.13 EC-GC的GO_CC（down）.csv")


#GO（biological process）
goBP <- enrichGO(GC.DEG_df$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                 pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,
                 maxGSSize = 500, readable = FALSE, pool = FALSE)
write.csv(goBP@result, "3.14 EC-GC的GO_BP（down）.csv")


#Reactome
Reactome <- enrichPathway(unique(GC.DEG_df$ENTREZID),
                          organism = "human",
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.2,
                          minGSSize = 10,
                          maxGSSize = 500,
                          readable = FALSE)
Reactome <- setReadable(Reactome,"org.Hs.eg.db","ENTREZID")  
##把ENTREZID转换为SYMBOL
write.csv(Reactome@result, "3.15 EC-GC的Reactome（down）.csv")



####EC-PTC（down）####################
#EC-PTC（down）
rownames(PTC.DEG) <- PTC.DEG[,1]
colnames(PTC.DEG)[1] <- "SYMBOL"
PTC.DEG_rich <- PTC.DEG %>%   
  filter(p_val_adj < 0.05) %>%  
  filter(avg_log2FC < 0) %>%
  mutate(abs_log2FC = abs(avg_log2FC)) %>%  
  arrange(desc(abs_log2FC)) %>% 
  slice(1:1000) %>% 
  rownames_to_column("gene_id") %>%
  select(gene_id) %>%
  unlist() 

PTC.DEG_df <- bitr(PTC.DEG_rich,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL"),
                   OrgDb = org.Hs.eg.db)
PTC.DEG_df <- PTC.DEG_df %>% distinct(SYMBOL, .keep_all = T)


##KEGG富集分析
kegg <- enrichKEGG(unique(PTC.DEG_df$ENTREZID), organism='hsa',
                   pvalueCutoff=0.05,pAdjustMethod='BH',qvalueCutoff=0.2,
                   minGSSize=10,maxGSSize=500,use_internal_data=F)
kegg <- setReadable(kegg,"org.Hs.eg.db","ENTREZID")  
##把ENTREZID转换为SYMBOL
write.csv(kegg@result, "4.11 EC-PTC的KEGG（down）.csv")


##GO（molecular function）
goMF <- enrichGO(PTC.DEG_df$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", ont = "MF",
                 pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,
                 maxGSSize = 500, readable = FALSE, pool = FALSE)
write.csv(goMF@result, "4.12 EC-PTC的GO_MF（down）.csv")


#GO（cell component）
goCC <- enrichGO(PTC.DEG_df$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", ont = "CC",
                 pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,
                 maxGSSize = 500, readable = FALSE, pool = FALSE)
write.csv(goCC@result, "4.13 EC-PTC的GO_CC（down）.csv")


#GO（biological process）
goBP <- enrichGO(PTC.DEG_df$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                 pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,
                 maxGSSize = 500, readable = FALSE, pool = FALSE)
write.csv(goBP@result, "4.14 EC-PTC的GO_BP（down）.csv")


#Reactome
Reactome <- enrichPathway(unique(PTC.DEG_df$ENTREZID),
                          organism = "human",
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.2,
                          minGSSize = 10,
                          maxGSSize = 500,
                          readable = FALSE)
Reactome <- setReadable(Reactome,"org.Hs.eg.db","ENTREZID")  
##把ENTREZID转换为SYMBOL
write.csv(Reactome@result, "4.15 EC-PTC的Reactome（down）.csv")
