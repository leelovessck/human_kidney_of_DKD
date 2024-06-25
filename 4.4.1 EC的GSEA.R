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
  if(!require(limma))install.packages("limma")
  if(!require(BiocParallel))install.packages("BiocParallel")
  if(!require(ReactomePA))BiocManager::install("ReactomePA")
  if(!require(DOSE))install.packages("DOSE")
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



####全部EC####################
#全部EC
##数据准备
rownames(ALL.DEG) <- ALL.DEG[,1]
colnames(ALL.DEG)[1] <- "SYMBOL"
colnames(ALL.DEG)[3] <- "LogFC"
ALL.DEG <- ALL.DEG %>%   
  filter(p_val_adj < 0.05) %>%  
  select(SYMBOL, LogFC)
ALL.DEG_rich <- ALL.DEG$SYMBOL
entrezID <- bitr(ALL.DEG_rich,fromType = "SYMBOL", toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)
ALL.DEG <- inner_join(ALL.DEG, entrezID,
                      by = "SYMBOL")
ALL.DEG<-ALL.DEG[,-1]
ALL.DEG<-na.omit(ALL.DEG)
ALL.DEG <- ALL.DEG %>%   
  arrange(desc(LogFC))

genelist = ALL.DEG[["LogFC"]]
names(genelist) = as.character(ALL.DEG[["ENTREZID"]])


#运算过程
KEGG_GSEA <- gseKEGG(geneList = genelist,
                     organism = "hsa",
                     minGSSize = 10,
                     maxGSSize = 500,
                     pvalueCutoff = 1,
                     pAdjustMethod = "BH",
                     verbose = TRUE,
                     eps = 1e-10)
KEGG_GSEA_result <- KEGG_GSEA@result
write.csv(KEGG_GSEA_result, "1.1 全部EC的GSEA-KEGG.csv")

significant_indices <- which(KEGG_GSEA@result$pvalue < 0.05) 
for (i in significant_indices) {  
  current_description <- KEGG_GSEA@result$Description[i]  
  plot1 <- gseaplot(KEGG_GSEA, geneSetID = i,  
                    title = current_description,  
                    pvalue_table = TRUE)  
  ggsave(filename = paste0(current_description, ".pdf"),  
         plot = plot1[[1]]/plot1[[2]], height = 10, width = 10)  
}
  ##得到1.1 全部EC的GSEA-KEGG（各折线图）


GO_GSEA <- gseGO(genelist,
                 ont = "ALL",
                 OrgDb = 'org.Hs.eg.db',
                 keyType = "ENTREZID",
                 exponent = 1,
                 minGSSize = 10,
                 maxGSSize = 500,
                 eps = 1e-10,
                 pvalueCutoff = 1,
                 pAdjustMethod = "BH",
                 verbose = TRUE)
GO_GSEA_result <- GO_GSEA@result
write.csv(GO_GSEA_result, "1.2 全部EC的GSEA-GO.csv")

significant_indices <- which(GO_GSEA@result$pvalue < 0.05) 
for (i in significant_indices) {  
  current_description <- GO_GSEA@result$Description[i]  
  plot1 <- gseaplot(GO_GSEA, geneSetID = i,  
                    title = current_description,  
                    pvalue_table = TRUE)  
  ggsave(filename = paste0(current_description, ".pdf"),  
         plot = plot1[[1]]/plot1[[2]], height = 10, width = 10)  
}
  ##得到1.2 全部EC的GSEA-GO（各折线图）


Reactome_GSEA <- gsePathway(genelist,
                            organism = "human",
                            exponent = 1,
                            minGSSize = 10,
                            maxGSSize = 500,
                            eps = 1e-10,
                            pvalueCutoff = 1,
                            pAdjustMethod = "BH",
                            verbose = TRUE)
Reactome_GSEA_result <- Reactome_GSEA@result
write.csv(Reactome_GSEA_result, "1.3 全部EC的GSEA-Reactome.csv")

significant_indices <- which(Reactome_GSEA@result$pvalue < 0.05) 
for (i in significant_indices) {  
  current_description <- Reactome_GSEA@result$Description[i]  
  plot1 <- gseaplot(Reactome_GSEA, geneSetID = i,  
                    title = current_description,  
                    pvalue_table = TRUE)  
  ggsave(filename = paste0(current_description, ".pdf"),  
         plot = plot1[[1]]/plot1[[2]], height = 10, width = 10)  
}
  ##得到1.3 全部EC的GSEA-Reactome（各折线图）


DO_GSEA <- gseDO(genelist,
                 organism = "hsa",
                 exponent = 1,
                 minGSSize = 10,
                 maxGSSize = 500,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 verbose = TRUE)
DO_GSEA_result <- DO_GSEA@result
write.csv(DO_GSEA_result, "1.4 全部EC的GSEA-DO.csv")

significant_indices <- which(DO_GSEA@result$pvalue < 0.05) 
for (i in significant_indices) {  
  current_description <- DO_GSEA@result$Description[i]  
  plot1 <- gseaplot(DO_GSEA, geneSetID = i,  
                    title = current_description,  
                    pvalue_table = TRUE)  
  ggsave(filename = paste0(current_description, ".pdf"),  
         plot = plot1[[1]]/plot1[[2]], height = 10, width = 10)  
}
  ##得到1.4 全部EC的GSEA-DO（各折线图）


DGN_GSEA <- gseDGN(genelist,
                   exponent = 1,
                   minGSSize = 10,
                   maxGSSize = 500,
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "BH",
                   verbose = TRUE)
DGN_GSEA_result <- DGN_GSEA@result
write.csv(DGN_GSEA_result, "1.5 全部EC的GSEA-DisGeNET.csv")

significant_indices <- which(DGN_GSEA@result$pvalue < 0.05) 
for (i in significant_indices) {  
  current_description <- DGN_GSEA@result$Description[i]  
  plot1 <- gseaplot(DGN_GSEA, geneSetID = i,  
                    title = current_description,  
                    pvalue_table = TRUE)  
  ggsave(filename = paste0(current_description, ".pdf"),  
         plot = plot1[[1]]/plot1[[2]], height = 10, width = 10)  
}
  ##得到1.5 全部EC的GSEA-DisGeNET（各折线图）

hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)
H_GSEA <- GSEA(genelist,
               TERM2GENE = hallmark,
               exponent = 1,
               minGSSize = 10,
               maxGSSize = 500,
               eps = 1e-10,
               pvalueCutoff = 1,
               pAdjustMethod = "BH")
H_GSEA_result <- H_GSEA@result
write.csv(H_GSEA_result, "1.6 全部EC的GSEA-hallmark.csv")

significant_indices <- which(H_GSEA@result$pvalue < 0.05) 
for (i in significant_indices) {  
  current_description <- H_GSEA@result$Description[i]  
  plot1 <- gseaplot(H_GSEA, geneSetID = i,  
                    title = current_description,  
                    pvalue_table = TRUE)  
  ggsave(filename = paste0(current_description, ".pdf"),  
         plot = plot1[[1]]/plot1[[2]], height = 10, width = 10)  
}
  ##得到1.6 全部EC的GSEA-hallmark（各折线图）



####EC-AEA####################
#EC-AEA
##数据准备
rownames(AEA.DEG) <- AEA.DEG[,1]
colnames(AEA.DEG)[1] <- "SYMBOL"
colnames(AEA.DEG)[3] <- "LogFC"
AEA.DEG <- AEA.DEG %>%   
  filter(p_val_adj < 0.05) %>%  
  select(SYMBOL, LogFC)
AEA.DEG_rich <- AEA.DEG$SYMBOL
entrezID <- bitr(AEA.DEG_rich,fromType = "SYMBOL", toType = "ENTREZID",
                 OrgDb = org.Hs.eg.db)
AEA.DEG <- inner_join(AEA.DEG, entrezID,
                      by = "SYMBOL")
AEA.DEG<-AEA.DEG[,-1]
AEA.DEG<-na.omit(AEA.DEG)
AEA.DEG <- AEA.DEG %>%   
  arrange(desc(LogFC))

genelist = AEA.DEG[["LogFC"]]
names(genelist) = as.character(AEA.DEG[["ENTREZID"]])


#运算过程
KEGG_GSEA <- gseKEGG(geneList = genelist,
                     organism = "hsa",
                     minGSSize = 10,
                     maxGSSize = 500,
                     pvalueCutoff = 1,
                     pAdjustMethod = "BH",
                     verbose = TRUE,
                     eps = 1e-10)
KEGG_GSEA_result <- KEGG_GSEA@result
write.csv(KEGG_GSEA_result, "2.1 EC-AEA的GSEA-KEGG.csv")


GO_GSEA <- gseGO(genelist,
                 ont = "ALL",
                 OrgDb = 'org.Hs.eg.db',
                 keyType = "ENTREZID",
                 exponent = 1,
                 minGSSize = 10,
                 maxGSSize = 500,
                 eps = 1e-10,
                 pvalueCutoff = 1,
                 pAdjustMethod = "BH",
                 verbose = TRUE)
GO_GSEA_result <- GO_GSEA@result
write.csv(GO_GSEA_result, "2.2 EC-AEA的GSEA-GO.csv")

significant_indices <- which(GO_GSEA@result$pvalue < 0.05) 
for (i in significant_indices) {  
  current_description <- GO_GSEA@result$Description[i]  
  plot1 <- gseaplot(GO_GSEA, geneSetID = i,  
                    title = current_description,  
                    pvalue_table = TRUE)  
  ggsave(filename = paste0(current_description, ".pdf"),  
         plot = plot1[[1]]/plot1[[2]], height = 10, width = 10)  
}
  ##得到2.2 EC-AEA的GSEA-GO（各折线图）


Reactome_GSEA <- gsePathway(genelist,
                            organism = "human",
                            exponent = 1,
                            minGSSize = 10,
                            maxGSSize = 500,
                            eps = 1e-10,
                            pvalueCutoff = 1,
                            pAdjustMethod = "BH",
                            verbose = TRUE)
Reactome_GSEA_result <- Reactome_GSEA@result
write.csv(Reactome_GSEA_result, "2.3 EC-AEA的GSEA-Reactome.csv")

significant_indices <- which(Reactome_GSEA@result$pvalue < 0.05) 


DO_GSEA <- gseDO(genelist,
                 organism = "hsa",
                 exponent = 1,
                 minGSSize = 10,
                 maxGSSize = 500,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 verbose = TRUE)
DO_GSEA_result <- DO_GSEA@result
write.csv(DO_GSEA_result, "2.4 EC-AEA的GSEA-DO.csv")


DGN_GSEA <- gseDGN(genelist,
                   exponent = 1,
                   minGSSize = 10,
                   maxGSSize = 500,
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "BH",
                   verbose = TRUE)
DGN_GSEA_result <- DGN_GSEA@result
write.csv(DGN_GSEA_result, "2.5 EC-AEA的GSEA-DisGeNET.csv")

hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)
H_GSEA <- GSEA(genelist,
               TERM2GENE = hallmark,
               exponent = 1,
               minGSSize = 10,
               maxGSSize = 500,
               eps = 1e-10,
               pvalueCutoff = 1,
               pAdjustMethod = "BH")
H_GSEA_result <- H_GSEA@result
write.csv(H_GSEA_result, "2.6 EC-AEA的GSEA-hallmark.csv")

significant_indices <- which(H_GSEA@result$pvalue < 0.05) 
for (i in significant_indices) {  
  current_description <- H_GSEA@result$Description[i]  
  plot1 <- gseaplot(H_GSEA, geneSetID = i,  
                    title = current_description,  
                    pvalue_table = TRUE)  
  ggsave(filename = paste0(current_description, ".pdf"),  
         plot = plot1[[1]]/plot1[[2]], height = 10, width = 10)  
}
  ##得到2.6 EC-AEA的GSEA-hallmark（各折线图）



####EC-GC####################
#EC-GC
##数据准备
rownames(GC.DEG) <- GC.DEG[,1]
colnames(GC.DEG)[1] <- "SYMBOL"
colnames(GC.DEG)[3] <- "LogFC"
GC.DEG <- GC.DEG %>%   
  filter(p_val_adj < 0.05) %>%  
  select(SYMBOL, LogFC)
GC.DEG_rich <- GC.DEG$SYMBOL
entrezID <- bitr(GC.DEG_rich,fromType = "SYMBOL", toType = "ENTREZID",
                 OrgDb = org.Hs.eg.db)
GC.DEG <- inner_join(GC.DEG, entrezID,
                      by = "SYMBOL")
GC.DEG<-GC.DEG[,-1]
GC.DEG<-na.omit(GC.DEG)
GC.DEG <- GC.DEG %>%   
  arrange(desc(LogFC))

genelist = GC.DEG[["LogFC"]]
names(genelist) = as.character(GC.DEG[["ENTREZID"]])


#运算过程
KEGG_GSEA <- gseKEGG(geneList = genelist,
                     organism = "hsa",
                     minGSSize = 10,
                     maxGSSize = 500,
                     pvalueCutoff = 1,
                     pAdjustMethod = "BH",
                     verbose = TRUE,
                     eps = 1e-10)
KEGG_GSEA_result <- KEGG_GSEA@result
write.csv(KEGG_GSEA_result, "3.1 EC-GC的GSEA-KEGG.csv")

significant_indices <- which(KEGG_GSEA@result$pvalue < 0.05) 
for (i in significant_indices) {  
  current_description <- KEGG_GSEA@result$Description[i]  
  plot1 <- gseaplot(KEGG_GSEA, geneSetID = i,  
                    title = current_description,  
                    pvalue_table = TRUE)  
  ggsave(filename = paste0(current_description, ".pdf"),  
         plot = plot1[[1]]/plot1[[2]], height = 10, width = 10)  
}
  ##得到3.1 EC-GC的GSEA-KEGG（各折线图）


GO_GSEA <- gseGO(genelist,
                 ont = "ALL",
                 OrgDb = 'org.Hs.eg.db',
                 keyType = "ENTREZID",
                 exponent = 1,
                 minGSSize = 10,
                 maxGSSize = 500,
                 eps = 1e-10,
                 pvalueCutoff = 1,
                 pAdjustMethod = "BH",
                 verbose = TRUE)
GO_GSEA_result <- GO_GSEA@result
write.csv(GO_GSEA_result, "3.2 EC-GC的GSEA-GO.csv")

significant_indices <- which(GO_GSEA@result$pvalue < 0.05) 
for (i in significant_indices) {  
  current_description <- GO_GSEA@result$Description[i]  
  plot1 <- gseaplot(GO_GSEA, geneSetID = i,  
                    title = current_description,  
                    pvalue_table = TRUE)  
  ggsave(filename = paste0(current_description, ".pdf"),  
         plot = plot1[[1]]/plot1[[2]], height = 10, width = 10)  
}
  ##得到3.2 EC-GC的GSEA-GO（各折线图）


Reactome_GSEA <- gsePathway(genelist,
                            organism = "human",
                            exponent = 1,
                            minGSSize = 10,
                            maxGSSize = 500,
                            eps = 1e-10,
                            pvalueCutoff = 1,
                            pAdjustMethod = "BH",
                            verbose = TRUE)
Reactome_GSEA_result <- Reactome_GSEA@result
write.csv(Reactome_GSEA_result, "3.3 EC-GC的GSEA-Reactome.csv")

significant_indices <- which(Reactome_GSEA@result$pvalue < 0.05) 
for (i in significant_indices) {  
  current_description <- Reactome_GSEA@result$Description[i]  
  plot1 <- gseaplot(Reactome_GSEA, geneSetID = i,  
                    title = current_description,  
                    pvalue_table = TRUE)  
  ggsave(filename = paste0(current_description, ".pdf"),  
         plot = plot1[[1]]/plot1[[2]], height = 10, width = 10)  
}
  ##得到3.3 EC-GC的GSEA-Reactome（各折线图）


DO_GSEA <- gseDO(genelist,
                 organism = "hsa",
                 exponent = 1,
                 minGSSize = 10,
                 maxGSSize = 500,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 verbose = TRUE)
DO_GSEA_result <- DO_GSEA@result
write.csv(DO_GSEA_result, "3.4 EC-GC的GSEA-DO.csv")


DGN_GSEA <- gseDGN(genelist,
                   exponent = 1,
                   minGSSize = 10,
                   maxGSSize = 500,
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "BH",
                   verbose = TRUE)
DGN_GSEA_result <- DGN_GSEA@result
write.csv(DGN_GSEA_result, "3.5 EC-GC的GSEA-DisGeNET.csv")


hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)
H_GSEA <- GSEA(genelist,
               TERM2GENE = hallmark,
               exponent = 1,
               minGSSize = 10,
               maxGSSize = 500,
               eps = 1e-10,
               pvalueCutoff = 1,
               pAdjustMethod = "BH")
H_GSEA_result <- H_GSEA@result
write.csv(H_GSEA_result, "3.6 EC-GC的GSEA-hallmark.csv")

significant_indices <- which(H_GSEA@result$pvalue < 0.05) 
for (i in significant_indices) {  
  current_description <- H_GSEA@result$Description[i]  
  plot1 <- gseaplot(H_GSEA, geneSetID = i,  
                    title = current_description,  
                    pvalue_table = TRUE)  
  ggsave(filename = paste0(current_description, ".pdf"),  
         plot = plot1[[1]]/plot1[[2]], height = 10, width = 10)  
}
  ##得到3.6 EC-GC的GSEA-hallmark（各折线图）



####EC-PTC####################
#EC-PTC
##数据准备
rownames(PTC.DEG) <- PTC.DEG[,1]
colnames(PTC.DEG)[1] <- "SYMBOL"
colnames(PTC.DEG)[3] <- "LogFC"
PTC.DEG <- PTC.DEG %>%   
  filter(p_val_adj < 0.05) %>%  
  select(SYMBOL, LogFC)
PTC.DEG_rich <- PTC.DEG$SYMBOL
entrezID <- bitr(PTC.DEG_rich,fromType = "SYMBOL", toType = "ENTREZID",
                 OrgDb = org.Hs.eg.db)
PTC.DEG <- inner_join(PTC.DEG, entrezID,
                     by = "SYMBOL")
PTC.DEG<-PTC.DEG[,-1]
PTC.DEG<-na.omit(PTC.DEG)
PTC.DEG <- PTC.DEG %>%   
  arrange(desc(LogFC))

genelist = PTC.DEG[["LogFC"]]
names(genelist) = as.character(PTC.DEG[["ENTREZID"]])


#运算过程
KEGG_GSEA <- gseKEGG(geneList = genelist,
                     organism = "hsa",
                     minGSSize = 10,
                     maxGSSize = 500,
                     pvalueCutoff = 1,
                     pAdjustMethod = "BH",
                     verbose = TRUE,
                     eps = 1e-10)
KEGG_GSEA_result <- KEGG_GSEA@result
write.csv(KEGG_GSEA_result, "4.1 EC-PTC的GSEA-KEGG.csv")

significant_indices <- which(KEGG_GSEA@result$pvalue < 0.05) 
for (i in significant_indices) {  
  current_description <- KEGG_GSEA@result$Description[i]  
  plot1 <- gseaplot(KEGG_GSEA, geneSetID = i,  
                    title = current_description,  
                    pvalue_table = TRUE)  
  ggsave(filename = paste0(current_description, ".pdf"),  
         plot = plot1[[1]]/plot1[[2]], height = 10, width = 10)  
}
  ##得到4.1 EC-AEA的GSEA-KEGG（各折线图）


GO_GSEA <- gseGO(genelist,
                 ont = "ALL",
                 OrgDb = 'org.Hs.eg.db',
                 keyType = "ENTREZID",
                 exponent = 1,
                 minGSSize = 10,
                 maxGSSize = 500,
                 eps = 1e-10,
                 pvalueCutoff = 1,
                 pAdjustMethod = "BH",
                 verbose = TRUE)
GO_GSEA_result <- GO_GSEA@result
write.csv(GO_GSEA_result, "4.2 EC-PTC的GSEA-GO.csv")

significant_indices <- which(GO_GSEA@result$pvalue < 0.05) 
for (i in significant_indices) {  
  current_description <- GO_GSEA@result$Description[i]  
  plot1 <- gseaplot(GO_GSEA, geneSetID = i,  
                    title = current_description,  
                    pvalue_table = TRUE)  
  ggsave(filename = paste0(current_description, ".pdf"),  
         plot = plot1[[1]]/plot1[[2]], height = 10, width = 10)  
}
  ##得到4.2 EC-AEA的GSEA-GO（各折线图）


Reactome_GSEA <- gsePathway(genelist,
                            organism = "human",
                            exponent = 1,
                            minGSSize = 10,
                            maxGSSize = 500,
                            eps = 1e-10,
                            pvalueCutoff = 1,
                            pAdjustMethod = "BH",
                            verbose = TRUE)
Reactome_GSEA_result <- Reactome_GSEA@result
write.csv(Reactome_GSEA_result, "4.3 EC-PTC的GSEA-Reactome.csv")

significant_indices <- which(Reactome_GSEA@result$pvalue < 0.05) 
for (i in significant_indices) {  
  current_description <- Reactome_GSEA@result$Description[i]  
  plot1 <- gseaplot(Reactome_GSEA, geneSetID = i,  
                    title = current_description,  
                    pvalue_table = TRUE)  
  ggsave(filename = paste0(current_description, ".pdf"),  
         plot = plot1[[1]]/plot1[[2]], height = 10, width = 10)  
}
  ##得到4.3 EC-AEA的GSEA-Reactome（各折线图）


DO_GSEA <- gseDO(genelist,
                 organism = "hsa",
                 exponent = 1,
                 minGSSize = 10,
                 maxGSSize = 500,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 verbose = TRUE)
DO_GSEA_result <- DO_GSEA@result
write.csv(DO_GSEA_result, "4.4 EC-PTC的GSEA-DO.csv")

significant_indices <- which(DO_GSEA@result$pvalue < 0.05) 
for (i in significant_indices) {  
  current_description <- DO_GSEA@result$Description[i]  
  plot1 <- gseaplot(DO_GSEA, geneSetID = i,  
                    title = current_description,  
                    pvalue_table = TRUE)  
  ggsave(filename = paste0(current_description, ".pdf"),  
         plot = plot1[[1]]/plot1[[2]], height = 10, width = 10)  
}
  ##得到4.4 EC-AEA的GSEA-DO（各折线图）


DGN_GSEA <- gseDGN(genelist,
                   exponent = 1,
                   minGSSize = 10,
                   maxGSSize = 500,
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "BH",
                   verbose = TRUE)
DGN_GSEA_result <- DGN_GSEA@result
write.csv(DGN_GSEA_result, "4.5 EC-PTC的GSEA-DisGeNET.csv")

significant_indices <- which(DGN_GSEA@result$pvalue < 0.05) 
for (i in significant_indices) {  
  current_description <- DGN_GSEA@result$Description[i]  
  plot1 <- gseaplot(DGN_GSEA, geneSetID = i,  
                    title = current_description,  
                    pvalue_table = TRUE)  
  ggsave(filename = paste0(current_description, ".pdf"),  
         plot = plot1[[1]]/plot1[[2]], height = 10, width = 10)  
}
  ##得到4.5 EC-AEA的GSEA-DGN（各折线图）


hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)
H_GSEA <- GSEA(genelist,
               TERM2GENE = hallmark,
               exponent = 1,
               minGSSize = 10,
               maxGSSize = 500,
               eps = 1e-10,
               pvalueCutoff = 1,
               pAdjustMethod = "BH")
H_GSEA_result <- H_GSEA@result
write.csv(H_GSEA_result, "4.6 EC-PTC的GSEA-hallmark.csv")

significant_indices <- which(H_GSEA@result$pvalue < 0.05) 
for (i in significant_indices) {  
  current_description <- H_GSEA@result$Description[i]  
  plot1 <- gseaplot(H_GSEA, geneSetID = i,  
                    title = current_description,  
                    pvalue_table = TRUE)  
  ggsave(filename = paste0(current_description, ".pdf"),  
         plot = plot1[[1]]/plot1[[2]], height = 10, width = 10)  
}
  ##得到4.6 EC-PTC的GSEA-hallmark（各折线图）
