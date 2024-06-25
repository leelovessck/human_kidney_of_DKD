rm(list = ls())
gc()



#打开必要的package
{
  if(!require(Seurat))install.packages("Seurat")
  if(!require(EnhancedVolcano))BiocManager::install("EnhancedVolcano")
  if(!require(presto))devtools::install_github("immunogenomics/presto")
  if(!require(tidyverse))install.packages("tidyverse")
  if(!require(ggrepel))install.packages("ggrepel")
  if(!require(VennDiagram))install.packages("VennDiagram")
}



setwd("D:/2023.10/调试/糖尿病肾病/1.人DKD分析")
getwd()
allDEG <- read.csv("./result/4.2.2 EC单细胞DEG韦恩图/allDEG.csv")
upDEG <- read.csv("./result/4.2.2 EC单细胞DEG韦恩图/upDEG.csv")
downDEG <- read.csv("./result/4.2.2 EC单细胞DEG韦恩图/downDEG.csv")



####全部DEG#######################
#全部DEG
venn_alllist <- list(all_EC = allDEG$all_EC, 
                     EC_AEA = allDEG$EC_AEA, 
                     EC_GC = allDEG$EC_GC,
                     EC_PTC = allDEG$EC_PTC)
venn.diagram(venn_alllist,
             filename = "1.全部DEG韦恩图.png", imagetype = "png", 
             fill = c("blue", "green", "orange", "grey"), alpha = 0.50, 
             cat.col = c("blue", "green", "orange", "grey"), cat.cex = 1.5, cat.fontfamily = "serif",
             col = c("blue", "green", "orange", "grey"), cex = 1.5, fontfamily = "serif")

inter <- get.venn.partitions(venn_alllist)
for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')
write.table(inter[-c(5, 6)], '内皮细胞全部DEG重合.txt', row.names = FALSE,
            sep = '\t', quote = FALSE)



####up DEG#######################
#up DEG
venn_uplist <- list(all_EC = upDEG$all_EC, 
                    EC_AEA = upDEG$EC_AEA, 
                    EC_GC = upDEG$EC_GC,
                    EC_PTC = upDEG$EC_PTC)
venn.diagram(venn_uplist,
             filename = "2.up DEG韦恩图.png", imagetype = "png", 
             fill = c("blue", "green", "orange", "grey"), alpha = 0.50, 
             cat.col = c("blue", "green", "orange", "grey"), cat.cex = 1.5, cat.fontfamily = "serif",
             col = c("blue", "green", "orange", "grey"), cex = 1.5, fontfamily = "serif")

inter <- get.venn.partitions(venn_uplist)
for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')
write.table(inter[-c(5, 6)], '内皮细胞up DEG重合.txt', row.names = FALSE,
            sep = '\t', quote = FALSE)



####down DEG#######################
#down DEG
venn_downlist <- list(all_EC = downDEG$all_EC, 
                      EC_AEA = downDEG$EC_AEA, 
                      EC_GC = downDEG$EC_GC,
                      EC_PTC = downDEG$EC_PTC)
venn.diagram(venn_downlist,
             filename = "3.down DEG韦恩图.png", imagetype = "png", 
             fill = c("blue", "green", "orange", "grey"), alpha = 0.50, 
             cat.col = c("blue", "green", "orange", "grey"), cat.cex = 1.5, cat.fontfamily = "serif",
             col = c("blue", "green", "orange", "grey"), cex = 1.5, fontfamily = "serif")

inter <- get.venn.partitions(venn_downlist)
for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')
write.table(inter[-c(5, 6)], '内皮细胞down DEG重合.txt', row.names = FALSE,
            sep = '\t', quote = FALSE)
