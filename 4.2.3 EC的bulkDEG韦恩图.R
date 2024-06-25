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
  if(!require(readxl))install.packages(readxl)
}



####载入数据#######################
#载入数据
rm(list = ls())
gc()
setwd("D:/2023.10/调试/糖尿病肾病/1.人DKD分析")
getwd()

allDEG <- read.csv("./result/4.2.3 EC的bulkDEG韦恩图/allDEG.csv")
upDEG <- read.csv("./result/4.2.3 EC的bulkDEG韦恩图/upDEG.csv")
downDEG <- read.csv("./result/4.2.3 EC的bulkDEG韦恩图/downDEG.csv")



####全部EC###########################
#全部EC
##all
venn_ECalllist <- list(all_EC = allDEG$all_EC,
                       group21 = allDEG$bulk_G12,
                       group32 = allDEG$bulk_G32)
venn.diagram(venn_ECalllist,
             filename = "1.1 ALL EC的全部DEG韦恩图.png", imagetype = "png", 
             fill = c("blue", "green", "grey"), alpha = 0.50, 
             cat.col = c("blue", "green", "grey"), cat.cex = 1.5, cat.fontfamily = "serif",
             col = c("blue", "green", "grey"), cex = 1.5, fontfamily = "serif")
inter <- get.venn.partitions(venn_ECalllist)
for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')
write.table(inter[-c(5, 6)], 'ALL EC的全部DEG重合.txt', row.names = FALSE,
            sep = '\t', quote = FALSE)


##up
venn_ECuplist <- list(all_EC = upDEG$all_EC,
                       group21 = upDEG$bulk_G12,
                       group32 = upDEG$bulk_G32)
venn.diagram(venn_ECuplist,
             filename = "1.2 ALL EC的up DEG韦恩图.png", imagetype = "png", 
             fill = c("blue", "green", "grey"), alpha = 0.50, 
             cat.col = c("blue", "green", "grey"), cat.cex = 1.5, cat.fontfamily = "serif",
             col = c("blue", "green", "grey"), cex = 1.5, fontfamily = "serif")
inter <- get.venn.partitions(venn_ECuplist)
for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')
write.table(inter[-c(5, 6)], 'ALL EC的up DEG重合.txt', row.names = FALSE,
            sep = '\t', quote = FALSE)


##down
venn_ECdownlist <- list(all_EC = downDEG$all_EC,
                      group21 = downDEG$bulk_G12,
                      group32 = downDEG$bulk_G32)
venn.diagram(venn_ECdownlist,
             filename = "1.3 ALL EC的down DEG韦恩图.png", imagetype = "png", 
             fill = c("blue", "green", "grey"), alpha = 0.50, 
             cat.col = c("blue", "green", "grey"), cat.cex = 1.5, cat.fontfamily = "serif",
             col = c("blue", "green", "grey"), cex = 1.5, fontfamily = "serif")
inter <- get.venn.partitions(venn_ECdownlist)
for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')
write.table(inter[-c(5, 6)], 'ALL EC的down DEG重合.txt', row.names = FALSE,
            sep = '\t', quote = FALSE)



####EC-AEA###########################
#EC-AEA
##all
venn_ECalllist <- list(EC_AEA = allDEG$EC_AEA,
                       group21 = allDEG$bulk_G12,
                       group32 = allDEG$bulk_G32)
venn.diagram(venn_ECalllist,
             filename = "2.1 EC-AEA的全部DEG韦恩图.png", imagetype = "png", 
             fill = c("blue", "green", "grey"), alpha = 0.50, 
             cat.col = c("blue", "green", "grey"), cat.cex = 1.5, cat.fontfamily = "serif",
             col = c("blue", "green", "grey"), cex = 1.5, fontfamily = "serif")
inter <- get.venn.partitions(venn_ECalllist)
for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')
write.table(inter[-c(5, 6)], 'EC-AEA的全部DEG重合.txt', row.names = FALSE,
            sep = '\t', quote = FALSE)


##up
venn_ECuplist <- list(EC_AEA = upDEG$EC_AEA,
                      group21 = upDEG$bulk_G12,
                      group32 = upDEG$bulk_G32)
venn.diagram(venn_ECuplist,
             filename = "2.2 EC-AEA的up DEG韦恩图.png", imagetype = "png", 
             fill = c("blue", "green", "grey"), alpha = 0.50, 
             cat.col = c("blue", "green", "grey"), cat.cex = 1.5, cat.fontfamily = "serif",
             col = c("blue", "green", "grey"), cex = 1.5, fontfamily = "serif")
inter <- get.venn.partitions(venn_ECuplist)
for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')
write.table(inter[-c(5, 6)], 'EC-AEA的up DEG重合.txt', row.names = FALSE,
            sep = '\t', quote = FALSE)


##down
venn_ECdownlist <- list(EC_AEA = downDEG$EC_AEA,
                        group21 = downDEG$bulk_G12,
                        group32 = downDEG$bulk_G32)
venn.diagram(venn_ECdownlist,
             filename = "2.3 EC_AEA的down DEG韦恩图.png", imagetype = "png", 
             fill = c("blue", "green", "grey"), alpha = 0.50, 
             cat.col = c("blue", "green", "grey"), cat.cex = 1.5, cat.fontfamily = "serif",
             col = c("blue", "green", "grey"), cex = 1.5, fontfamily = "serif")
inter <- get.venn.partitions(venn_ECdownlist)
for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')
write.table(inter[-c(5, 6)], 'EC_AEA的down DEG重合.txt', row.names = FALSE,
            sep = '\t', quote = FALSE)



####EC-PTC###########################
#EC-PTC
##all
venn_ECalllist <- list(EC_PTC = allDEG$EC_PTC,
                       group21 = allDEG$bulk_G12,
                       group32 = allDEG$bulk_G32)
venn.diagram(venn_ECalllist,
             filename = "3.1 EC-PTC的全部DEG韦恩图.png", imagetype = "png", 
             fill = c("blue", "green", "grey"), alpha = 0.50, 
             cat.col = c("blue", "green", "grey"), cat.cex = 1.5, cat.fontfamily = "serif",
             col = c("blue", "green", "grey"), cex = 1.5, fontfamily = "serif")
inter <- get.venn.partitions(venn_ECalllist)
for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')
write.table(inter[-c(5, 6)], 'EC-PTC的全部DEG重合.txt', row.names = FALSE,
            sep = '\t', quote = FALSE)


##up
venn_ECuplist <- list(EC_PTC = upDEG$EC_PTC,
                      group21 = upDEG$bulk_G12,
                      group32 = upDEG$bulk_G32)
venn.diagram(venn_ECuplist,
             filename = "3.2 EC-PTC的up DEG韦恩图.png", imagetype = "png", 
             fill = c("blue", "green", "grey"), alpha = 0.50, 
             cat.col = c("blue", "green", "grey"), cat.cex = 1.5, cat.fontfamily = "serif",
             col = c("blue", "green", "grey"), cex = 1.5, fontfamily = "serif")
inter <- get.venn.partitions(venn_ECuplist)
for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')
write.table(inter[-c(5, 6)], 'EC-PTC的up DEG重合.txt', row.names = FALSE,
            sep = '\t', quote = FALSE)


##down
venn_ECdownlist <- list(EC_PTC = downDEG$EC_PTC,
                        group21 = downDEG$bulk_G12,
                        group32 = downDEG$bulk_G32)
venn.diagram(venn_ECdownlist,
             filename = "3.3 EC_PTC的down DEG韦恩图.png", imagetype = "png", 
             fill = c("blue", "green", "grey"), alpha = 0.50, 
             cat.col = c("blue", "green", "grey"), cat.cex = 1.5, cat.fontfamily = "serif",
             col = c("blue", "green", "grey"), cex = 1.5, fontfamily = "serif")
inter <- get.venn.partitions(venn_ECdownlist)
for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')
write.table(inter[-c(5, 6)], 'EC_PTC的down DEG重合.txt', row.names = FALSE,
            sep = '\t', quote = FALSE)



####EC-GC###########################
#EC-GC
##all
venn_ECalllist <- list(EC_GC = allDEG$EC_GC,
                       group21 = allDEG$bulk_G12,
                       group32 = allDEG$bulk_G32)
venn.diagram(venn_ECalllist,
             filename = "4.1 EC-GC的全部DEG韦恩图.png", imagetype = "png", 
             fill = c("blue", "green", "grey"), alpha = 0.50, 
             cat.col = c("blue", "green", "grey"), cat.cex = 1.5, cat.fontfamily = "serif",
             col = c("blue", "green", "grey"), cex = 1.5, fontfamily = "serif")
inter <- get.venn.partitions(venn_ECalllist)
for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')
write.table(inter[-c(5, 6)], 'EC-GC的全部DEG重合.txt', row.names = FALSE,
            sep = '\t', quote = FALSE)


##up
venn_ECuplist <- list(EC_GC = upDEG$EC_GC,
                      group21 = upDEG$bulk_G12,
                      group32 = upDEG$bulk_G32)
venn.diagram(venn_ECuplist,
             filename = "4.2 EC-GC的up DEG韦恩图.png", imagetype = "png", 
             fill = c("blue", "green", "grey"), alpha = 0.50, 
             cat.col = c("blue", "green", "grey"), cat.cex = 1.5, cat.fontfamily = "serif",
             col = c("blue", "green", "grey"), cex = 1.5, fontfamily = "serif")
inter <- get.venn.partitions(venn_ECuplist)
for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')
write.table(inter[-c(5, 6)], 'EC-GC的up DEG重合.txt', row.names = FALSE,
            sep = '\t', quote = FALSE)


##down
venn_ECdownlist <- list(EC_GC = downDEG$EC_GC,
                        group21 = downDEG$bulk_G12,
                        group32 = downDEG$bulk_G32)
venn.diagram(venn_ECdownlist,
             filename = "4.3 EC_GC的down DEG韦恩图.png", imagetype = "png", 
             fill = c("blue", "green", "grey"), alpha = 0.50, 
             cat.col = c("blue", "green", "grey"), cat.cex = 1.5, cat.fontfamily = "serif",
             col = c("blue", "green", "grey"), cex = 1.5, fontfamily = "serif")
inter <- get.venn.partitions(venn_ECdownlist)
for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')
write.table(inter[-c(5, 6)], 'EC_GC的down DEG重合.txt', row.names = FALSE,
            sep = '\t', quote = FALSE)
