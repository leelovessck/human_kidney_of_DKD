rm(list = ls())
gc()

#打开必要的package
{
  if(!require(Seurat))install.packages("Seurat")
  if(!require(EnhancedVolcano))BiocManager::install("EnhancedVolcano")
  if(!require(presto))devtools::install_github("immunogenomics/presto")
  if(!require(tidyverse))install.packages("tidyverse")
  if(!require(ggrepel))install.packages("ggrepel")
}

setwd("D:/2023.10/调试/糖尿病肾病/1.人DKD分析")
getwd()
human_all <- readRDS("./data source/merge_rds/全部样本（自定义注释）.rds")



####提取基因（EC-GC）############################
#提取基因（EC-GC）
EC <- subset(human_all, usetype == "EC-GC")
EC$usetype <- droplevels(EC$usetype,
                         exclude = setdiff(
                           levels(EC$usetype),
                           unique(EC$usetype)))
table(EC$usetype)

CON_EC <- subset(EC, sampletype == "LD")
table(CON_EC$usetype)
DKD_EC <- subset(EC, sampletype == "DKD")
table(DKD_EC$usetype)

express_CONEC <- GetAssayData(CON_EC, layer = "data")
expression_CONEC <- as.data.frame(express_CONEC)
expression_CONEC <- expression_CONEC %>% 
  rownames_to_column("gene_id")
LIF_CONEC <- expression_CONEC %>% 
  filter(gene_id == "LIF" | gene_id == "LIFR" | gene_id == "LIFR-AS1")
write.csv(LIF_CONEC, "LIF_CONEC-GC.csv")

express_DKDEC <- GetAssayData(DKD_EC, layer = "data")
expression_DKDEC <- as.data.frame(express_DKDEC)
expression_DKDEC <- expression_DKDEC %>% 
  rownames_to_column("gene_id")
LIF_DKDEC <- expression_DKDEC %>% 
  filter(gene_id == "LIF" | gene_id == "LIFR" | gene_id == "LIFR-AS1")
write.csv(LIF_DKDEC, "LIF_DKDEC-GC.csv")



####读取长数据基因（EC-GC）###############
#读取长数据基因（EC-GC）
human <- read.csv("./result/8.1.1 LIF-基因/长-LIF_EC-GC.csv")

p1 <- ggplot(human,
             aes(x = group,
                 y = express,
                 color = group,
                 fill = group))+
  geom_bar(stat = "summary")+
  stat_summary(geom = "errorbar", width = 0.3)+
  theme_bw()+
  facet_wrap(~gene,scales = "free")
p1
filename <- "1.1 EC-GC的Lif表达量总览.png"  
png(filename, width = 1000, height = 500)  
print(p1)  
dev.off()

p1 <- ggplot(human,
             aes(x = group,
                 y = express,
                 color = group,
                 fill = group))+
  geom_bar(stat = "summary")+
  stat_summary(geom = "errorbar", width = 0.3)+
  theme_bw()+
  facet_wrap(~gene,scales = "fixed")
p1
filename <- "1.2 EC-GC的Lif表达量总览（固定刻度）.png"  
png(filename, width = 1000, height = 500)  
print(p1)  
dev.off()



####提取基因（EC-AEA）############################
#提取基因（EC-AEA）
rm(list = ls()[!ls() %in% "human_all"])

EC <- subset(human_all, usetype == "EC-AEA")
EC$usetype <- droplevels(EC$usetype,
                         exclude = setdiff(
                           levels(EC$usetype),
                           unique(EC$usetype)))
table(EC$usetype)

CON_EC <- subset(EC, sampletype == "LD")
table(CON_EC$usetype)
DKD_EC <- subset(EC, sampletype == "DKD")
table(DKD_EC$usetype)

express_CONEC <- GetAssayData(CON_EC, layer = "data")
expression_CONEC <- as.data.frame(express_CONEC)
expression_CONEC <- expression_CONEC %>% 
  rownames_to_column("gene_id")
LIF_CONEC <- expression_CONEC %>% 
  filter(gene_id == "LIF" | gene_id == "LIFR" | gene_id == "LIFR-AS1")
write.csv(LIF_CONEC, "LIF_CONEC-AEA.csv")

express_DKDEC <- GetAssayData(DKD_EC, layer = "data")
expression_DKDEC <- as.data.frame(express_DKDEC)
expression_DKDEC <- expression_DKDEC %>% 
  rownames_to_column("gene_id")
LIF_DKDEC <- expression_DKDEC %>% 
  filter(gene_id == "LIF" | gene_id == "LIFR" | gene_id == "LIFR-AS1")
write.csv(LIF_DKDEC, "LIF_DKDEC-AEA.csv")



####读取长数据基因（EC-AEA）###############
#读取长数据基因（EC-AEA）
human <- read.csv("./result/8.1.1 LIF-基因/长-LIF_EC-AEA.csv")

p1 <- ggplot(human,
             aes(x = group,
                 y = express,
                 color = group,
                 fill = group))+
  geom_bar(stat = "summary")+
  stat_summary(geom = "errorbar", width = 0.3)+
  theme_bw()+
  facet_wrap(~gene,scales = "free")
p1
filename <- "2.1 EC-AEA的Lif表达量总览.png"  
png(filename, width = 1000, height = 500)  
print(p1)  
dev.off()

p1 <- ggplot(human,
             aes(x = group,
                 y = express,
                 color = group,
                 fill = group))+
  geom_bar(stat = "summary")+
  stat_summary(geom = "errorbar", width = 0.3)+
  theme_bw()+
  facet_wrap(~gene,scales = "fixed")
p1
filename <- "2.2 EC-GC的Lif表达量总览（固定刻度）.png"  
png(filename, width = 1000, height = 500)  
print(p1)  
dev.off()
