rm(list = ls())
gc()
setwd("D:/2023.10/调试/糖尿病肾病/1.人DKD分析")
getwd()

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
if(!require(ggalluvial))install.packages('ggalluvial')
if(!require(igraph))install.packages('igraph')
if(!require(forcats))install.packages("forcats")
if(!require(aplot))install.packages("aplot")
}



####数据准备############################
#数据准备
con.chat <- readRDS("./data source/cellchat/control（全部）.rds")
dkd.chat <- readRDS("./data source/cellchat/DKD（全部）.rds")
object.list <- list(CON = con.chat, DKD = dkd.chat)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
levels(object.list[[1]]@idents) 
levels(object.list[[2]]@idents) 
group.cellType <- levels(object.list[[1]]@idents)
group.cellType <- factor(group.cellType, levels = levels(object.list[[1]]@idents))
object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})#合并细胞类型
con.chat <- netAnalysis_computeCentrality(con.chat) 
dkd.chat <- netAnalysis_computeCentrality(dkd.chat) 
object.list <- list(CON = con.chat, DKD = dkd.chat)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))



####差异计算#########################
#差异计算
pos.dataset = "DKD"
features.name = pos.dataset
cellchat <- identifyOverExpressedGenes(cellchat, 
                                       group.dataset = "datasets", 
                                       pos.dataset = pos.dataset, 
                                       features.name = features.name, 
                                       only.pos = FALSE, 
                                       thresh.pc = 0, 
                                       thresh.fc = 0, 
                                       thresh.p = 1)



####提取DKD中改变的配受体对（logFC = 0.2）########################
#提取DKD中改变的配受体对（logFC = 0.2）
net <- netMappingDEG(cellchat, 
                     features.name = features.name)




gg1 <- netVisual_bubble(cellchat,
                        signaling = "CSPG4",
                        comparison = c(1, 2), 
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "CSPG4")
gg1
  ##得到1.CSPG4的气泡图A


net.cspg4 <- net %>%  
  filter(pathway_name == "CSPG4") %>% 
  mutate(p="")
net.cspg4$signway <- paste(net.cspg4$source,
                     "  ->  ", 
                     net.cspg4$target) 
p1 <- ggplot(net.cspg4,
             aes(x = signway,
                 y = interaction_name_2)) + 
  geom_point(aes(size = pval, 
                 color = ligand.logFC)) +
  geom_point(shape = 1,
             aes(size = pval), 
             color="black")+
  scale_size(rang = c(1.5,6)) +
  labs(x=NULL,y=NULL,size="P value",
       title = "CSPG4") +
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.text =element_text(size = 8, color = "black"),
        axis.text.y = element_text(face="italic"),
        axis.text.x=element_text(angle = 90,
                                 hjust = 1,
                                 vjust = 0.5)) +
  scale_color_gradient2(low='navy',
                        high='RED',
                        midpoint = 0)
p1
  ##得到2.CSPG4的气泡图B



gg1 <- netVisual_bubble(cellchat,
                        signaling = "ADGRA",
                        comparison = c(1, 2), 
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "ADGRA")
gg1
  ##得到3.ADGRA的气泡图A

net.adgra <- net %>%  
  filter(pathway_name == "ADGRA") %>% 
  mutate(p="")
net.adgra$signway <- paste(net.adgra$source,
                           "  ->  ", 
                           net.adgra$target) 
net.adgra <- net %>%  
  filter(pathway_name == "ADGRA")
p1 <- ggplot(net.adgra,
             aes(x = signway,
                 y = interaction_name_2)) + 
  geom_point(aes(size = pval, 
                 color = ligand.logFC)) +
  geom_point(shape = 1,
             aes(size = pval), 
             color="black")+
  scale_size(rang = c(1.5,6)) +
  labs(x=NULL,y=NULL,size="P value",
       title = "CSPG4") +
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.text =element_text(size = 8, color = "black"),
        axis.text.y = element_text(face="italic"),
        axis.text.x=element_text(angle = 90,
                                 hjust = 1,
                                 vjust = 0.5)) +
  scale_color_gradient2(low='navy',
                        high='RED',
                        midpoint = 0)
p1
  ##得到4.ADGRA的气泡图B



gg1 <- netVisual_bubble(cellchat,
                        signaling = "ADGRL",
                        comparison = c(1, 2), 
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "ADGRL")
gg1
##得到3.ADGRL的气泡图A


net.adgrl <- net %>%  
  filter(pathway_name == "ADGRL")
p1 <- ggplot(net.adgrl,
             aes(x = signway,
                 y = interaction_name_2)) + 
  geom_point(aes(size = pval, 
                 color = ligand.logFC)) +
  geom_point(shape = 1,
             aes(size = pval), 
             color="black")+
  scale_size(rang = c(1.5,6)) +
  labs(x=NULL,y=NULL,size="P value",
       title = "CSPG4") +
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.text =element_text(size = 8, color = "black"),
        axis.text.y = element_text(face="italic"),
        axis.text.x=element_text(angle = 90,
                                 hjust = 1,
                                 vjust = 0.5)) +
  scale_color_gradient2(low='navy',
                        high='RED',
                        midpoint = 0)
p1
##得到5.ADGRL的气泡图B
