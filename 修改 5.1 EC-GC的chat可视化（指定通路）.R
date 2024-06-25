rm(list = ls())
gc()
setwd("D:/2023.10/调试/糖尿病肾病/1.人DKD分析/cellchat详细")
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
con.chat <- readRDS("./cellchat的rds/control（全部）.rds")
dkd.chat <- readRDS("./cellchat的rds/DKD（全部）.rds")
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
                                       thresh.pc = 0.1, 
                                       thresh.fc = 0.1, 
                                       thresh.p = 1)



####提取DKD中改变的配受体对########################
#提取DKD中改变的配受体对
net <- netMappingDEG(cellchat, 
                     features.name = features.name)



####TGFb通路###############################
#TGFb通路
net_TGFb <- net %>%  
  filter(pathway_name == "TGFb" & (source %in% "EC-GC" | target %in% "EC-GC")) 

exclude_values <- c("EC-PTC", "VSMC", "LOH", "PT", "DT", "CT", "EC-AEA")
net_TGFb <- net_TGFb %>%  
  filter(!(source %in% exclude_values | target %in% exclude_values))  

unique(net_TGFb$pathway_name)
net_TGFb.use = net_TGFb[, "interaction_name", drop = F]

net_TGFb$interaction_name <- as.factor(net_TGFb$interaction_name)
net_TGFb$interaction_name <- fct_inorder(net_TGFb$interaction_name)
net_TGFb$interaction_name_2 <- as.factor(net_TGFb$interaction_name_2)
net_TGFb$interaction_name_2 <- fct_inorder(net_TGFb$interaction_name_2)
net_TGFb$pathway_name <- as.factor(net_TGFb$pathway_name)
net_TGFb$pathway_name <- fct_inorder(net_TGFb$pathway_name)
net_TGFb <- net_TGFb %>% 
  mutate(p="")
net_TGFb$signway <- paste(net_TGFb$source,
                          "  ->  ", 
                          net_TGFb$target) 

p1 <- ggplot(net_TGFb,
             aes(x = signway,
                 y = interaction_name_2)) + 
  geom_point(aes(size = pval, 
                 color = ligand.logFC)) +
  geom_point(shape = 1,
             aes(size = pval), 
             color="black")+
  scale_size(rang = c(1.5,6)) +
  labs(x=NULL,y=NULL,size="P value",
       title = "TGFb pathway") +
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.text =element_text(size = 8, color = "black"),
        axis.text.y = element_text(face="italic"),
        axis.text.x=element_text(angle = 90,
                                 hjust = 1,
                                 vjust = 0.5)) +
  scale_color_gradient2(low = "navy",
                        high = "#FF0000",
                        midpoint = 0) +
  theme(plot.margin = margin(t = 5, 
                             r = 5, 
                             b = 5, 
                             l = 5, 
                             unit = "pt"))
p1
  ##得到1.1 EC-GC的TGFb气泡图A


gg1 <- netVisual_bubble(cellchat, 
                        pairLR.use = net_TGFb.use, 
                        sources.use = "EC-GC",
                        targets.use = c("B", "T", "FIB", "MC", 
                                        "POD", "EC-GC", "DC", 
                                        "MAC"),
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "Up-regulated Ligand of EC-GC in DKD")
gg1
  ##得到1.2 EC-GC（source）的TGFb气泡图B

gg2 <- netVisual_bubble(cellchat, 
                        pairLR.use = net_TGFb.use, 
                        targets.use = "EC-GC",
                        sources.use = c("B", "T", "FIB", "MC", 
                                        "POD", "EC-GC", "DC", 
                                        "MAC"),
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "Up-regulated Ligand of EC-GC in DKD")
gg2
  ##得到1.3 EC-GC（target）的TGFb气泡图B



####TRAIL通路###############################
#TRAIL通路
net_TRAIL <- net %>%  
  filter(pathway_name == "TRAIL" & (source %in% "EC-GC" | target %in% "EC-GC")) 

net_TRAIL <- net_TRAIL %>%  
  filter(!(source %in% exclude_values | target %in% exclude_values))  

unique(net_TRAIL$pathway_name)
net_TRAIL.use = net_TRAIL[, "interaction_name", drop = F]

net_TRAIL$interaction_name <- as.factor(net_TRAIL$interaction_name)
net_TRAIL$interaction_name <- fct_inorder(net_TRAIL$interaction_name)
net_TRAIL$interaction_name_2 <- as.factor(net_TRAIL$interaction_name_2)
net_TRAIL$interaction_name_2 <- fct_inorder(net_TRAIL$interaction_name_2)
net_TRAIL$pathway_name <- as.factor(net_TRAIL$pathway_name)
net_TRAIL$pathway_name <- fct_inorder(net_TRAIL$pathway_name)
net_TRAIL <- net_TRAIL %>% 
  mutate(p="")
net_TRAIL$signway <- paste(net_TRAIL$source,
                          "  ->  ", 
                          net_TRAIL$target) 

p1 <- ggplot(net_TRAIL,
             aes(x = signway,
                 y = interaction_name_2)) + 
  geom_point(aes(size = pval, 
                 color = ligand.logFC)) +
  geom_point(shape = 1,
             aes(size = pval), 
             color="black")+
  scale_size(rang = c(1.5,6)) +
  labs(x=NULL,y=NULL,size="P value",
       title = "TRAIL pathway") +
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.text =element_text(size = 8, color = "black"),
        axis.text.y = element_text(face="italic"),
        axis.text.x=element_text(angle = 90,
                                 hjust = 1,
                                 vjust = 0.5)) +
  scale_color_gradient2(low = "navy",
                        high = "#FF0000",
                        midpoint = 0) +
  theme(plot.margin = margin(t = 5, 
                             r = 5, 
                             b = 5, 
                             l = 5, 
                             unit = "pt"))
p1
  ##得到2.1 EC-GC的TRAIL气泡图A


gg1 <- netVisual_bubble(cellchat, 
                        pairLR.use = net_TRAIL.use, 
                        sources.use = "EC-GC",
                        targets.use = c("B", "T", "FIB", "MC", 
                                        "POD", "EC-GC", "DC", 
                                        "MAC"),
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "Up-regulated Ligand of EC-GC in DKD")
gg1
  ##得到2.2 EC-GC（source）的TRAIL气泡图B

gg2 <- netVisual_bubble(cellchat, 
                        pairLR.use = net_TRAIL.use, 
                        targets.use = "EC-GC",
                        sources.use = c("B", "T", "FIB", "MC", 
                                        "POD", "EC-GC", "DC", 
                                        "MAC"),
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "Up-regulated Ligand of EC-GC in DKD")
gg2
  ##得到2.3 EC-GC（target）的TRAIL气泡图B



####VISFATIN通路###############################
#VISFATIN通路
net_VISFATIN <- net %>%  
  filter(pathway_name == "VISFATIN" & (source %in% "EC-GC" | target %in% "EC-GC")) 

net_VISFATIN <- net_VISFATIN %>%  
  filter(!(source %in% exclude_values | target %in% exclude_values))  

unique(net_VISFATIN$pathway_name)
net_VISFATIN.use = net_VISFATIN[, "interaction_name", drop = F]

net_VISFATIN$interaction_name <- as.factor(net_VISFATIN$interaction_name)
net_VISFATIN$interaction_name <- fct_inorder(net_VISFATIN$interaction_name)
net_VISFATIN$interaction_name_2 <- as.factor(net_VISFATIN$interaction_name_2)
net_VISFATIN$interaction_name_2 <- fct_inorder(net_VISFATIN$interaction_name_2)
net_VISFATIN$pathway_name <- as.factor(net_VISFATIN$pathway_name)
net_VISFATIN$pathway_name <- fct_inorder(net_VISFATIN$pathway_name)
net_VISFATIN <- net_VISFATIN %>% 
  mutate(p="")
net_VISFATIN$signway <- paste(net_VISFATIN$source,
                           "  ->  ", 
                           net_VISFATIN$target) 

p1 <- ggplot(net_VISFATIN,
             aes(x = signway,
                 y = interaction_name_2)) + 
  geom_point(aes(size = pval, 
                 color = ligand.logFC)) +
  geom_point(shape = 1,
             aes(size = pval), 
             color="black")+
  scale_size(rang = c(1.5,6)) +
  labs(x=NULL,y=NULL,size="P value",
       title = "VISFATIN pathway") +
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.text =element_text(size = 8, color = "black"),
        axis.text.y = element_text(face="italic"),
        axis.text.x=element_text(angle = 90,
                                 hjust = 1,
                                 vjust = 0.5)) +
  scale_color_gradient2(low = "navy",
                        high = "#FF0000",
                        midpoint = 0) +
  theme(plot.margin = margin(t = 5, 
                             r = 5, 
                             b = 5, 
                             l = 5, 
                             unit = "pt"))
p1
  ##得到3.1 EC-GC的VISFATIN气泡图A


gg1 <- netVisual_bubble(cellchat, 
                        pairLR.use = net_VISFATIN.use, 
                        sources.use = "EC-GC",
                        targets.use = c("B", "T", "FIB", "MC", 
                                        "POD", "EC-GC", "DC", 
                                        "MAC"),
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "Up-regulated Ligand of EC-GC in DKD")
gg1
  ##得到3.2 EC-GC（source）的VISFATIN气泡图B--这里无结果

gg2 <- netVisual_bubble(cellchat, 
                        pairLR.use = net_VISFATIN.use, 
                        targets.use = "EC-GC",
                        sources.use = c("B", "T", "FIB", "MC", 
                                        "POD", "EC-GC", "DC", 
                                        "MAC"),
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "Up-regulated Ligand of EC-GC in DKD")
gg2
  ##得到3.3 EC-GC（target）的VISFATIN气泡图B



####PTN通路###############################
#PTN通路
net_PTN <- net %>%  
  filter(pathway_name == "PTN" & (source %in% "EC-GC" | target %in% "EC-GC")) 

net_PTN <- net_PTN %>%  
  filter(!(source %in% exclude_values | target %in% exclude_values))  

unique(net_PTN$pathway_name)
net_PTN.use = net_PTN[, "interaction_name", drop = F]

net_PTN$interaction_name <- as.factor(net_PTN$interaction_name)
net_PTN$interaction_name <- fct_inorder(net_PTN$interaction_name)
net_PTN$interaction_name_2 <- as.factor(net_PTN$interaction_name_2)
net_PTN$interaction_name_2 <- fct_inorder(net_PTN$interaction_name_2)
net_PTN$pathway_name <- as.factor(net_PTN$pathway_name)
net_PTN$pathway_name <- fct_inorder(net_PTN$pathway_name)
net_PTN <- net_PTN %>% 
  mutate(p="")
net_PTN$signway <- paste(net_PTN$source,
                              "  ->  ", 
                              net_PTN$target) 

p1 <- ggplot(net_PTN,
             aes(x = signway,
                 y = interaction_name_2)) + 
  geom_point(aes(size = pval, 
                 color = ligand.logFC)) +
  geom_point(shape = 1,
             aes(size = pval), 
             color="black")+
  scale_size(rang = c(1.5,6)) +
  labs(x=NULL,y=NULL,size="P value",
       title = "PTN pathway") +
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.text =element_text(size = 8, color = "black"),
        axis.text.y = element_text(face="italic"),
        axis.text.x=element_text(angle = 90,
                                 hjust = 1,
                                 vjust = 0.5)) +
  scale_color_gradient2(low = "navy",
                        high = "#FF0000",
                        midpoint = 0) +
  theme(plot.margin = margin(t = 5, 
                             r = 5, 
                             b = 5, 
                             l = 5, 
                             unit = "pt"))
p1
  ##得到4.1 EC-GC的PTN气泡图A


gg1 <- netVisual_bubble(cellchat, 
                        pairLR.use = net_PTN.use, 
                        sources.use = "EC-GC",
                        targets.use = c("B", "T", "FIB", "MC", 
                                        "POD", "EC-GC", "DC", 
                                        "MAC"),
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "Up-regulated Ligand of EC-GC in DKD")
gg1
  ##得到4.2 EC-GC（source）的PTN气泡图B--这里无结果

gg2 <- netVisual_bubble(cellchat, 
                        pairLR.use = net_PTN.use, 
                        targets.use = "EC-GC",
                        sources.use = c("B", "T", "FIB", "MC", 
                                        "POD", "EC-GC", "DC", 
                                        "MAC"),
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "Up-regulated Ligand of EC-GC in DKD")
gg2
  ##得到4.3 EC-GC（target）的PTN气泡图B



####ANGPT通路###############################
#ANGPT通路
net_ANGPT <- net %>%  
  filter(pathway_name == "ANGPT" & (source %in% "EC-GC" | target %in% "EC-GC")) 

net_ANGPT <- net_ANGPT %>%  
  filter(!(source %in% exclude_values | target %in% exclude_values))  

unique(net_ANGPT$pathway_name)
net_ANGPT.use = net_ANGPT[, "interaction_name", drop = F]

net_ANGPT$interaction_name <- as.factor(net_ANGPT$interaction_name)
net_ANGPT$interaction_name <- fct_inorder(net_ANGPT$interaction_name)
net_ANGPT$interaction_name_2 <- as.factor(net_ANGPT$interaction_name_2)
net_ANGPT$interaction_name_2 <- fct_inorder(net_ANGPT$interaction_name_2)
net_ANGPT$pathway_name <- as.factor(net_ANGPT$pathway_name)
net_ANGPT$pathway_name <- fct_inorder(net_ANGPT$pathway_name)
net_ANGPT <- net_ANGPT %>% 
  mutate(p="")
net_ANGPT$signway <- paste(net_ANGPT$source,
                         "  ->  ", 
                         net_ANGPT$target) 

p1 <- ggplot(net_ANGPT,
             aes(x = signway,
                 y = interaction_name_2)) + 
  geom_point(aes(size = pval, 
                 color = ligand.logFC)) +
  geom_point(shape = 1,
             aes(size = pval), 
             color="black")+
  scale_size(rang = c(1.5,6)) +
  labs(x=NULL,y=NULL,size="P value",
       title = "ANGPT pathway") +
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.text =element_text(size = 8, color = "black"),
        axis.text.y = element_text(face="italic"),
        axis.text.x=element_text(angle = 90,
                                 hjust = 1,
                                 vjust = 0.5)) +
  scale_color_gradient2(low = "navy",
                        high = "#FF0000",
                        midpoint = 0) +
  theme(plot.margin = margin(t = 5, 
                             r = 5, 
                             b = 5, 
                             l = 5, 
                             unit = "pt"))
p1
  ##得到5.1 EC-GC的ANGPT气泡图A


gg1 <- netVisual_bubble(cellchat, 
                        pairLR.use = net_ANGPT.use, 
                        sources.use = "EC-GC",
                        targets.use = c("B", "T", "FIB", "MC", 
                                        "POD", "EC-GC", "DC", 
                                        "MAC"),
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "Up-regulated Ligand of EC-GC in DKD")
gg1
  ##得到5.2 EC-GC（source）的ANGPT气泡图B

gg2 <- netVisual_bubble(cellchat, 
                        pairLR.use = net_ANGPT.use, 
                        targets.use = "EC-GC",
                        sources.use = c("B", "T", "FIB", "MC", 
                                        "POD", "EC-GC", "DC", 
                                        "MAC"),
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "Up-regulated Ligand of EC-GC in DKD")
gg2
  ##得到5.3 EC-GC（target）的ANGPT气泡图B
