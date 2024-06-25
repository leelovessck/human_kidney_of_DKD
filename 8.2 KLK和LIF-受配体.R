rm(list = ls())
gc()
setwd("D:/2023.10/调试/糖尿病肾病/2.鼠DKD分析")
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
                                       thresh.pc = 0.05, 
                                       thresh.fc = 0.05, 
                                       thresh.p = 1)



####提取DKD中改变的配受体对########################
#提取DKD中改变的配受体对
net <- netMappingDEG(cellchat, 
                     features.name = features.name)



####KLK#############
#KLK
net_KLK <- net %>%  
  filter(pathway_name == "KLK" & (source %in% "EC" | target %in% "EC")) 

unique(net_KLK$pathway_name)
net_KLK.use = net_KLK[, "interaction_name", drop = F]

net_KLK$interaction_name <- as.factor(net_KLK$interaction_name)
net_KLK$interaction_name <- fct_inorder(net_KLK$interaction_name)
net_KLK$interaction_name_2 <- as.factor(net_KLK$interaction_name_2)
net_KLK$interaction_name_2 <- fct_inorder(net_KLK$interaction_name_2)
net_KLK$pathway_name <- as.factor(net_KLK$pathway_name)
net_KLK$pathway_name <- fct_inorder(net_KLK$pathway_name)
net_KLK <- net_KLK %>% 
  mutate(p="")
net_KLK$signway <- paste(net_KLK$source,
                          "  ->  ", 
                          net_KLK$target) 


##配体
p1 <- ggplot(net_KLK,
             aes(x = signway,
                 y = interaction_name_2)) + 
  geom_point(aes(size = pval, 
                 color = ligand.logFC)) +
  geom_point(shape = 1,
             aes(size = pval), 
             color="black")+
  scale_size(rang = c(1.5,6)) +
  labs(x=NULL,y=NULL,size="P value",
       title = "KLK pathway") +
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
filename <- "1.1 EC的KLK气泡图A（配体）.png"
png(filename, width = 1000, height = 500)
print(p1)
dev.off()
  ##得到1.1 EC的KLK气泡图A（配体）


##受体
p1 <- ggplot(net_KLK,
             aes(x = signway,
                 y = interaction_name_2)) + 
  geom_point(aes(size = pval, 
                 color = receptor.logFC)) +
  geom_point(shape = 1,
             aes(size = pval), 
             color="black")+
  scale_size(rang = c(1.5,6)) +
  labs(x=NULL,y=NULL,size="P value",
       title = "KLK pathway") +
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
filename <- "1.2 EC的KLK气泡图A（受体）.png"
png(filename, width = 1000, height = 500)
print(p1)
dev.off()
  ##得到1.2 EC的KLK气泡图A（受体）



####LIFR#############
#LIFR
net_LIFR <- net %>%  
  filter(pathway_name == "LIFR" & (source %in% "EC" | target %in% "EC")) 

unique(net_LIFR$pathway_name)
net_LIFR.use = net_LIFR[, "interaction_name", drop = F]

net_LIFR$interaction_name <- as.factor(net_LIFR$interaction_name)
net_LIFR$interaction_name <- fct_inorder(net_LIFR$interaction_name)
net_LIFR$interaction_name_2 <- as.factor(net_LIFR$interaction_name_2)
net_LIFR$interaction_name_2 <- fct_inorder(net_LIFR$interaction_name_2)
net_LIFR$pathway_name <- as.factor(net_LIFR$pathway_name)
net_LIFR$pathway_name <- fct_inorder(net_LIFR$pathway_name)
net_LIFR <- net_LIFR %>% 
  mutate(p="")
net_LIFR$signway <- paste(net_LIFR$source,
                         "  ->  ", 
                         net_LIFR$target) 


##配体
p1 <- ggplot(net_LIFR,
             aes(x = signway,
                 y = interaction_name_2)) + 
  geom_point(aes(size = pval, 
                 color = ligand.logFC)) +
  geom_point(shape = 1,
             aes(size = pval), 
             color="black")+
  scale_size(rang = c(1.5,6)) +
  labs(x=NULL,y=NULL,size="P value",
       title = "LIFR pathway") +
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
filename <- "1.1 EC的LIFR气泡图A（配体）.png"
png(filename, width = 1000, height = 500)
print(p1)
dev.off()
##得到1.1 EC的LIFR气泡图A（配体）


##受体
p1 <- ggplot(net_LIFR,
             aes(x = signway,
                 y = interaction_name_2)) + 
  geom_point(aes(size = pval, 
                 color = receptor.logFC)) +
  geom_point(shape = 1,
             aes(size = pval), 
             color="black")+
  scale_size(rang = c(1.5,6)) +
  labs(x=NULL,y=NULL,size="P value",
       title = "LIFR pathway") +
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
filename <- "1.2 EC的LIFR气泡图A（受体）.png"
png(filename, width = 1000, height = 500)
print(p1)
dev.off()
##得到1.2 EC的LIFR气泡图A（受体）
