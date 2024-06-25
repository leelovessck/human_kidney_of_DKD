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
                                       thresh.pc = 0.1, 
                                       thresh.fc = 0.1, 
                                       thresh.p = 1)



####提取DKD中改变的配受体对（logFC = 0.2）########################
#提取DKD中改变的配受体对（logFC = 0.2）
net <- netMappingDEG(cellchat, 
                     features.name = features.name)
ligand.up <- subsetCommunication(cellchat, 
                              net = net, 
                              datasets = "DKD",
                              ligand.logFC = 0.2, 
                              receptor.logFC = NULL)
receptor.up <- subsetCommunication(cellchat, 
                                net = net, 
                                datasets = "DKD",
                                ligand.logFC = NULL, 
                                receptor.logFC = 0.2)
ligand.down <- subsetCommunication(cellchat, 
                                 net = net, 
                                 datasets = "DKD",
                                 ligand.logFC = -0.2, 
                                 receptor.logFC = NULL)
receptor.down <- subsetCommunication(cellchat, 
                                   net = net, 
                                   datasets = "DKD",
                                   ligand.logFC = NULL, 
                                   receptor.logFC = -0.2)



####提取EC-AEA的上调配体（logFC = 0.2）#######################
#提取EC-AEA的上调配体（logFC = 0.2）
AEA.s.lup <- ligand.up %>%  
  filter(source == "EC-AEA")
AEA.t.lup <- ligand.up %>%  
  filter(target == "EC-AEA")



####EC-AEA的上调配体气泡图（logFC = 0.2）#####################
#EC-AEA的上调配体气泡图（logFC = 0.2）
##EC-AEA作为source
unique(AEA.s.lup$pathway_name)
AEA.s.lup.use = AEA.s.lup[, "interaction_name", drop = F]
AEA.s.lup$interaction_name <- as.factor(AEA.s.lup$interaction_name)
AEA.s.lup$interaction_name <- fct_inorder(AEA.s.lup$interaction_name)
AEA.s.lup$interaction_name_2 <- as.factor(AEA.s.lup$interaction_name_2)
AEA.s.lup$interaction_name_2 <- fct_inorder(AEA.s.lup$interaction_name_2)
AEA.s.lup$pathway_name <- as.factor(AEA.s.lup$pathway_name)
AEA.s.lup$pathway_name <- fct_inorder(AEA.s.lup$pathway_name)
AEA.s.lup <- AEA.s.lup %>% 
  mutate(p="")
AEA.s.lup$signway <- paste(AEA.s.lup$source,
                          "  ->  ", 
                          AEA.s.lup$target) 

table(unique(AEA.s.lup$interaction_name_2))
summary_data <- AEA.s.lup %>%  
  group_by(pathway_name) %>%
  summarise(interaction_count = n_distinct(interaction_name_2))
print(summary_data)
write.csv(summary_data, "data.csv")

p1 <- ggplot(AEA.s.lup,
             aes(x = signway,
                 y = interaction_name_2)) + 
  geom_point(aes(size = pval, 
                 color = ligand.logFC)) +
  geom_point(shape = 1,
             aes(size = pval), 
             color="black")+
  scale_size(rang = c(1.5,6)) +
  labs(x=NULL,y=NULL,size="P value",
       title = "Up-regulated signaling of EC-AEA in DKD") +
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.text =element_text(size = 8, color = "black"),
        axis.text.y = element_text(face="italic"),
        axis.text.x=element_text(angle = 90,
                                 hjust = 1,
                                 vjust = 0.5)) +
  scale_color_gradient2(low='#FFB6C1',
                        high='#FF0000') +
  geom_hline(yintercept=c(3.5, 4.5, 5.5, 11.5, 12.5))+
  theme(plot.margin = margin(t = 5, 
                             r = 5, 
                             b = 5, 
                             l = 5, 
                             unit = "pt"))
p1

labels_and_y <- list("MIF" = 2, "SEMA3" = 4, 
                     "PSAP" = 5, "EPHA" = 8.5, 
                     "MHC-II" = 12, "PECAM1" = 13) 
gene_group <- ggplot(AEA.s.lup,
                     aes(x = p,
                         y = interaction_name_2,
                         fill = pathway_name))+
  geom_tile() + 
  scale_y_discrete(position="left") +
  theme_minimal()+xlab(NULL) + ylab(NULL) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x =element_text(
          angle =90,hjust =0.5,vjust = 0.5),
        legend.position = 'none')+
  lapply(names(labels_and_y), function(label) {  
    y_val <- labels_and_y[[label]]  
    annotate("text", label=label, x=1.2, y=y_val)})+  
  theme(plot.margin = margin(t = 0, 
                             r = 0, 
                             b = 0, 
                             l = 0, 
                             unit = "pt"))
gene_group

p1%>%insert_left(gene_group, width = 0.3)
  ##得到1.1 EC-AEA的source上调配体气泡图A（logFC = 0.2）


gg1 <- netVisual_bubble(cellchat, 
                        pairLR.use = AEA.s.lup.use, 
                        sources.use = "EC-AEA", 
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "Up-regulated signaling of EC-AEA in DKD")
gg1 <- gg1 + 
  geom_hline(yintercept=c(3.5, 4.5, 5.5, 11.5, 12.5))
gg1 %>% insert_left(gene_group, width = 0.3)
  ##得到1.2 EC-AEA的source上调配体气泡图B（logFC = 0.2）


##EC-AEA作为target
AEA.t.lup <- AEA.t.lup %>%  
  arrange(substr(pathway_name, 1, 1)) 

unique(AEA.t.lup$pathway_name)
AEA.t.lup.use = AEA.t.lup[, "interaction_name", drop = F]
AEA.t.lup$interaction_name <- as.factor(AEA.t.lup$interaction_name)
AEA.t.lup$interaction_name <- fct_inorder(AEA.t.lup$interaction_name)
AEA.t.lup$interaction_name_2 <- as.factor(AEA.t.lup$interaction_name_2)
AEA.t.lup$interaction_name_2 <- fct_inorder(AEA.t.lup$interaction_name_2)
AEA.t.lup$pathway_name <- as.factor(AEA.t.lup$pathway_name)
AEA.t.lup$pathway_name <- fct_inorder(AEA.t.lup$pathway_name)
AEA.t.lup <- AEA.t.lup %>% 
  mutate(p="")

AEA.t.lup$signway <- paste(AEA.t.lup$source,
                          "  ->  ", 
                          AEA.t.lup$target) 

table(unique(AEA.t.lup$interaction_name_2))
summary_data <- AEA.t.lup %>%  
  group_by(pathway_name) %>%
  summarise(interaction_count = n_distinct(interaction_name_2))
print(summary_data)
write.csv(summary_data, "data.csv")

p1 <- ggplot(AEA.t.lup,
             aes(x = signway,
                 y = interaction_name_2)) + 
  geom_point(aes(size = pval, 
                 color = ligand.logFC)) +
  geom_point(shape = 1,
             aes(size = pval), 
             color="black")+
  scale_size(rang = c(1.5,6)) +
  labs(x=NULL,y=NULL,size="P value",
       title = "Up-regulated signaling of EC-AEA in DKD") +
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.text =element_text(size = 8, color = "black"),
        axis.text.y = element_text(face="italic"),
        axis.text.x=element_text(angle = 90,
                                 hjust = 1,
                                 vjust = 0.5)) +
  scale_color_gradient2(low='#FFB6C1',
                        high='#FF0000',
                        midpoint = 0.8) +
  geom_hline(yintercept=c(1.5, 2.5, 3.5, 6.5, 24.5, 25.5, 26.5, 27.5, 
                          29.5, 31.5, 44.5, 47.5, 53.5, 54.5, 55.5, 
                          61.5, 64.5, 68.5, 71.5, 77.5, 78.5, 79.5, 
                          80.5, 83.5, 85.5, 86.5, 87.5, 93.5, 96.5, 
                          97.5, 98.5, 100.5, 102.5, 104.5))+
  theme(plot.margin = margin(t = 5, 
                             r = 5, 
                             b = 5, 
                             l = 5, 
                             unit = "pt"))
p1

labels_and_y <- list("APP" = 1, "CCL" = 2, 
                     "CXCL" = 3, "COMPLEMENT" = 5, 
                     "COLLAGEN" = 15.5, "CD45" = 25, 
                     "CD46" = 26, "CD48" = 27, 
                     "CDH" = 28.5, "CDH1" = 30.5, 
                     "EPHA" = 38, "EPHB" = 46, 
                     "FN1" = 50.5, "GAS" = 54, 
                     "LT" = 55, "LAMININ" = 58.5, 
                     "MIF" = 63, "MK" = 66.5, 
                     "MHC-I" = 70, "MHC-II" = 74.5, 
                     "NRG" = 78, "Netrin" = 79, 
                     "OCLN" = 80, "PDGF" = 82, 
                     "PTN" = 84.5, "PSAP" = 86, 
                     "PECAM1" = 87, "SPP1" = 90.5, 
                     "SEMA3" = 95, "TRAIL" = 97, 
                     "TENASCIN" = 98, "THBS" = 99.5, 
                     "THY1" = 101.5, "VCAM" = 103.5, 
                     "VISTA" = 105)  
gene_group <- ggplot(AEA.t.lup,
                     aes(x = p,
                         y = interaction_name_2,
                         fill = pathway_name))+
  geom_tile() + 
  scale_y_discrete(position="left") +
  theme_minimal()+xlab(NULL) + ylab(NULL) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x =element_text(
          angle =90,hjust =0.5,vjust = 0.5),
        legend.position = 'none')+
  lapply(names(labels_and_y), function(label) {  
    y_val <- labels_and_y[[label]]  
    annotate("text", label=label, x=1.2, y=y_val)})+
  theme(plot.margin = margin(t = 0, 
                             r = 0, 
                             b = 0, 
                             l = 0, 
                             unit = "pt"))
gene_group
p1 %>% insert_left(gene_group, width = 0.3)
  ##得到1.3 EC-AEA的target上调配体气泡图A（logFC = 0.2）

gg1 <- netVisual_bubble(cellchat, 
                        pairLR.use = AEA.t.lup.use, 
                        targets.use = "EC-AEA", 
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "Up-regulated signaling of EC-AEA in DKD")
gg1 <- gg1 + 
  geom_hline(yintercept=c(1.5, 2.5, 3.5, 6.5, 24.5, 25.5, 26.5, 27.5, 
                          29.5, 31.5, 44.5, 47.5, 53.5, 54.5, 55.5, 
                          61.5, 64.5, 68.5, 71.5, 77.5, 78.5, 79.5, 
                          80.5, 83.5, 85.5, 86.5, 87.5, 93.5, 96.5, 
                          97.5, 98.5, 100.5, 102.5, 104.5))
gg1 %>% insert_left(gene_group, width = 0.3)
  ##得到1.4 EC-AEA的target上调配体气泡图B（logFC = 0.2）



####提取EC-AEA的上调受体（logFC = 0.2）#######################
#提取EC-AEA的上调受体（logFC = 0.2）
AEA.s.rup <- receptor.up %>%  
  filter(source == "EC-AEA")
AEA.t.rup <- receptor.up %>%  
  filter(target == "EC-AEA")



####EC-AEA的上调受体气泡图（logFC = 0.2）#####################
#EC-AEA的上调受体气泡图（logFC = 0.2）
##EC-AEA作为source
AEA.s.rup <- AEA.s.rup %>%  
  arrange(substr(pathway_name, 1, 1)) 

unique(AEA.s.rup$pathway_name)
AEA.s.rup.use = AEA.s.rup[, "interaction_name", drop = F]
AEA.s.rup$interaction_name <- as.factor(AEA.s.rup$interaction_name)
AEA.s.rup$interaction_name <- fct_inorder(AEA.s.rup$interaction_name)
AEA.s.rup$interaction_name_2 <- as.factor(AEA.s.rup$interaction_name_2)
AEA.s.rup$interaction_name_2 <- fct_inorder(AEA.s.rup$interaction_name_2)
AEA.s.rup$pathway_name <- as.factor(AEA.s.rup$pathway_name)
AEA.s.rup$pathway_name <- fct_inorder(AEA.s.rup$pathway_name)
AEA.s.rup <- AEA.s.rup %>% 
  mutate(p="")
AEA.s.rup$signway <- paste(AEA.s.rup$source,
                          "  ->  ", 
                          AEA.s.rup$target) 

table(unique(AEA.s.rup$interaction_name_2))
summary_data <- AEA.s.rup %>%  
  group_by(pathway_name) %>%
  summarise(interaction_count = n_distinct(interaction_name_2))
print(summary_data)
write.csv(summary_data, "data.csv")

p1 <- ggplot(AEA.s.rup,
             aes(x = signway,
                 y = interaction_name_2)) + 
  geom_point(aes(size = pval, 
                 color = receptor.logFC)) +
  geom_point(shape = 1,
             aes(size = pval), 
             color="black")+
  scale_size(rang = c(1.5,6)) +
  labs(x=NULL,y=NULL,size="P value",
       title = "Up-regulated signaling of EC-AEA in DKD") +
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.text =element_text(size = 8, color = "black"),
        axis.text.y = element_text(face="italic"),
        axis.text.x=element_text(angle = 90,
                                 hjust = 1,
                                 vjust = 0.5)) +
  scale_color_gradient2(low='#FFB6C1',
                        high='#FF0000') +
  geom_hline(yintercept=c(1.5, 2.5, 6.5, 9.5, 10.5, 11.5, 12.5, 28.5,
                          29.5, 31.5, 32.5, 33.5, 43.5, 44.5, 45.5, 
                          47.5, 50.5, 51.5, 52.5, 53.5, 55.5, 57.5,
                          78.5, 79.5, 80.5, 85.5, 87.5, 93.5, 94.5,
                          95.5, 96.5, 97.5, 98.5, 107.5, 110.5, 114.5,
                          121.5, 123.5))+
  theme(plot.margin = margin(t = 5, 
                             r = 5, 
                             b = 5, 
                             l = 5, 
                             unit = "pt"))
p1

labels_and_y <- list("ANGPT" = 1, "APP" = 2, 
                     "ADGRG" = 4.5, "BMP" = 8, 
                     "BAFF" = 10, "CD40" = 11, 
                     "CALCR" = 12, "COLLAGEN" = 20.5, 
                     "CD34" = 29, "CDH" = 30.5, 
                     "CNTN" = 32, "EDN" = 33, 
                     "FN1" = 38.5, "GRN" = 44, 
                     "HGF" = 45, "IGF" = 46.5, 
                     "IL2" = 49, "IL6" = 51, 
                     "IGFBP" = 52, "ICAM" = 53, 
                     "JAM" = 54.5, "LIFR" = 56.5, 
                     "LAMININ" = 68, "L1CAM" = 79, 
                     "MIF" = 80, "MK" = 83, 
                     "NT" = 86.5, "NOTCH" = 90.5, 
                     "OSM" = 94, "OCLN" = 95, 
                     "PECAM1" = 96, "RELN" = 97, 
                     "SIRP" = 98, "TGFb" = 103, 
                     "THBS" = 109, "TENASCIN" = 112.5, 
                     "VEGF" = 118, "VCAM" = 122.5, 
                     "WNT" = 124.5) 
gene_group <- ggplot(AEA.s.rup,
                     aes(x = p,
                         y = interaction_name_2,
                         fill = pathway_name))+
  geom_tile() + 
  scale_y_discrete(position="left") +
  theme_minimal()+xlab(NULL) + ylab(NULL) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x =element_text(
          angle =90,hjust =0.5,vjust = 0.5),
        legend.position = 'none')+
  lapply(names(labels_and_y), function(label) {  
    y_val <- labels_and_y[[label]]  
    annotate("text", label=label, x=1.2, y=y_val)})+  
  theme(plot.margin = margin(t = 0, 
                             r = 0, 
                             b = 0, 
                             l = 0, 
                             unit = "pt"))
gene_group

p1 %>% insert_left(gene_group, width = 0.3)
  ##得到2.1 EC-AEA的source上调受体气泡图A（logFC = 0.2）


gg1 <- netVisual_bubble(cellchat, 
                        pairLR.use = AEA.s.rup.use, 
                        sources.use = "EC-AEA", 
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "Up-regulated signaling of EC-AEA in DKD")
gg1 <- gg1 + 
  geom_hline(yintercept=c(1.5, 2.5, 6.5, 9.5, 10.5, 11.5, 12.5, 28.5,
                          29.5, 31.5, 32.5, 33.5, 43.5, 44.5, 45.5, 
                          47.5, 50.5, 51.5, 52.5, 53.5, 55.5, 57.5,
                          78.5, 79.5, 80.5, 85.5, 87.5, 93.5, 94.5,
                          95.5, 96.5, 97.5, 98.5, 107.5, 110.5, 114.5,
                          121.5, 123.5))
gg1 %>% insert_left(gene_group, width = 0.3)
  ##得到2.2 EC-AEA的source上调受体气泡图B（logFC = 0.2）


##EC-AEA作为target
AEA.t.rup <- AEA.t.rup %>%  
  arrange(substr(pathway_name, 1, 1)) 

unique(AEA.t.rup$pathway_name)
AEA.t.rup.use = AEA.t.rup[, "interaction_name", drop = F]
AEA.t.rup$interaction_name <- as.factor(AEA.t.rup$interaction_name)
AEA.t.rup$interaction_name <- fct_inorder(AEA.t.rup$interaction_name)
AEA.t.rup$interaction_name_2 <- as.factor(AEA.t.rup$interaction_name_2)
AEA.t.rup$interaction_name_2 <- fct_inorder(AEA.t.rup$interaction_name_2)
AEA.t.rup$pathway_name <- as.factor(AEA.t.rup$pathway_name)
AEA.t.rup$pathway_name <- fct_inorder(AEA.t.rup$pathway_name)
AEA.t.rup <- AEA.t.rup %>% 
  mutate(p="")

AEA.t.rup$signway <- paste(AEA.t.rup$source,
                          "  ->  ", 
                          AEA.t.rup$target) 

table(unique(AEA.t.rup$interaction_name_2))
summary_data <- AEA.t.rup %>%  
  group_by(pathway_name) %>%
  summarise(interaction_count = n_distinct(interaction_name_2))
print(summary_data)
write.csv(summary_data, "data.csv")

p1 <- ggplot(AEA.t.rup,
             aes(x = signway,
                 y = interaction_name_2)) + 
  geom_point(aes(size = pval, 
                 color = receptor.logFC)) +
  geom_point(shape = 1,
             aes(size = pval), 
             color="black")+
  scale_size(rang = c(1.5,6)) +
  labs(x=NULL,y=NULL,size="P value",
       title = "Up-regulated signaling of EC-AEA in DKD") +
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.text =element_text(size = 8, color = "black"),
        axis.text.y = element_text(face="italic"),
        axis.text.x=element_text(angle = 90,
                                 hjust = 1,
                                 vjust = 0.5)) +
  scale_color_gradient2(low='#FFB6C1',
                        high='#FF0000') +
  geom_hline(yintercept=c(1.5, 2.5, 6.5))+
  theme(plot.margin = margin(t = 5, 
                             r = 5, 
                             b = 5, 
                             l = 5, 
                             unit = "pt"))
p1

labels_and_y <- list("PECAM1" = 1, "THBS" = 2, 
                     "VEGF" = 4.5, "VCAM" = 7.5)  
gene_group <- ggplot(AEA.t.rup,
                     aes(x = p,
                         y = interaction_name_2,
                         fill = pathway_name))+
  geom_tile() + 
  scale_y_discrete(position="left") +
  theme_minimal()+xlab(NULL) + ylab(NULL) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x =element_text(
          angle =90,hjust =0.5,vjust = 0.5),
        legend.position = 'none')+
  lapply(names(labels_and_y), function(label) {  
    y_val <- labels_and_y[[label]]  
    annotate("text", label=label, x=1.2, y=y_val)})+
  theme(plot.margin = margin(t = 0, 
                             r = 0, 
                             b = 0, 
                             l = 0, 
                             unit = "pt"))
gene_group
p1 %>% insert_left(gene_group, width = 0.3)
  ##得到2.3 EC-AEA的target上调受体气泡图A（logFC = 0.2）

gg1 <- netVisual_bubble(cellchat, 
                        pairLR.use = AEA.t.rup.use, 
                        targets.use = "EC-AEA", 
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "Up-regulated signaling of EC-AEA in DKD")
gg1 <- gg1 + 
  geom_hline(yintercept=c(1.5, 2.5, 6.5))
gg1 %>% insert_left(gene_group, width = 0.3)
  ##得到2.4 EC-AEA的target上调受体气泡图B（logFC = 0.2）



####提取EC-AEA的下调配体（logFC = 0.2）#######################
#提取EC-AEA的下调配体（logFC = 0.2）
AEA.s.ldown <- ligand.down %>%  
  filter(source == "EC-AEA")
AEA.t.ldown <- ligand.down %>%  
  filter(target == "EC-AEA")



####EC-AEA的下调配体气泡图（logFC = 0.2）#####################
#EC-AEA的下调配体气泡图（logFC = 0.2）
##EC-AEA作为source
AEA.s.ldown <- ligand.down %>%  
  filter(source == "EC-AEA")
AEA.s.ldown <- AEA.s.ldown %>%  
  arrange(substr(pathway_name, 1, 1)) 

unique(AEA.s.ldown$pathway_name)
AEA.s.ldown.use = AEA.s.ldown[, "interaction_name", drop = F]
AEA.s.ldown$interaction_name <- as.factor(AEA.s.ldown$interaction_name)
AEA.s.ldown$interaction_name <- fct_inorder(AEA.s.ldown$interaction_name)
AEA.s.ldown$interaction_name_2 <- as.factor(AEA.s.ldown$interaction_name_2)
AEA.s.ldown$interaction_name_2 <- fct_inorder(AEA.s.ldown$interaction_name_2)
AEA.s.ldown$pathway_name <- as.factor(AEA.s.ldown$pathway_name)
AEA.s.ldown$pathway_name <- fct_inorder(AEA.s.ldown$pathway_name)
AEA.s.ldown <- AEA.s.ldown %>% 
  mutate(p="")
AEA.s.ldown$signway <- paste(AEA.s.ldown$source,
                          "  ->  ", 
                          AEA.s.ldown$target) 

table(unique(AEA.s.ldown$interaction_name_2))
summary_data <- AEA.s.ldown %>%  
  group_by(pathway_name) %>%
  summarise(interaction_count = n_distinct(interaction_name_2))
print(summary_data)
write.csv(summary_data, "data.csv")

p1 <- ggplot(AEA.s.ldown,
             aes(x = signway,
                 y = interaction_name_2)) + 
  geom_point(aes(size = pval, 
                 color = ligand.logFC)) +
  geom_point(shape = 1,
             aes(size = pval), 
             color="black")+
  scale_size(rang = c(1.5,6)) +
  labs(x=NULL,y=NULL,size="P value",
       title = "Down-regulated signaling of EC-AEA in DKD") +
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.text =element_text(size = 8, color = "black"),
        axis.text.y = element_text(face="italic"),
        axis.text.x=element_text(angle = 90,
                                 hjust = 1,
                                 vjust = 0.5)) +
  scale_color_gradient2(low='navy',
                        high='#ADD8E6') +
  geom_hline(yintercept=c(1.5, 3.5, 4.5, 5.5, 8.5, 10.5, 12.5, 14.5, 
                          18.5))+
  theme(plot.margin = margin(t = 5, 
                             r = 5, 
                             b = 5, 
                             l = 5, 
                             unit = "pt"))
p1

labels_and_y <- list("CD46" = 1, "CD99" = 2.5, 
                     "ESAM" = 4, "GAS" = 5, 
                     "IGF" = 7, "ICAM" = 9.5, 
                     "MHC-I" = 11.5, "MHC-II" = 13.5, 
                     "NOTCH" = 16.5, "VEGF" = 20) 
gene_grodown <- ggplot(AEA.s.ldown,
                     aes(x = p,
                         y = interaction_name_2,
                         fill = pathway_name))+
  geom_tile() + 
  scale_y_discrete(position="left") +
  theme_minimal()+xlab(NULL) + ylab(NULL) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x =element_text(
          angle =90,hjust =0.5,vjust = 0.5),
        legend.position = 'none')+
  lapply(names(labels_and_y), function(label) {  
    y_val <- labels_and_y[[label]]  
    annotate("text", label=label, x=1.2, y=y_val)})+  
  theme(plot.margin = margin(t = 0, 
                             r = 0, 
                             b = 0, 
                             l = 0, 
                             unit = "pt"))
gene_grodown

p1 %>% insert_left(gene_grodown, width = 0.3)
  ##得到3.1 EC-AEA的source下调配体气泡图A（logFC = 0.2）


gg1 <- netVisual_bubble(cellchat, 
                        pairLR.use = AEA.s.ldown.use, 
                        sources.use = "EC-AEA", 
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "Down-regulated signaling of EC-AEA in DKD")
gg1 <- gg1 + 
  geom_hline(yintercept=c(1.5, 3.5, 4.5, 5.5, 8.5, 10.5, 12.5, 14.5, 
                          18.5))
gg1 %>% insert_left(gene_grodown, width = 0.3)
  ##得到3.2 EC-AEA的source下调配体气泡图B（logFC = 0.2）


##EC-AEA作为target
AEA.t.ldown <- AEA.t.ldown %>%  
  arrange(substr(pathway_name, 1, 1)) 

unique(AEA.t.ldown$pathway_name)
AEA.t.ldown.use = AEA.t.ldown[, "interaction_name", drop = F]
AEA.t.ldown$interaction_name <- as.factor(AEA.t.ldown$interaction_name)
AEA.t.ldown$interaction_name <- fct_inorder(AEA.t.ldown$interaction_name)
AEA.t.ldown$interaction_name_2 <- as.factor(AEA.t.ldown$interaction_name_2)
AEA.t.ldown$interaction_name_2 <- fct_inorder(AEA.t.ldown$interaction_name_2)
AEA.t.ldown$pathway_name <- as.factor(AEA.t.ldown$pathway_name)
AEA.t.ldown$pathway_name <- fct_inorder(AEA.t.ldown$pathway_name)
AEA.t.ldown <- AEA.t.ldown %>% 
  mutate(p="")

AEA.t.ldown$signway <- paste(AEA.t.ldown$source,
                          "  ->  ", 
                          AEA.t.ldown$target) 

table(unique(AEA.t.ldown$interaction_name_2))
summary_data <- AEA.t.ldown %>%  
  group_by(pathway_name) %>%
  summarise(interaction_count = n_distinct(interaction_name_2))
print(summary_data)
write.csv(summary_data, "data.csv")

p1 <- ggplot(AEA.t.ldown,
             aes(x = signway,
                 y = interaction_name_2)) + 
  geom_point(aes(size = pval, 
                 color = ligand.logFC)) +
  geom_point(shape = 1,
             aes(size = pval), 
             color="black")+
  scale_size(rang = c(1.5,6)) +
  labs(x=NULL,y=NULL,size="P value",
       title = "Down-regulated signaling of EC-AEA in DKD") +
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.text =element_text(size = 8, color = "black"),
        axis.text.y = element_text(face="italic"),
        axis.text.x=element_text(angle = 90,
                                 hjust = 1,
                                 vjust = 0.5)) +
  scale_color_gradient2(low='navy',
                        high='#ADD8E6') +
  geom_hline(yintercept=c(2.5, 3.5, 5.5, 6.5, 12.5, 15.5, 24.5, 25.5,
                          26.5, 28.5, 29.5, 30.5, 32.5, 34.5, 35.5, 
                          39.5, 40.5, 41.5, 44.5, 46.5, 48.5, 49.5, 
                          61.5, 66.5, 71.5, 75.5, 78.5, 79.5, 80.5, 
                          81.5, 82.5, 84.5, 89.5))+
  theme(plot.margin = margin(t = 5, 
                             r = 5, 
                             b = 5, 
                             l = 5, 
                             unit = "pt"))
p1

labels_and_y <- list("ANGPT" = 1.5, "APP" = 3, 
                     "BMP" = 4.5, "BTLA" = 6, 
                     "CCL" = 9.5, "COMPLEMENT" = 14, 
                     "COLLAGEN" = 20, "CD46" = 25, 
                     "CD48" = 26, "CD86" = 27.5, 
                     "CD99" = 29, "CDH5" = 30, 
                     "DHEAS" = 31.5, "EGF" = 33.5, 
                     "ESAM" = 35, "FGF" = 37.5, 
                     "GDF" = 40, "GAS" = 41, 
                     "IGF" = 43, "IL1" = 45.5, 
                     "ICAM" = 47.5, "KLK" = 49, 
                     "LAMININ" = 55.5, "MHC-I" = 64, 
                     "MHC-II" = 69, "NT" = 73.5, 
                     "NOTCH" = 77, "Netrin" = 79, 
                     "PDGF" = 80, "SEMA6" = 81, 
                     "THBS" = 82, "THY1" = 83.5, 
                     "VEGF" = 87, "VCAM" = 90)  
gene_grodown <- ggplot(AEA.t.ldown,
                     aes(x = p,
                         y = interaction_name_2,
                         fill = pathway_name))+
  geom_tile() + 
  scale_y_discrete(position="left") +
  theme_minimal()+xlab(NULL) + ylab(NULL) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x =element_text(
          angle =90,hjust =0.5,vjust = 0.5),
        legend.position = 'none')+
  lapply(names(labels_and_y), function(label) {  
    y_val <- labels_and_y[[label]]  
    annotate("text", label=label, x=1.2, y=y_val)})+
  theme(plot.margin = margin(t = 0, 
                             r = 0, 
                             b = 0, 
                             l = 0, 
                             unit = "pt"))
gene_grodown
p1 %>% insert_left(gene_grodown, width = 0.3)
  ##得到3.3 EC-AEA的target下调配体气泡图A（logFC = 0.2）

gg1 <- netVisual_bubble(cellchat, 
                        pairLR.use = AEA.t.ldown.use, 
                        targets.use = "EC-AEA", 
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "Down-regulated signaling of EC-AEA in DKD")
gg1 <- gg1 + 
  geom_hline(yintercept=c(2.5, 3.5, 5.5, 6.5, 12.5, 15.5, 24.5, 25.5,
                          26.5, 28.5, 29.5, 30.5, 32.5, 34.5, 35.5, 
                          39.5, 40.5, 41.5, 44.5, 46.5, 48.5, 49.5, 
                          61.5, 66.5, 71.5, 75.5, 78.5, 79.5, 80.5, 
                          81.5, 82.5, 84.5, 89.5))
gg1 %>% insert_left(gene_grodown, width = 0.3)
  ##得到3.4 EC-AEA的target下调配体气泡图B（logFC = 0.2）



####提取EC-AEA的下调受体（logFC = 0.2）#######################
#提取EC-AEA的下调受体（logFC = 0.2）
AEA.s.rdown <- receptor.down %>%  
  filter(source == "EC-AEA")
AEA.t.rdown <- receptor.down %>%  
  filter(target == "EC-AEA")



####EC-AEA的下调受体气泡图（logFC = 0.2）#####################
#EC-AEA的下调受体气泡图（logFC = 0.2）
##EC-AEA作为source
AEA.s.rdown <- AEA.s.rdown %>%  
  arrange(substr(pathway_name, 1, 1)) 

unique(AEA.s.rdown$pathway_name)
AEA.s.rdown.use = AEA.s.rdown[, "interaction_name", drop = F]
AEA.s.rdown$interaction_name <- as.factor(AEA.s.rdown$interaction_name)
AEA.s.rdown$interaction_name <- fct_inorder(AEA.s.rdown$interaction_name)
AEA.s.rdown$interaction_name_2 <- as.factor(AEA.s.rdown$interaction_name_2)
AEA.s.rdown$interaction_name_2 <- fct_inorder(AEA.s.rdown$interaction_name_2)
AEA.s.rdown$pathway_name <- as.factor(AEA.s.rdown$pathway_name)
AEA.s.rdown$pathway_name <- fct_inorder(AEA.s.rdown$pathway_name)
AEA.s.rdown <- AEA.s.rdown %>% 
  mutate(p="")
AEA.s.rdown$signway <- paste(AEA.s.rdown$source,
                          "  ->  ", 
                          AEA.s.rdown$target) 

table(unique(AEA.s.rdown$interaction_name_2))
summary_data <- AEA.s.rdown %>%  
  group_by(pathway_name) %>%
  summarise(interaction_count = n_distinct(interaction_name_2))
print(summary_data)
write.csv(summary_data, "data.csv")

p1 <- ggplot(AEA.s.rdown,
             aes(x = signway,
                 y = interaction_name_2)) + 
  geom_point(aes(size = pval, 
                 color = receptor.logFC)) +
  geom_point(shape = 1,
             aes(size = pval), 
             color="black")+
  scale_size(rang = c(1.5,6)) +
  labs(x=NULL,y=NULL,size="P value",
       title = "Down-regulated signaling of EC-AEA in DKD") +
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.text =element_text(size = 8, color = "black"),
        axis.text.y = element_text(face="italic"),
        axis.text.x=element_text(angle = 90,
                                 hjust = 1,
                                 vjust = 0.5)) +
  scale_color_gradient2(low='navy',
                        high='#ADD8E6') +
  geom_hline(yintercept=c(3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 24.5, 25.5, 
                          26.5, 27.5, 31.5, 33.5, 34.5, 38.5, 46.5, 
                          47.5, 49.5, 50.5, 51.5, 53.5, 54.5, 56.5, 
                          58.5, 59.5, 78.5, 79.5, 83.5, 86.5, 89.5, 
                          90.5, 93.5, 94.5, 95.5, 104.5, 105.5, 107.5, 
                          109.5))+
  theme(plot.margin = margin(t = 5, 
                             r = 5, 
                             b = 5, 
                             l = 5, 
                             unit = "pt"))
p1

labels_and_y <- list("ANGPT" = 2, "ApoA" = 4, 
                     "BAFF" = 5, "CXCL" = 6, 
                     "CD40" = 7, "CALCR" = 8, 
                     "COLLAGEN" = 16.5, "CD46" = 25, 
                     "CD99" = 26, "CDH5" = 27, 
                     "EPHA" = 29.5, "EPHB" = 32.5, 
                     "ESAM" = 34, "FGF" = 36.5, 
                     "FN1" = 42.5, "GH" = 47, 
                     "IGF" = 48.5, "IL4" = 50, 
                     "IL6" = 51, "IL1" = 52.5, 
                     "IFN-II" = 54, "JAM" = 55.5, 
                     "LIFR" = 57.5, "LT" = 59, 
                     "LAMININ" = 69, "MIF" = 79, 
                     "MK" = 81.5, "NT" = 85, 
                     "NOTCH" = 88, "OSM" = 90, 
                     "PDGF" = 92, "PARs" = 94, 
                     "Prostaglandin" = 95, "TGFb" = 100, 
                     "THBS" = 105, "TENASCIN" = 106.5, 
                     "VTN" = 108.5, "WNT" = 110) 
gene_grodown <- ggplot(AEA.s.rdown,
                     aes(x = p,
                         y = interaction_name_2,
                         fill = pathway_name))+
  geom_tile() + 
  scale_y_discrete(position="left") +
  theme_minimal()+xlab(NULL) + ylab(NULL) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x =element_text(
          angle =90,hjust =0.5,vjust = 0.5),
        legend.position = 'none')+
  lapply(names(labels_and_y), function(label) {  
    y_val <- labels_and_y[[label]]  
    annotate("text", label=label, x=1.2, y=y_val)})+  
  theme(plot.margin = margin(t = 0, 
                             r = 0, 
                             b = 0, 
                             l = 0, 
                             unit = "pt"))
gene_grodown

p1 %>% insert_left(gene_grodown, width = 0.3)
  ##得到4.1 EC-AEA的source下调受体气泡图A（logFC = 0.2）


gg1 <- netVisual_bubble(cellchat, 
                        pairLR.use = AEA.s.rdown.use, 
                        sources.use = "EC-AEA", 
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "Down-regulated signaling of EC-AEA in DKD")
gg1 <- gg1 + 
  geom_hline(yintercept=c(3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 24.5, 25.5, 
                          26.5, 27.5, 31.5, 33.5, 34.5, 38.5, 46.5, 
                          47.5, 49.5, 50.5, 51.5, 53.5, 54.5, 56.5, 
                          58.5, 59.5, 78.5, 79.5, 83.5, 86.5, 89.5, 
                          90.5, 93.5, 94.5, 95.5, 104.5, 105.5, 107.5, 
                          109.5))
gg1 %>% insert_left(gene_grodown, width = 0.3)
  ##得到4.2 EC-AEA的source下调受体气泡图B（logFC = 0.2）


##EC-AEA作为target
AEA.t.rdown <- AEA.t.rdown %>%  
  arrange(substr(pathway_name, 1, 1)) 

unique(AEA.t.rdown$pathway_name)
AEA.t.rdown.use = AEA.t.rdown[, "interaction_name", drop = F]
AEA.t.rdown$interaction_name <- as.factor(AEA.t.rdown$interaction_name)
AEA.t.rdown$interaction_name <- fct_inorder(AEA.t.rdown$interaction_name)
AEA.t.rdown$interaction_name_2 <- as.factor(AEA.t.rdown$interaction_name_2)
AEA.t.rdown$interaction_name_2 <- fct_inorder(AEA.t.rdown$interaction_name_2)
AEA.t.rdown$pathway_name <- as.factor(AEA.t.rdown$pathway_name)
AEA.t.rdown$pathway_name <- fct_inorder(AEA.t.rdown$pathway_name)
AEA.t.rdown <- AEA.t.rdown %>% 
  mutate(p="")

AEA.t.rdown$signway <- paste(AEA.t.rdown$source,
                          "  ->  ", 
                          AEA.t.rdown$target) 

table(unique(AEA.t.rdown$interaction_name_2))
summary_data <- AEA.t.rdown %>%  
  group_by(pathway_name) %>%
  summarise(interaction_count = n_distinct(interaction_name_2))
print(summary_data)
write.csv(summary_data, "data.csv")

p1 <- ggplot(AEA.t.rdown,
             aes(x = signway,
                 y = interaction_name_2)) + 
  geom_point(aes(size = pval, 
                 color = receptor.logFC)) +
  geom_point(shape = 1,
             aes(size = pval), 
             color="black")+
  scale_size(rang = c(1.5,6)) +
  labs(x=NULL,y=NULL,size="P value",
       title = "Down-regulated signaling of EC-AEA in DKD") +
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.text =element_text(size = 8, color = "black"),
        axis.text.y = element_text(face="italic"),
        axis.text.x=element_text(angle = 90,
                                 hjust = 1,
                                 vjust = 0.5)) +
  scale_color_gradient2(low='navy',
                        high='#ADD8E6') +
  geom_hline(yintercept=c(1.5, 2.5, 6.5, 9.5, 10.5, 16.5, 17.5, 
                          18.5, 19.5))+
  theme(plot.margin = margin(t = 5, 
                             r = 5, 
                             b = 5, 
                             l = 5, 
                             unit = "pt"))
p1

labels_and_y <- list("CD46" = 1, "CD99" = 2, 
                     "EPHA" = 4.5, "EPHB" = 8, 
                     "ESAM" = 10, "FGF" = 13.5, 
                     "LT" = 17, "MK" = 18, 
                     "PTN" = 19, "TWEAK" = 20)  
gene_grodown <- ggplot(AEA.t.rdown,
                     aes(x = p,
                         y = interaction_name_2,
                         fill = pathway_name))+
  geom_tile() + 
  scale_y_discrete(position="left") +
  theme_minimal()+xlab(NULL) + ylab(NULL) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x =element_text(
          angle =90,hjust =0.5,vjust = 0.5),
        legend.position = 'none')+
  lapply(names(labels_and_y), function(label) {  
    y_val <- labels_and_y[[label]]  
    annotate("text", label=label, x=1.2, y=y_val)})+
  theme(plot.margin = margin(t = 0, 
                             r = 0, 
                             b = 0, 
                             l = 0, 
                             unit = "pt"))
gene_grodown
p1 %>% insert_left(gene_grodown, width = 0.3)
  ##得到4.3 EC-AEA的target下调受体气泡图A（logFC = 0.2）

gg1 <- netVisual_bubble(cellchat, 
                        pairLR.use = AEA.t.rdown.use, 
                        targets.use = "EC-AEA", 
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "Down-regulated signaling of EC-AEA in DKD")
gg1 <- gg1 + 
  geom_hline(yintercept=c(1.5, 2.5, 6.5, 9.5, 10.5, 16.5, 17.5, 
                          18.5, 19.5))
gg1 %>% insert_left(gene_grodown, width = 0.3)
  ##得到4.4 EC-AEA的target下调受体气泡图B（logFC = 0.2）
