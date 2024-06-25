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



####提取EC-PTC的上调配体（logFC = 0.2）#######################
#提取EC-PTC的上调配体（logFC = 0.2）
PTC.s.lup <- ligand.up %>%  
  filter(source == "EC-PTC")
PTC.t.lup <- ligand.up %>%  
  filter(target == "EC-PTC")



####EC-PTC的上调配体气泡图（logFC = 0.2）#####################
#EC-PTC的上调配体气泡图（logFC = 0.2）
##EC-PTC作为source
unique(PTC.s.lup$pathway_name)
PTC.s.lup.use = PTC.s.lup[, "interaction_name", drop = F]
PTC.s.lup$interaction_name <- as.factor(PTC.s.lup$interaction_name)
PTC.s.lup$interaction_name <- fct_inorder(PTC.s.lup$interaction_name)
PTC.s.lup$interaction_name_2 <- as.factor(PTC.s.lup$interaction_name_2)
PTC.s.lup$interaction_name_2 <- fct_inorder(PTC.s.lup$interaction_name_2)
PTC.s.lup$pathway_name <- as.factor(PTC.s.lup$pathway_name)
PTC.s.lup$pathway_name <- fct_inorder(PTC.s.lup$pathway_name)
PTC.s.lup <- PTC.s.lup %>% 
  mutate(p="")
PTC.s.lup$signway <- paste(PTC.s.lup$source,
                          "  ->  ", 
                          PTC.s.lup$target) 

table(unique(PTC.s.lup$interaction_name_2))
summary_data <- PTC.s.lup %>%  
  group_by(pathway_name) %>%
  summarise(interaction_count = n_distinct(interaction_name_2))
print(summary_data)
write.csv(summary_data, "data.csv")

p1 <- ggplot(PTC.s.lup,
             aes(x = signway,
                 y = interaction_name_2)) + 
  geom_point(aes(size = pval, 
                 color = ligand.logFC)) +
  geom_point(shape = 1,
             aes(size = pval), 
             color="black")+
  scale_size(rang = c(1.5,6)) +
  labs(x=NULL,y=NULL,size="P value",
       title = "Up-regulated signaling of EC-PTC in DKD") +
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
  geom_hline(yintercept=c(1.5, 4.5, 5.5, 16.5, 17.5, 20.5,
                          21.5, 22.5))+
  theme(plot.margin = margin(t = 5, 
                             r = 5, 
                             b = 5, 
                             l = 5, 
                             unit = "pt"))
p1

labels_and_y <- list("CCL" = 1, "MIF" = 3, 
                     "TRAIL" = 5, "FN1" = 11, 
                     "CD34" = 17, "EPHB" = 19, 
                     "PECAM1" = 21, "VCAM" = 22, 
                     "ADGRG" = 23) 
gene_group <- ggplot(PTC.s.lup,
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
  ##得到1.1 EC-PTC的source上调配体气泡图A（logFC = 0.2）


gg1 <- netVisual_bubble(cellchat, 
                        pairLR.use = PTC.s.lup.use, 
                        sources.use = "EC-PTC", 
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "Up-regulated signaling of EC-PTC in DKD")
gg1 <- gg1 + 
  geom_hline(yintercept=c(1.5, 4.5, 5.5, 16.5, 17.5, 20.5,
                          21.5, 22.5))
gg1 %>% insert_left(gene_group, width = 0.3)
  ##得到1.2 EC-PTC的source上调配体气泡图B（logFC = 0.2）


##EC-PTC作为target
PTC.t.lup <- PTC.t.lup %>%  
  arrange(substr(pathway_name, 1, 1)) 

unique(PTC.t.lup$pathway_name)
PTC.t.lup.use = PTC.t.lup[, "interaction_name", drop = F]
PTC.t.lup$interaction_name <- as.factor(PTC.t.lup$interaction_name)
PTC.t.lup$interaction_name <- fct_inorder(PTC.t.lup$interaction_name)
PTC.t.lup$interaction_name_2 <- as.factor(PTC.t.lup$interaction_name_2)
PTC.t.lup$interaction_name_2 <- fct_inorder(PTC.t.lup$interaction_name_2)
PTC.t.lup$pathway_name <- as.factor(PTC.t.lup$pathway_name)
PTC.t.lup$pathway_name <- fct_inorder(PTC.t.lup$pathway_name)
PTC.t.lup <- PTC.t.lup %>% 
  mutate(p="")

PTC.t.lup$signway <- paste(PTC.t.lup$source,
                          "  ->  ", 
                          PTC.t.lup$target) 

table(unique(PTC.t.lup$interaction_name_2))
summary_data <- PTC.t.lup %>%  
  group_by(pathway_name) %>%
  summarise(interaction_count = n_distinct(interaction_name_2))
print(summary_data)
write.csv(summary_data, "data.csv")

p1 <- ggplot(PTC.t.lup,
             aes(x = signway,
                 y = interaction_name_2)) + 
  geom_point(aes(size = pval, 
                 color = ligand.logFC)) +
  geom_point(shape = 1,
             aes(size = pval), 
             color="black")+
  scale_size(rang = c(1.5,6)) +
  labs(x=NULL,y=NULL,size="P value",
       title = "Up-regulated signaling of EC-PTC in DKD") +
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
  geom_hline(yintercept=c(1.5, 2.5, 3.5, 6.5, 24.5, 25.5, 26.5, 
                          27.5, 28.5, 30.5, 32.5, 45.5, 48.5, 54.5,
                          55.5, 56.5, 62.5, 65.5, 69.5, 72.5, 78.5,
                          79.5, 80.5, 81.5, 84.5, 86.5, 87.5, 88.5,
                          94.5, 97.5, 98.5, 99.5, 100.5, 102.5))+
  theme(plot.margin = margin(t = 5, 
                             r = 5, 
                             b = 5, 
                             l = 5, 
                             unit = "pt"))
p1

labels_and_y <- list("APP" = 1, "CCL" = 2, 
                     "CXCL" = 3, "COMPLEMENT" = 5, 
                     "COLLAGEN" = 15.5, "CD34" = 25, 
                     "CD45" = 26, "CD46" = 27, 
                     "CD48" = 28, "CDH" = 29.5, 
                     "CDH1" = 31.5, "EPHA" = 39, 
                     "EPHB" = 47, "FN1" = 51.5, 
                     "GAS" = 55, "LT" = 56, 
                     "LAMININ" = 59.5, "MIF" = 64, 
                     "MK" = 67.5, "MHC-I" = 71, 
                     "MHC-II" = 75.5, "NRG" = 79, 
                     "Netrin" = 80, "OCLN" = 81, 
                     "PDGF" = 83, "PTN" = 85.5, 
                     "PSAP" = 87, "PECAM1" = 88, 
                     "SPP1" = 91.5, "SEMA3" = 96, 
                     "TRAIL" = 98, "TENASCIN" = 99, 
                     "THBS" = 100, "THY1" = 101.5, 
                     "VISTA" = 103)  
gene_group <- ggplot(PTC.t.lup,
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
  ##得到1.3 EC-PTC的target上调配体气泡图A（logFC = 0.2）

gg1 <- netVisual_bubble(cellchat, 
                        pairLR.use = PTC.t.lup.use, 
                        targets.use = "EC-PTC", 
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "Up-regulated signaling of EC-PTC in DKD")
gg1 <- gg1 + 
  geom_hline(yintercept=c(1.5, 2.5, 3.5, 6.5, 24.5, 25.5, 26.5, 
                          27.5, 28.5, 30.5, 32.5, 45.5, 48.5, 54.5,
                          55.5, 56.5, 62.5, 65.5, 69.5, 72.5, 78.5,
                          79.5, 80.5, 81.5, 84.5, 86.5, 87.5, 88.5,
                          94.5, 97.5, 98.5, 99.5, 100.5, 102.5))
gg1 %>% insert_left(gene_group, width = 0.3)
  ##得到1.4 EC-PTC的target上调配体气泡图B（logFC = 0.2）



####提取EC-PTC的上调受体（logFC = 0.2）#######################
#提取EC-PTC的上调受体（logFC = 0.2）
PTC.s.rup <- receptor.up %>%  
  filter(source == "EC-PTC")
PTC.t.rup <- receptor.up %>%  
  filter(target == "EC-PTC")



####EC-PTC的上调受体气泡图（logFC = 0.2）#####################
#EC-PTC的上调受体气泡图（logFC = 0.2）
##EC-PTC作为source
PTC.s.rup <- PTC.s.rup %>%  
  arrange(substr(pathway_name, 1, 1)) 

unique(PTC.s.rup$pathway_name)
PTC.s.rup.use = PTC.s.rup[, "interaction_name", drop = F]
PTC.s.rup$interaction_name <- as.factor(PTC.s.rup$interaction_name)
PTC.s.rup$interaction_name <- fct_inorder(PTC.s.rup$interaction_name)
PTC.s.rup$interaction_name_2 <- as.factor(PTC.s.rup$interaction_name_2)
PTC.s.rup$interaction_name_2 <- fct_inorder(PTC.s.rup$interaction_name_2)
PTC.s.rup$pathway_name <- as.factor(PTC.s.rup$pathway_name)
PTC.s.rup$pathway_name <- fct_inorder(PTC.s.rup$pathway_name)
PTC.s.rup <- PTC.s.rup %>% 
  mutate(p="")
PTC.s.rup$signway <- paste(PTC.s.rup$source,
                          "  ->  ", 
                          PTC.s.rup$target) 

table(unique(PTC.s.rup$interaction_name_2))
summary_data <- PTC.s.rup %>%  
  group_by(pathway_name) %>%
  summarise(interaction_count = n_distinct(interaction_name_2))
print(summary_data)
write.csv(summary_data, "data.csv")

p1 <- ggplot(PTC.s.rup,
             aes(x = signway,
                 y = interaction_name_2)) + 
  geom_point(aes(size = pval, 
                 color = receptor.logFC)) +
  geom_point(shape = 1,
             aes(size = pval), 
             color="black")+
  scale_size(rang = c(1.5,6)) +
  labs(x=NULL,y=NULL,size="P value",
       title = "Up-regulated signaling of EC-PTC in DKD") +
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
  geom_hline(yintercept=c(1.5, 2.5, 6.5, 9.5, 10.5, 11.5, 12.5, 30.5, 
                          31.5, 33.5, 34.5, 44.5, 45.5, 46.5, 49.5, 
                          50.5, 51.5, 52.5, 54.5, 61.5, 62.5, 63.5, 
                          65.5, 67.5, 68.5, 69.5, 70.5, 71.5, 72.5, 
                          73.5, 82.5, 85.5, 92.5, 95.5, 97.5))+
  theme(plot.margin = margin(t = 5, 
                             r = 5, 
                             b = 5, 
                             l = 5, 
                             unit = "pt"))
p1

labels_and_y <- list("ANGPT" = 1, "APP" = 2, 
                     "ADGRG" = 4.5, "BMP" = 8, 
                     "BAFF" = 10, "CD40" = 11, 
                     "CALCR" = 12, "COLLAGEN" = 21.5, 
                     "CD34" = 31, "CDH" = 32.5, 
                     "CNTN" = 34, "FN1" = 39.5, 
                     "GRN" = 45, "IGF" = 46, 
                     "IL2" = 48, "IL6" = 50, 
                     "IGFBP" = 51, "ICAM" = 52, 
                     "LIFR" = 53.5, "LAMININ" = 58, 
                     "L1CAM" = 62, "MIF" = 63, 
                     "NT" = 64.5, "NOTCH" = 66.5, 
                     "OSM" = 68, "OCLN" = 69, 
                     "PECAM1" = 70, "RELN" = 71, 
                     "SELPLG" = 72, "SIRP" = 73, 
                     "TGFb" = 78, "THBS" = 84, 
                     "VEGF" = 89, "VTN" = 94, 
                     "VCAM" = 96.5, "WNT" = 98.5) 
gene_group <- ggplot(PTC.s.rup,
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
  ##得到2.1 EC-PTC的source上调受体气泡图A（logFC = 0.2）


gg1 <- netVisual_bubble(cellchat, 
                        pairLR.use = PTC.s.rup.use, 
                        sources.use = "EC-PTC", 
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "Up-regulated signaling of EC-PTC in DKD")
gg1 <- gg1 + 
  geom_hline(yintercept=c(1.5, 2.5, 6.5, 9.5, 10.5, 11.5, 12.5, 30.5, 
                          31.5, 33.5, 34.5, 44.5, 45.5, 46.5, 49.5, 
                          50.5, 51.5, 52.5, 54.5, 61.5, 62.5, 63.5, 
                          65.5, 67.5, 68.5, 69.5, 70.5, 71.5, 72.5, 
                          73.5, 82.5, 85.5, 92.5, 95.5, 97.5))
gg1 %>% insert_left(gene_group, width = 0.3)
  ##得到2.2 EC-PTC的source上调受体气泡图B（logFC = 0.2）


##EC-PTC作为target
PTC.t.rup <- PTC.t.rup %>%  
  arrange(substr(pathway_name, 1, 1)) 

unique(PTC.t.rup$pathway_name)
PTC.t.rup.use = PTC.t.rup[, "interaction_name", drop = F]
PTC.t.rup$interaction_name <- as.factor(PTC.t.rup$interaction_name)
PTC.t.rup$interaction_name <- fct_inorder(PTC.t.rup$interaction_name)
PTC.t.rup$interaction_name_2 <- as.factor(PTC.t.rup$interaction_name_2)
PTC.t.rup$interaction_name_2 <- fct_inorder(PTC.t.rup$interaction_name_2)
PTC.t.rup$pathway_name <- as.factor(PTC.t.rup$pathway_name)
PTC.t.rup$pathway_name <- fct_inorder(PTC.t.rup$pathway_name)
PTC.t.rup <- PTC.t.rup %>% 
  mutate(p="")

PTC.t.rup$signway <- paste(PTC.t.rup$source,
                          "  ->  ", 
                          PTC.t.rup$target) 

table(unique(PTC.t.rup$interaction_name_2))
summary_data <- PTC.t.rup %>%  
  group_by(pathway_name) %>%
  summarise(interaction_count = n_distinct(interaction_name_2))
print(summary_data)
write.csv(summary_data, "data.csv")

p1 <- ggplot(PTC.t.rup,
             aes(x = signway,
                 y = interaction_name_2)) + 
  geom_point(aes(size = pval, 
                 color = receptor.logFC)) +
  geom_point(shape = 1,
             aes(size = pval), 
             color="black")+
  scale_size(rang = c(1.5,6)) +
  labs(x=NULL,y=NULL,size="P value",
       title = "Up-regulated signaling of EC-PTC in DKD") +
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
  geom_hline(yintercept=c(5.5, 6.5, 22.5, 23.5, 24.5, 27.5, 28.5, 
                          30.5, 31.5, 51.5, 53.5, 54.5, 55.5, 58.5, 
                          60.5, 61.5, 70.5, 72.5))+
  theme(plot.margin = margin(t = 5, 
                             r = 5, 
                             b = 5, 
                             l = 5, 
                             unit = "pt"))
p1

labels_and_y <- list("BMP" = 3, "CALCR" = 6, 
                     "COLLAGEN" = 14.5, "CD34" = 23, 
                     "CDH1" = 24, "FN1" = 26, 
                     "GDF" = 28, "IGF" = 29.5, 
                     "JAM" = 31, "LAMININ" = 41.5, 
                     "MK" = 52.5, "NRG" = 54, 
                     "PECAM1" = 55, "SPP1" = 57, 
                     "SELL" = 59.5, "SELPLG" = 61, 
                     "TGFb" = 66, "TENASCIN" = 71.5, 
                     "VTN" = 73.5)  
gene_group <- ggplot(PTC.t.rup,
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
  ##得到2.3 EC-PTC的target上调受体气泡图A（logFC = 0.2）

gg1 <- netVisual_bubble(cellchat, 
                        pairLR.use = PTC.t.rup.use, 
                        targets.use = "EC-PTC", 
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "Up-regulated signaling of EC-PTC in DKD")
gg1 <- gg1 + 
  geom_hline(yintercept=c(5.5, 6.5, 22.5, 23.5, 24.5, 27.5, 28.5, 
                          30.5, 31.5, 51.5, 53.5, 54.5, 55.5, 58.5, 
                          60.5, 61.5, 70.5, 72.5))
gg1 %>% insert_left(gene_group, width = 0.3)
  ##得到2.4 EC-PTC的target上调受体气泡图B（logFC = 0.2）



####提取EC-PTC的下调配体（logFC = 0.2）#######################
#提取EC-PTC的下调配体（logFC = 0.2）
PTC.s.ldown <- ligand.down %>%  
  filter(source == "EC-PTC")
PTC.t.ldown <- ligand.down %>%  
  filter(target == "EC-PTC")



####EC-PTC的下调配体气泡图（logFC = 0.2）#####################
#EC-PTC的下调配体气泡图（logFC = 0.2）
##EC-PTC作为source
PTC.s.ldown <- ligand.down %>%  
  filter(source == "EC-PTC")
PTC.s.ldown <- PTC.s.ldown %>%  
  arrange(substr(pathway_name, 1, 1)) 

unique(PTC.s.ldown$pathway_name)
PTC.s.ldown.use = PTC.s.ldown[, "interaction_name", drop = F]
PTC.s.ldown$interaction_name <- as.factor(PTC.s.ldown$interaction_name)
PTC.s.ldown$interaction_name <- fct_inorder(PTC.s.ldown$interaction_name)
PTC.s.ldown$interaction_name_2 <- as.factor(PTC.s.ldown$interaction_name_2)
PTC.s.ldown$interaction_name_2 <- fct_inorder(PTC.s.ldown$interaction_name_2)
PTC.s.ldown$pathway_name <- as.factor(PTC.s.ldown$pathway_name)
PTC.s.ldown$pathway_name <- fct_inorder(PTC.s.ldown$pathway_name)
PTC.s.ldown <- PTC.s.ldown %>% 
  mutate(p="")
PTC.s.ldown$signway <- paste(PTC.s.ldown$source,
                          "  ->  ", 
                          PTC.s.ldown$target) 

table(unique(PTC.s.ldown$interaction_name_2))
summary_data <- PTC.s.ldown %>%  
  group_by(pathway_name) %>%
  summarise(interaction_count = n_distinct(interaction_name_2))
print(summary_data)
write.csv(summary_data, "data.csv")

p1 <- ggplot(PTC.s.ldown,
             aes(x = signway,
                 y = interaction_name_2)) + 
  geom_point(aes(size = pval, 
                 color = ligand.logFC)) +
  geom_point(shape = 1,
             aes(size = pval), 
             color="black")+
  scale_size(rang = c(1.5,6)) +
  labs(x=NULL,y=NULL,size="P value",
       title = "Down-regulated signaling of EC-PTC in DKD") +
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
  geom_hline(yintercept=c(2.5, 3.5, 5.5))+
  theme(plot.margin = margin(t = 5, 
                             r = 5, 
                             b = 5, 
                             l = 5, 
                             unit = "pt"))
p1

labels_and_y <- list("APP" = 1.5, "GAS" = 3, 
                     "MHC-I" = 4.5, "MHC-II" = 6.5) 
gene_grodown <- ggplot(PTC.s.ldown,
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
  ##得到3.1 EC-PTC的source下调配体气泡图A（logFC = 0.2）


gg1 <- netVisual_bubble(cellchat, 
                        pairLR.use = PTC.s.ldown.use, 
                        sources.use = "EC-PTC", 
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "Down-regulated signaling of EC-PTC in DKD")
gg1 <- gg1 + 
  geom_hline(yintercept=c(2.5, 3.5, 5.5))
gg1 %>% insert_left(gene_grodown, width = 0.3)
  ##得到3.2 EC-PTC的source下调配体气泡图B（logFC = 0.2）


##EC-PTC作为target
PTC.t.ldown <- PTC.t.ldown %>%  
  arrange(substr(pathway_name, 1, 1)) 

unique(PTC.t.ldown$pathway_name)
PTC.t.ldown.use = PTC.t.ldown[, "interaction_name", drop = F]
PTC.t.ldown$interaction_name <- as.factor(PTC.t.ldown$interaction_name)
PTC.t.ldown$interaction_name <- fct_inorder(PTC.t.ldown$interaction_name)
PTC.t.ldown$interaction_name_2 <- as.factor(PTC.t.ldown$interaction_name_2)
PTC.t.ldown$interaction_name_2 <- fct_inorder(PTC.t.ldown$interaction_name_2)
PTC.t.ldown$pathway_name <- as.factor(PTC.t.ldown$pathway_name)
PTC.t.ldown$pathway_name <- fct_inorder(PTC.t.ldown$pathway_name)
PTC.t.ldown <- PTC.t.ldown %>% 
  mutate(p="")

PTC.t.ldown$signway <- paste(PTC.t.ldown$source,
                          "  ->  ", 
                          PTC.t.ldown$target) 

table(unique(PTC.t.ldown$interaction_name_2))
summary_data <- PTC.t.ldown %>%  
  group_by(pathway_name) %>%
  summarise(interaction_count = n_distinct(interaction_name_2))
print(summary_data)
write.csv(summary_data, "data.csv")

p1 <- ggplot(PTC.t.ldown,
             aes(x = signway,
                 y = interaction_name_2)) + 
  geom_point(aes(size = pval, 
                 color = ligand.logFC)) +
  geom_point(shape = 1,
             aes(size = pval), 
             color="black")+
  scale_size(rang = c(1.5,6)) +
  labs(x=NULL,y=NULL,size="P value",
       title = "Down-regulated signaling of EC-PTC in DKD") +
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
                          26.5, 27.5, 29.5, 30.5, 31.5, 33.5, 35.5, 
                          36.5, 39.5, 40.5, 41.5, 43.5, 45.5, 46.5,
                          48.5, 49.5, 61.5, 66.5, 71.5, 75.5, 78.5,
                          79.5, 80.5, 81.5, 83.5))+
  theme(plot.margin = margin(t = 5, 
                             r = 5, 
                             b = 5, 
                             l = 5, 
                             unit = "pt"))
p1

labels_and_y <- list("ANGPT" = 1.5, "APP" = 3, 
                     "BMP" = 4.5, "BTLA" = 6, 
                     "CCL" = 9.5, "COMPLEMENT" = 14, 
                     "COLLAGEN" = 20, "CD34" = 25, 
                     "CD46" = 26, "CD48" = 27, 
                     "CD86" = 28.5, "CD99" = 30, 
                     "CDH5" = 31, "DHEAS" = 32.5, 
                     "EGF" = 34.5, "ESAM" = 36, 
                     "FGF" = 38, "GDF" = 40, 
                     "GAS" = 41, "IGF" = 42.5, 
                     "IL1" = 44.5, "IGFBP" = 46, 
                     "ICAM" = 47.5, "KLK" = 49, 
                     "LAMININ" = 55.5, "MHC-I" = 64, 
                     "MHC-II" = 69, "NT" = 73.5, 
                     "NOTCH" = 77, "Netrin" = 79, 
                     "PDGF" = 80, "SEMA6" = 81, 
                     "THY1" = 82.5, "VEGF" = 86)  
gene_grodown <- ggplot(PTC.t.ldown,
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
  ##得到3.3 EC-PTC的target下调配体气泡图A（logFC = 0.2）

gg1 <- netVisual_bubble(cellchat, 
                        pairLR.use = PTC.t.ldown.use, 
                        targets.use = "EC-PTC", 
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "Down-regulated signaling of EC-PTC in DKD")
gg1 <- gg1 + 
  geom_hline(yintercept=c(2.5, 3.5, 5.5, 6.5, 12.5, 15.5, 24.5, 25.5,
                          26.5, 27.5, 29.5, 30.5, 31.5, 33.5, 35.5, 
                          36.5, 39.5, 40.5, 41.5, 43.5, 45.5, 46.5,
                          48.5, 49.5, 61.5, 66.5, 71.5, 75.5, 78.5,
                          79.5, 80.5, 81.5, 83.5))
gg1 %>% insert_left(gene_grodown, width = 0.3)
  ##得到3.4 EC-PTC的target下调配体气泡图B（logFC = 0.2）



####提取EC-PTC的下调受体（logFC = 0.2）#######################
#提取EC-PTC的下调受体（logFC = 0.2）
PTC.s.rdown <- receptor.down %>%  
  filter(source == "EC-PTC")
PTC.t.rdown <- receptor.down %>%  
  filter(target == "EC-PTC")



####EC-PTC的下调受体气泡图（logFC = 0.2）#####################
#EC-PTC的下调受体气泡图（logFC = 0.2）
##EC-PTC作为source
PTC.s.rdown <- PTC.s.rdown %>%  
  arrange(substr(pathway_name, 1, 1)) 

unique(PTC.s.rdown$pathway_name)
PTC.s.rdown.use = PTC.s.rdown[, "interaction_name", drop = F]
PTC.s.rdown$interaction_name <- as.factor(PTC.s.rdown$interaction_name)
PTC.s.rdown$interaction_name <- fct_inorder(PTC.s.rdown$interaction_name)
PTC.s.rdown$interaction_name_2 <- as.factor(PTC.s.rdown$interaction_name_2)
PTC.s.rdown$interaction_name_2 <- fct_inorder(PTC.s.rdown$interaction_name_2)
PTC.s.rdown$pathway_name <- as.factor(PTC.s.rdown$pathway_name)
PTC.s.rdown$pathway_name <- fct_inorder(PTC.s.rdown$pathway_name)
PTC.s.rdown <- PTC.s.rdown %>% 
  mutate(p="")
PTC.s.rdown$signway <- paste(PTC.s.rdown$source,
                          "  ->  ", 
                          PTC.s.rdown$target) 

table(unique(PTC.s.rdown$interaction_name_2))
summary_data <- PTC.s.rdown %>%  
  group_by(pathway_name) %>%
  summarise(interaction_count = n_distinct(interaction_name_2))
print(summary_data)
write.csv(summary_data, "data.csv")

p1 <- ggplot(PTC.s.rdown,
             aes(x = signway,
                 y = interaction_name_2)) + 
  geom_point(aes(size = pval, 
                 color = receptor.logFC)) +
  geom_point(shape = 1,
             aes(size = pval), 
             color="black")+
  scale_size(rang = c(1.5,6)) +
  labs(x=NULL,y=NULL,size="P value",
       title = "Down-regulated signaling of EC-PTC in DKD") +
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
  geom_hline(yintercept=c(3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 26.5, 27.5,
                          28.5, 29.5, 33.5, 36.5, 37.5, 39.5, 47.5,
                          48.5, 49.5, 50.5, 51.5, 53.5, 54.5, 56.5,
                          57.5, 63.5, 64.5, 67.5, 68.5, 69.5, 71.5,
                          72.5, 73.5, 74.5, 83.5, 85.5, 89.5))+
  theme(plot.margin = margin(t = 5, 
                             r = 5, 
                             b = 5, 
                             l = 5, 
                             unit = "pt"))
p1

labels_and_y <- list("ANGPT" = 2, "ApoA" = 4, 
                     "BAFF" = 5, "CXCL" = 6, 
                     "CD40" = 7, "CALCR" = 8, 
                     "COLLAGEN" = 17.5, "CD46" = 27, 
                     "CD99" = 28, "CDH5" = 29, 
                     "EPHA" = 31.5, "EPHB" = 35, 
                     "ESAM" = 37, "FGF" = 38.5, 
                     "FN1" = 43.5, "GH" = 48, 
                     "IGF" = 49, "IL4" = 50, 
                     "IL6" = 51, "IL1" = 52.5, 
                     "IFN-II" = 54, "LIFR" = 55.5, 
                     "LT" = 57, "LAMININ" = 60.5, 
                     "MIF" = 64, "NT" = 66, 
                     "NOTCH" = 68, "OSM" = 69, 
                     "PDGF" = 70.5, "PARs" = 72, 
                     "Prostaglandin" = 73, "SIRP" = 74, 
                     "TGFb" = 79, "THBS" = 84.5, 
                     "VTN" = 87.5, "WNT" = 90) 
gene_grodown <- ggplot(PTC.s.rdown,
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
  ##得到4.1 EC-PTC的source下调受体气泡图A（logFC = 0.2）


gg1 <- netVisual_bubble(cellchat, 
                        pairLR.use = PTC.s.rdown.use, 
                        sources.use = "EC-PTC", 
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "Down-regulated signaling of EC-PTC in DKD")
gg1 <- gg1 + 
  geom_hline(yintercept=c(3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 26.5, 27.5,
                          28.5, 29.5, 33.5, 36.5, 37.5, 39.5, 47.5,
                          48.5, 49.5, 50.5, 51.5, 53.5, 54.5, 56.5,
                          57.5, 63.5, 64.5, 67.5, 68.5, 69.5, 71.5,
                          72.5, 73.5, 74.5, 83.5, 85.5, 89.5))
gg1 %>% insert_left(gene_grodown, width = 0.3)
  ##得到4.2 EC-PTC的source下调受体气泡图B（logFC = 0.2）


##EC-PTC作为target
PTC.t.rdown <- PTC.t.rdown %>%  
  arrange(substr(pathway_name, 1, 1)) 

unique(PTC.t.rdown$pathway_name)
PTC.t.rdown.use = PTC.t.rdown[, "interaction_name", drop = F]
PTC.t.rdown$interaction_name <- as.factor(PTC.t.rdown$interaction_name)
PTC.t.rdown$interaction_name <- fct_inorder(PTC.t.rdown$interaction_name)
PTC.t.rdown$interaction_name_2 <- as.factor(PTC.t.rdown$interaction_name_2)
PTC.t.rdown$interaction_name_2 <- fct_inorder(PTC.t.rdown$interaction_name_2)
PTC.t.rdown$pathway_name <- as.factor(PTC.t.rdown$pathway_name)
PTC.t.rdown$pathway_name <- fct_inorder(PTC.t.rdown$pathway_name)
PTC.t.rdown <- PTC.t.rdown %>% 
  mutate(p="")

PTC.t.rdown$signway <- paste(PTC.t.rdown$source,
                          "  ->  ", 
                          PTC.t.rdown$target) 

table(unique(PTC.t.rdown$interaction_name_2))
summary_data <- PTC.t.rdown %>%  
  group_by(pathway_name) %>%
  summarise(interaction_count = n_distinct(interaction_name_2))
print(summary_data)
write.csv(summary_data, "data.csv")

p1 <- ggplot(PTC.t.rdown,
             aes(x = signway,
                 y = interaction_name_2)) + 
  geom_point(aes(size = pval, 
                 color = receptor.logFC)) +
  geom_point(shape = 1,
             aes(size = pval), 
             color="black")+
  scale_size(rang = c(1.5,6)) +
  labs(x=NULL,y=NULL,size="P value",
       title = "Down-regulated signaling of EC-PTC in DKD") +
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
  geom_hline(yintercept=c(6.5, 8.5, 9.5))+
  theme(plot.margin = margin(t = 5, 
                             r = 5, 
                             b = 5, 
                             l = 5, 
                             unit = "pt"))
p1

labels_and_y <- list("FGF" = 3.5, "LIFR" = 7.5, 
                     "OSM" = 9, "WNT" = 10)  
gene_grodown <- ggplot(PTC.t.rdown,
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
  ##得到4.3 EC-PTC的target下调受体气泡图A（logFC = 0.2）

gg1 <- netVisual_bubble(cellchat, 
                        pairLR.use = PTC.t.rdown.use, 
                        targets.use = "EC-PTC", 
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "Down-regulated signaling of EC-PTC in DKD")
gg1 <- gg1 + 
  geom_hline(yintercept=c(6.5, 8.5, 9.5))
gg1 %>% insert_left(gene_grodown, width = 0.3)
  ##得到4.4 EC-PTC的target下调受体气泡图B（logFC = 0.2）


