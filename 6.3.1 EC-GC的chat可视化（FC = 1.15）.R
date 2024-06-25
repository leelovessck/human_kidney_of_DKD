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



####提取EC-GC的上调配体（logFC = 0.2）#######################
#提取EC-GC的上调配体（logFC = 0.2）
gc.s.lup <- ligand.up %>%  
  filter(source == "EC-GC")
gc.t.lup <- ligand.up %>%  
  filter(target == "EC-GC")



####EC-GC的上调配体气泡图（logFC = 0.2）#####################
#EC-GC的上调配体气泡图（logFC = 0.2）
##EC-GC作为source
unique(gc.s.lup$pathway_name)
gc.s.lup.use = gc.s.lup[, "interaction_name", drop = F]
gc.s.lup$interaction_name <- as.factor(gc.s.lup$interaction_name)
gc.s.lup$interaction_name <- fct_inorder(gc.s.lup$interaction_name)
gc.s.lup$interaction_name_2 <- as.factor(gc.s.lup$interaction_name_2)
gc.s.lup$interaction_name_2 <- fct_inorder(gc.s.lup$interaction_name_2)
gc.s.lup$pathway_name <- as.factor(gc.s.lup$pathway_name)
gc.s.lup$pathway_name <- fct_inorder(gc.s.lup$pathway_name)
gc.s.lup <- gc.s.lup %>% 
  mutate(p="")
gc.s.lup$signway <- paste(gc.s.lup$source,
                          "  ->  ", 
                          gc.s.lup$target) 

table(unique(gc.s.lup$interaction_name_2))
summary_data <- gc.s.lup %>%  
  group_by(pathway_name) %>%
  summarise(interaction_count = n_distinct(interaction_name_2))
print(summary_data)

p1 <- ggplot(gc.s.lup,
             aes(x = signway,
                 y = interaction_name_2)) + 
  geom_point(aes(size = pval, 
                 color = ligand.logFC)) +
  geom_point(shape = 1,
             aes(size = pval), 
             color="black")+
  scale_size(rang = c(1.5,6)) +
  labs(x=NULL,y=NULL,size="P value",
       title = "Up-regulated signaling of EC-GC in DKD") +
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.text =element_text(size = 8, color = "black"),
        axis.text.y = element_text(face="italic"),
        axis.text.x=element_text(angle = 90,
                                 hjust = 1,
                                 vjust = 0.5)) +
  scale_color_gradient2(low = "#FFB6C1",
                        high = "#FF0000",
                        midpoint = 0) +
  geom_hline(yintercept=c(2.5, 5.5, 6.5))+
  theme(plot.margin = margin(t = 5, 
                             r = 5, 
                             b = 5, 
                             l = 5, 
                             unit = "pt"))
p1

labels_and_y <- list("MIF" = 1.5,  
                     "EPHB" = 4,  
                     "MHC-I" = 6,  
                     "PECAM1" = 7) 
gene_group <- ggplot(gc.s.lup,
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
  ##得到1.1 EC-GC的source上调配体气泡图A（logFC = 0.2）


gg1 <- netVisual_bubble(cellchat, 
                        pairLR.use = gc.s.lup.use, 
                        sources.use = "EC-GC", 
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "Up-regulated signaling of EC-GC in DKD")
gg1 <- gg1 + 
  geom_hline(yintercept=c(2.5, 5.5, 6.5))
gg1 %>% insert_left(gene_group, width = 0.3)
  ##得到1.2 EC-GC的source上调配体气泡图B（logFC = 0.2）


##EC-GC作为target
order_annotation <- c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact")  
gc.t.lup <- gc.t.lup %>%  
  mutate(sort_key = case_when(  
    annotation %in% order_annotation ~ factor(annotation, 
                                              levels = order_annotation),  
    TRUE ~ as.factor(annotation)  
  )) %>%  
  arrange(sort_key, substr(pathway_name, 1, 1)) %>%  
  select(-sort_key)  

unique(gc.t.lup$pathway_name)
gc.t.lup.use = gc.t.lup[, "interaction_name", drop = F]
gc.t.lup$interaction_name <- at.factor(g.t.lup$interaction_name)
gc.t.lup$interaction_name <- fct_inorder(gc.t.lup$interaction_name)
gc.t.lup$interaction_name_2 <- at.factor(g.t.lup$interaction_name_2)
gc.t.lup$interaction_name_2 <- fct_inorder(gc.t.lup$interaction_name_2)
gc.t.lup$pathway_name <- at.factor(g.t.lup$pathway_name)
gc.t.lup$pathway_name <- fct_inorder(gc.t.lup$pathway_name)
gc.t.lup <- gc.t.lup %>% 
  mutate(p="")

gc.t.lup$signway <- paste(gc.t.lup$source,
                          "  ->  ", 
                          gc.t.lup$target) 

table(unique(gc.t.lup$interaction_name_2))
summary_data <- gc.t.lup %>%  
  group_by(pathway_name) %>%
  summarise(interaction_count = n_distinct(interaction_name_2))
print(summary_data)

p1 <- ggplot(gc.t.lup,
             aes(x = signway,
                 y = interaction_name_2)) + 
  geom_point(aes(size = pval, 
                 color = ligand.logFC)) +
  geom_point(shape = 1,
             aes(size = pval), 
             color="black")+
  scale_size(rang = c(1.5,6)) +
  labs(x=NULL,y=NULL,size="P value",
       title = "Up-regulated signaling of EC-GC in DKD") +
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.text =element_text(size = 8, color = "black"),
        axis.text.y = element_text(face="italic"),
        axis.text.x=element_text(angle = 90,
                                 hjust = 1,
                                 vjust = 0.5)) +
  scale_color_gradient2(low="#FFB6C1",
                        high="#FF0000") +
  geom_hline(yintercept=c(1.5, 2.5, 5.5, 6.5, 7.5, 9.5, 13.5, 14.5, 
                          17.5, 19.5, 25.5, 28.5, 29.5, 47.5, 53.5,
                          59.5, 60.5, 61.5, 62.5, 63.5, 64.5, 66.5,
                          68.5, 77.5, 80.5, 83.5, 89.5, 90.5, 91.5,
                          92.5, 94.5))+
  theme(plot.margin = margin(t = 5, 
                             r = 5, 
                             b = 5, 
                             l = 5, 
                             unit = "pt"))
p1

labels_and_y <- list("CCL" = 1, "CXCL" = 2,  
                     "COMPLEMENT" = 4, "GAS" = 6,  
                     "LT" = 7, "MIF" = 8.5,
                     "MK" = 11.5, "NRG" = 14,
                     "PDGF" = 16, "PTN" = 18.5,
                     "SPP1" = 22.5, "SEMA3" = 27,
                     "TRAIL" = 29, "COLLAGEN" = 39.5,
                     "FN1" = 50.5, "LAMININ" = 56.5,
                     "TENASCIN" = 60, "APP" = 61,
                     "CD45" = 62, "CD46" = 63,
                     "CD48" = 64, "CDH" = 65.5,
                     "CDH1" = 67.5, "EPHA" = 74,
                     "EPHB" = 79, "MHC-I" = 82,
                     "MHC-II" = 87.5,
                     "Netrin" = 90,
                     "OCLN" = 91,
                     "PECAM1" = 92,
                     "THY1" = 93.5,
                     "VISTA" = 95)  
gene_group <- ggplot(gc.t.lup,
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
  ##得到1.3 EC-GC的target上调配体气泡图A（logFC = 0.2）

gg1 <- netVisual_bubble(cellchat, 
                        pairLR.use = gc.t.lup.use, 
                        targets.use = "EC-GC", 
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "Up-regulated signaling of EC-GC in DKD")
gg1 <- gg1 + 
  geom_hline(yintercept=c(1.5, 2.5, 5.5, 6.5, 7.5, 9.5, 13.5, 14.5, 
                          17.5, 19.5, 25.5, 28.5, 29.5, 47.5, 53.5,
                          59.5, 60.5, 61.5, 62.5, 63.5, 64.5, 66.5,
                          68.5, 77.5, 80.5, 83.5, 89.5, 90.5, 91.5,
                          92.5, 94.5))
gg1 %>% insert_left(gene_group, width = 0.3)
  ##得到1.4 EC-GC的target上调配体气泡图B（logFC = 0.2）



####提取EC-GC的上调受体（logFC = 0.2）#######################
#提取EC-GC的上调受体（logFC = 0.2）
gc.s.rup <- receptor.up %>%  
  filter(source == "EC-GC")
gc.t.rup <- receptor.up %>%  
  filter(target == "EC-GC")



####EC-GC的上调受体气泡图（logFC = 0.2）#####################
#EC-GC的上调受体气泡图（logFC = 0.2）
##EC-GC作为source
gc.s.rup <- gc.s.rup %>%  
  arrange(substr(pathway_name, 1, 1)) 

unique(gc.s.rup$pathway_name)
gc.s.rup.use = gc.s.rup[, "interaction_name", drop = F]
gc.s.rup$interaction_name <- as.factor(gc.s.rup$interaction_name)
gc.s.rup$interaction_name <- fct_inorder(gc.s.rup$interaction_name)
gc.s.rup$interaction_name_2 <- as.factor(gc.s.rup$interaction_name_2)
gc.s.rup$interaction_name_2 <- fct_inorder(gc.s.rup$interaction_name_2)
gc.s.rup$pathway_name <- as.factor(gc.s.rup$pathway_name)
gc.s.rup$pathway_name <- fct_inorder(gc.s.rup$pathway_name)
gc.s.rup <- gc.s.rup %>% 
  mutate(p="")
gc.s.rup$signway <- paste(gc.s.rup$source,
                          "  ->  ", 
                          gc.s.rup$target) 

table(unique(gc.s.rup$interaction_name_2))
summary_data <- gc.s.rup %>%  
  group_by(pathway_name) %>%
  summarise(interaction_count = n_distinct(interaction_name_2))
print(summary_data)

p1 <- ggplot(gc.s.rup,
             aes(x = signway,
                 y = interaction_name_2)) + 
  geom_point(aes(size = pval, 
                 color = receptor.logFC)) +
  geom_point(shape = 1,
             aes(size = pval), 
             color="black")+
  scale_size(rang = c(1.5,6)) +
  labs(x=NULL,y=NULL,size="P value",
       title = "Up-regulated signaling of EC-GC in DKD") +
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.text =element_text(size = 8, color = "black"),
        axis.text.y = element_text(face="italic"),
        axis.text.x=element_text(angle = 90,
                                 hjust = 1,
                                 vjust = 0.5)) +
  scale_color_gradient2(low = "#FFB6C1",
                        high = "#FF0000") +
  geom_hline(yintercept=c(1.5, 2.5, 6.5, 9.5, 10.5, 11.5, 12.5, 13.5,
                          30.5, 31.5, 33.5, 34.5, 43.5, 44.5, 45.5,
                          48.5, 49.5, 50.5, 51.5, 53.5, 60.5, 61.5,
                          63.5, 65.5, 66.5, 67.5, 68.5, 69.5, 70.5,
                          79.5, 82.5, 89.5, 91.5))+
  theme(plot.margin = margin(t = 5, 
                             r = 5, 
                             b = 5, 
                             l = 5, 
                             unit = "pt"))
p1

labels_and_y <- list("ANGPT" = 1, "APP" = 2,  
                     "ADGRG" = 4.5, "BMP" = 8,
                     "BAFF" = 10, "CX3C" = 11,
                     "CD40" = 12, "CALCR" = 13,  
                     "COLLAGEN" = 22, "CD34" = 31,
                     "CDH" = 32.5, "CNTN" = 34,  
                     "FN1" = 39, "GRN" = 44,
                     "IGF" = 45, "IL2" = 47,
                     "IL6" = 49, "IGFBP" = 50,
                     "ICAM" = 51, "LIFR" = 52.5,
                     "LAMININ" = 57, "L1CAM" = 61,
                     "NT" = 62.5, "NOTCH" = 64.5,
                     "OSM" = 66, "OCLN" = 67,
                     "PECAM1" = 68, "RELN" = 69,
                     "SIRP" = 70, "TGFb" = 75,
                     "THBS" = 81, "VEGF" = 86,
                     "VCAM" = 90.5, "WNT" = 92.5) 
gene_group <- ggplot(gc.s.rup,
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
  ##得到2.1 EC-GC的source上调受体气泡图A（logFC = 0.2）


gg1 <- netVisual_bubble(cellchat, 
                        pairLR.use = gc.s.rup.use, 
                        sources.use = "EC-GC", 
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "Up-regulated signaling of EC-GC in DKD")
gg1 <- gg1 + 
  geom_hline(yintercept=c(1.5, 2.5, 6.5, 9.5, 10.5, 11.5, 12.5, 13.5,
                          30.5, 31.5, 33.5, 34.5, 43.5, 44.5, 45.5,
                          48.5, 49.5, 50.5, 51.5, 53.5, 60.5, 61.5,
                          63.5, 65.5, 66.5, 67.5, 68.5, 69.5, 70.5,
                          79.5, 82.5, 89.5, 91.5))
gg1 %>% insert_left(gene_group, width = 0.3)
  ##得到2.2 EC-GC的source上调受体气泡图B（logFC = 0.2）


##EC-GC作为target
gc.t.rup <- gc.t.rup %>%  
  arrange(substr(pathway_name, 1, 1)) 

unique(gc.t.rup$pathway_name)
gc.t.rup.use = gc.t.rup[, "interaction_name", drop = F]
gc.t.rup$interaction_name <- as.factor(gc.t.rup$interaction_name)
gc.t.rup$interaction_name <- fct_inorder(gc.t.rup$interaction_name)
gc.t.rup$interaction_name_2 <- as.factor(gc.t.rup$interaction_name_2)
gc.t.rup$interaction_name_2 <- fct_inorder(gc.t.rup$interaction_name_2)
gc.t.rup$pathway_name <- as.factor(gc.t.rup$pathway_name)
gc.t.rup$pathway_name <- fct_inorder(gc.t.rup$pathway_name)
gc.t.rup <- gc.t.rup %>% 
  mutate(p="")

gc.t.rup$signway <- paste(gc.t.rup$source,
                          "  ->  ", 
                          gc.t.rup$target) 

table(unique(gc.t.rup$interaction_name_2))
summary_data <- gc.t.rup %>%  
  group_by(pathway_name) %>%
  summarise(interaction_count = n_distinct(interaction_name_2))
print(summary_data)

p1 <- ggplot(gc.t.rup,
             aes(x = signway,
                 y = interaction_name_2)) + 
  geom_point(aes(size = pval, 
                 color = receptor.logFC)) +
  geom_point(shape = 1,
             aes(size = pval), 
             color="black")+
  scale_size(rang = c(1.5,6)) +
  labs(x=NULL,y=NULL,size="P value",
       title = "Up-regulated signaling of EC-GC in DKD") +
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
  geom_hline(yintercept=c(1.5, 4.5, 5.5))+
  theme(plot.margin = margin(t = 5, 
                             r = 5, 
                             b = 5, 
                             l = 5, 
                             unit = "pt"))
p1

labels_and_y <- list("CNTN" = 1, "NOTCH" = 3,  
                     "PECAM1" = 5, "VEGF" = 7)  
gene_group <- ggplot(gc.t.rup,
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
  ##得到2.3 EC-GC的target上调受体气泡图A（logFC = 0.2）

gg1 <- netVisual_bubble(cellchat, 
                        pairLR.use = gc.t.rup.use, 
                        targets.use = "EC-GC", 
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "Up-regulated signaling of EC-GC in DKD")
gg1 <- gg1 + 
  geom_hline(yintercept=c(1.5, 4.5, 5.5))
gg1 %>% insert_left(gene_group, width = 0.3)
  ##得到2.4 EC-GC的target上调受体气泡图B（logFC = 0.2）



####提取EC-GC的下调配体（logFC = 0.2）#######################
#提取EC-GC的下调配体（logFC = 0.2）
gc.s.ldown <- ligand.down %>%  
  filter(source == "EC-GC")
gc.t.ldown <- ligand.down %>%  
  filter(target == "EC-GC")



####EC-GC的下调配体气泡图（logFC = 0.2）#####################
#EC-GC的下调配体气泡图（logFC = 0.2）
##EC-GC作为source
gc.s.ldown <- ligand.down %>%  
  filter(source == "EC-GC")
gc.s.ldown <- gc.s.ldown %>%  
  arrange(substr(pathway_name, 1, 1)) 
rows_to_move <- gc.s.ldown[9:16, ]  
remaining_rows <- gc.s.ldown[-(9:16), ]  
gc.s.ldown <- rbind(rows_to_move, remaining_rows)  
 

unique(gc.s.ldown$pathway_name)
gc.s.ldown.use = gc.s.ldown[, "interaction_name", drop = F]
gc.s.ldown$interaction_name <- as.factor(gc.s.ldown$interaction_name)
gc.s.ldown$interaction_name <- fct_inorder(gc.s.ldown$interaction_name)
gc.s.ldown$interaction_name_2 <- as.factor(gc.s.ldown$interaction_name_2)
gc.s.ldown$interaction_name_2 <- fct_inorder(gc.s.ldown$interaction_name_2)
gc.s.ldown$pathway_name <- as.factor(gc.s.ldown$pathway_name)
gc.s.ldown$pathway_name <- fct_inorder(gc.s.ldown$pathway_name)
gc.s.ldown <- gc.s.ldown %>% 
  mutate(p="")
gc.s.ldown$signway <- paste(gc.s.ldown$source,
                          "  ->  ", 
                          gc.s.ldown$target) 

table(unique(gc.s.ldown$interaction_name_2))
summary_data <- gc.s.ldown %>%  
  group_by(pathway_name) %>%
  summarise(interaction_count = n_distinct(interaction_name_2))
print(summary_data)

p1 <- ggplot(gc.s.ldown,
             aes(x = signway,
                 y = interaction_name_2)) + 
  geom_point(aes(size = pval, 
                 color = ligand.logFC)) +
  geom_point(shape = 1,
             aes(size = pval), 
             color="black")+
  scale_size(rang = c(1.5,6)) +
  labs(x=NULL,y=NULL,size="P value",
       title = "Down-regulated signaling of EC-GC in DKD") +
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.text =element_text(size = 8, color = "black"),
        axis.text.y = element_text(face="italic"),
        axis.text.x=element_text(angle = 90,
                                 hjust = 1,
                                 vjust = 0.5)) +
  scale_color_gradient2(low="#000080",
                        high="#ADD8E6") +
  geom_hline(yintercept=c(3.5, 5.5, 6.5, 7.5, 21.5, 22.5, 24.5,
                          25.5, 26.5, 27.5, 28.5, 30.5, 38.5,
                          39.5, 43.5, 44.5, 45.5, 46.5))+
  theme(plot.margin = margin(t = 5, 
                             r = 5, 
                             b = 5, 
                             l = 5, 
                             unit = "pt"))
p1

labels_and_y <- list("APP" = 4.5, "ADGRG" = 2,
                     "BTLA" = 6, "CX3C" = 7,
                     "COLLAGEN" = 15.5, "CD34" = 22,  
                     "CD99" = 23.5, "CDH5" = 25,
                     "ESAM" = 26, "GAS" = 27,  
                     "GRN" = 28, "ICAM" = 29.5,
                     "LAMININ" = 34.5, "MHC-I" = 39,  
                     "MHC-II" = 41.5, "Netrin" = 44,
                     "PDGF" = 45, "SEMA6" = 46,
                     "VCAM" = 47) 
gene_grodown <- ggplot(gc.s.ldown,
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
  ##得到3.1 EC-GC的source下调配体气泡图A（logFC = 0.2）


gg1 <- netVisual_bubble(cellchat, 
                        pairLR.use = gc.s.ldown.use, 
                        sources.use = "EC-GC", 
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "Down-regulated signaling of EC-GC in DKD")
gg1 <- gg1 + 
  geom_hline(yintercept=c(3.5, 5.5, 6.5, 7.5, 21.5, 22.5, 24.5,
                          25.5, 26.5, 27.5, 28.5, 30.5, 38.5,
                          39.5, 43.5, 44.5, 45.5, 46.5))
gg1 %>% insert_left(gene_grodown, width = 0.3)
  ##得到3.2 EC-GC的source下调配体气泡图B（logFC = 0.2）


##EC-GC作为target
gc.t.ldown <- gc.t.ldown %>%  
  arrange(substr(pathway_name, 1, 1)) 

unique(gc.t.ldown$pathway_name)
gc.t.ldown.use = gc.t.ldown[, "interaction_name", drop = F]
gc.t.ldown$interaction_name <- as.factor(gc.t.ldown$interaction_name)
gc.t.ldown$interaction_name <- fct_inorder(gc.t.ldown$interaction_name)
gc.t.ldown$interaction_name_2 <- as.factor(gc.t.ldown$interaction_name_2)
gc.t.ldown$interaction_name_2 <- fct_inorder(gc.t.ldown$interaction_name_2)
gc.t.ldown$pathway_name <- as.factor(gc.t.ldown$pathway_name)
gc.t.ldown$pathway_name <- fct_inorder(gc.t.ldown$pathway_name)
gc.t.ldown <- gc.t.ldown %>% 
  mutate(p="")

gc.t.ldown$signway <- paste(gc.t.ldown$source,
                          "  ->  ", 
                          gc.t.ldown$target) 

table(unique(gc.t.ldown$interaction_name_2))
summary_data <- gc.t.ldown %>%  
  group_by(pathway_name) %>%
  summarise(interaction_count = n_distinct(interaction_name_2))
print(summary_data)

p1 <- ggplot(gc.t.ldown,
             aes(x = signway,
                 y = interaction_name_2)) + 
  geom_point(aes(size = pval, 
                 color = ligand.logFC)) +
  geom_point(shape = 1,
             aes(size = pval), 
             color="black")+
  scale_size(rang = c(1.5,6)) +
  labs(x=NULL,y=NULL,size="P value",
       title = "Down-regulated signaling of EC-GC in DKD") +
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.text =element_text(size = 8, color = "black"),
        axis.text.y = element_text(face="italic"),
        axis.text.x=element_text(angle = 90,
                                 hjust = 1,
                                 vjust = 0.5)) +
  scale_color_gradient2(low='navy',
                        high="#ADD8E6") +
  geom_hline(yintercept=c(2.5, 3.5, 5.5, 6.5, 12.5, 13.5, 16.5,
                          25.5, 26.5, 27.5, 29.5, 30.5, 31.5, 33.5,
                          35.5, 36.5, 40.5, 41.5, 42.5, 44.5, 46.5,
                          47.5, 59.5, 64.5, 69.5, 73.5, 76.5, 77.5,
                          78.5, 79.5, 81.5))+
  theme(plot.margin = margin(t = 5, 
                             r = 5, 
                             b = 5, 
                             l = 5, 
                             unit = "pt"))
p1

labels_and_y <- list("ANGPT" = 1.5, "APP" = 3,  
                     "BMP" = 4.5, "BTLA" = 6, 
                     "CCL" = 9.5, "CX3C" = 13, 
                     "COMPLEMENT" = 15, "COLLAGEN" = 21, 
                     "CD46" = 26, "CD48" = 27, 
                     "CD86" = 28.5, "CD99" = 30, 
                     "CDH5" = 31, "DHEAS" = 32.5,
                     "EGF" = 34.5, "ESAM" = 36, 
                     "FGF" = 38.5, "GDF" = 41, 
                     "GAS" = 42, "IL1" = 43.5, 
                     "ICAM" = 45.5, "KLK" = 47, 
                     "LAMININ" = 53.5, "MHC-I" = 62, 
                     "MHC-II" = 67, "NT" = 71.5, 
                     "NOTCH" = 75, "Netrin" = 77, 
                     "PDGF" = 78, "SEMA6" = 79, 
                     "THY1" = 80.5, "VEGF" = 84)  
gene_grodown <- ggplot(gc.t.ldown,
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
  ##得到3.3 EC-GC的target下调配体气泡图A（logFC = 0.2）

gg1 <- netVisual_bubble(cellchat, 
                        pairLR.use = gc.t.ldown.use, 
                        targets.use = "EC-GC", 
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "Down-regulated signaling of EC-GC in DKD")
gg1 <- gg1 + 
  geom_hline(yintercept=c(2.5, 3.5, 5.5, 6.5, 12.5, 13.5, 16.5,
                          25.5, 26.5, 27.5, 29.5, 30.5, 31.5, 33.5,
                          35.5, 36.5, 40.5, 41.5, 42.5, 44.5, 46.5,
                          47.5, 59.5, 64.5, 69.5, 73.5, 76.5, 77.5,
                          78.5, 79.5, 81.5))
gg1 %>% insert_left(gene_grodown, width = 0.3)
  ##得到3.4 EC-GC的target下调配体气泡图B（logFC = 0.2）



####提取EC-GC的下调受体（logFC = 0.2）#######################
#提取EC-GC的下调受体（logFC = 0.2）
gc.s.rdown <- receptor.down %>%  
  filter(source == "EC-GC")
gc.t.rdown <- receptor.down %>%  
  filter(target == "EC-GC")



####EC-GC的下调受体气泡图（logFC = 0.2）#####################
#EC-GC的下调受体气泡图（logFC = 0.2）
##EC-GC作为source
gc.s.rdown <- gc.s.rdown %>%  
  arrange(substr(pathway_name, 1, 1)) 

unique(gc.s.rdown$pathway_name)
gc.s.rdown.use = gc.s.rdown[, "interaction_name", drop = F]
gc.s.rdown$interaction_name <- as.factor(gc.s.rdown$interaction_name)
gc.s.rdown$interaction_name <- fct_inorder(gc.s.rdown$interaction_name)
gc.s.rdown$interaction_name_2 <- as.factor(gc.s.rdown$interaction_name_2)
gc.s.rdown$interaction_name_2 <- fct_inorder(gc.s.rdown$interaction_name_2)
gc.s.rdown$pathway_name <- as.factor(gc.s.rdown$pathway_name)
gc.s.rdown$pathway_name <- fct_inorder(gc.s.rdown$pathway_name)
gc.s.rdown <- gc.s.rdown %>% 
  mutate(p="")
gc.s.rdown$signway <- paste(gc.s.rdown$source,
                          "  ->  ", 
                          gc.s.rdown$target) 

table(unique(gc.s.rdown$interaction_name_2))
summary_data <- gc.s.rdown %>%  
  group_by(pathway_name) %>%
  summarise(interaction_count = n_distinct(interaction_name_2))
print(summary_data)

p1 <- ggplot(gc.s.rdown,
             aes(x = signway,
                 y = interaction_name_2)) + 
  geom_point(aes(size = pval, 
                 color = receptor.logFC)) +
  geom_point(shape = 1,
             aes(size = pval), 
             color="black")+
  scale_size(rang = c(1.5,6)) +
  labs(x=NULL,y=NULL,size="P value",
       title = "Down-regulated signaling of EC-GC in DKD") +
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.text =element_text(size = 8, color = "black"),
        axis.text.y = element_text(face="italic"),
        axis.text.x=element_text(angle = 90,
                                 hjust = 1,
                                 vjust = 0.5)) +
  scale_color_gradient2(low='navy',
                        high="#ADD8E6") +
  geom_hline(yintercept=c(3.5, 4.5, 5.5, 6.5, 7.5, 24.5, 25.5, 26.5,
                          27.5, 31.5, 34.5, 35.5, 39.5, 46.5, 47.5,
                          48.5, 49.5, 51.5, 52.5, 54.5, 55.5, 61.5,
                          62.5, 65.5, 66.5, 67.5, 68.5, 71.5, 72.5, 
                          73.5, 82.5, 84.5, 85.5))+
  theme(plot.margin = margin(t = 5, 
                             r = 5, 
                             b = 5, 
                             l = 5, 
                             unit = "pt"))
p1

labels_and_y <- list("ANGPT" = 2, "ApoA" = 4,  
                     "BAFF" = 5, "CD40" = 6, 
                     "CALCR" = 7, "COLLAGEN" = 16,  
                     "CD46" = 25, "CD99" = 26, 
                     "CDH5" = 27, "EPHA" = 29.5,  
                     "EPHB" = 33, "ESAM" = 35, 
                     "FGF" = 37, "FN1" = 43,  
                     "IGF" = 47, "IL4" = 48, 
                     "IL6" = 49, "IL1" = 50.5,  
                     "IFN-II" = 52, "LIFR" = 53.5, 
                     "LT" = 55, "LAMININ" = 58.5,  
                     "MIF" = 62, "NT" = 64, 
                     "NOTCH" = 66, "Netrin" = 67,  
                     "OSM" = 68, "PDGF" = 70, 
                     "PARs" = 72, "Prostaglandin" = 73,  
                     "TGFb" = 78, "THBS" = 83.5, 
                     "VTN" = 85, "WNT" = 86) 
gene_grodown <- ggplot(gc.s.rdown,
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
  ##得到4.1 EC-GC的source下调受体气泡图A（logFC = 0.2）


gg1 <- netVisual_bubble(cellchat, 
                        pairLR.use = gc.s.rdown.use, 
                        sources.use = "EC-GC", 
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "Down-regulated signaling of EC-GC in DKD")
gg1 <- gg1 + 
  geom_hline(yintercept=c(3.5, 4.5, 5.5, 6.5, 7.5, 24.5, 25.5, 26.5,
                          27.5, 31.5, 34.5, 35.5, 39.5, 46.5, 47.5,
                          48.5, 49.5, 51.5, 52.5, 54.5, 55.5, 61.5,
                          62.5, 65.5, 66.5, 67.5, 68.5, 71.5, 72.5, 
                          73.5, 82.5, 84.5, 85.5))
gg1 %>% insert_left(gene_grodown, width = 0.3)
  ##得到4.2 EC-GC的source下调受体气泡图B（logFC = 0.2）


##EC-GC作为target
gc.t.rdown <- gc.t.rdown %>%  
  arrange(substr(pathway_name, 1, 1)) 

unique(gc.t.rdown$pathway_name)
gc.t.rdown.use = gc.t.rdown[, "interaction_name", drop = F]
gc.t.rdown$interaction_name <- as.factor(gc.t.rdown$interaction_name)
gc.t.rdown$interaction_name <- fct_inorder(gc.t.rdown$interaction_name)
gc.t.rdown$interaction_name_2 <- as.factor(gc.t.rdown$interaction_name_2)
gc.t.rdown$interaction_name_2 <- fct_inorder(gc.t.rdown$interaction_name_2)
gc.t.rdown$pathway_name <- as.factor(gc.t.rdown$pathway_name)
gc.t.rdown$pathway_name <- fct_inorder(gc.t.rdown$pathway_name)
gc.t.rdown <- gc.t.rdown %>% 
  mutate(p="")

gc.t.rdown$signway <- paste(gc.t.rdown$source,
                          "  ->  ", 
                          gc.t.rdown$target) 

table(unique(gc.t.rdown$interaction_name_2))
summary_data <- gc.t.rdown %>%  
  group_by(pathway_name) %>%
  summarise(interaction_count = n_distinct(interaction_name_2))
print(summary_data)

p1 <- ggplot(gc.t.rdown,
             aes(x = signway,
                 y = interaction_name_2)) + 
  geom_point(aes(size = pval, 
                 color = receptor.logFC)) +
  geom_point(shape = 1,
             aes(size = pval), 
             color="black")+
  scale_size(rang = c(1.5,6)) +
  labs(x=NULL,y=NULL,size="P value",
       title = "Down-regulated signaling of EC-GC in DKD") +
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.text =element_text(size = 8, color = "black"),
        axis.text.y = element_text(face="italic"),
        axis.text.x=element_text(angle = 90,
                                 hjust = 1,
                                 vjust = 0.5)) +
  scale_color_gradient2(low='navy',
                        high="#ADD8E6") +
  geom_hline(yintercept=c(3.5, 4.5, 5.5, 21.5, 22.5, 23.5, 24.5,
                          28.5, 31.5, 32.5, 38.5, 42.5, 43.5, 44.5, 
                          45.5, 47.5, 63.5, 66.5, 69.5, 70.5, 71.5,
                          75.5, 77.5, 79.5, 81.5))+
  theme(plot.margin = margin(t = 5, 
                             r = 5, 
                             b = 5, 
                             l = 5, 
                             unit = "pt"))
p1

labels_and_y <- list("ANGPT" = 2, "CD40" = 4,  
                     "CALCR" = 5, "COLLAGEN" = 13.5,
                     "CD99" = 22, "CDH1" = 23,  
                     "CDH5" = 24, "EPHA" = 26.5,
                     "EPHB" = 30, "ESAM" = 32,  
                     "FGF" = 35.5, "FN1" = 40.5,
                     "IL6" = 43, "IL1" = 44,  
                     "JAM" = 45, "LIFR" = 46.5,
                     "LAMININ" = 55.5, "MK" = 65,  
                     "NOTCH" = 68, "OSM" = 70,
                     "PTN" = 71, "SPP1" = 73.5,  
                     "SELL" = 76.5, "TENASCIN" = 78.5,
                     "VTN" = 80.5)  
gene_grodown <- ggplot(gc.t.rdown,
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
##得到4.3 EC-GC的target下调受体气泡图A（logFC = 0.2）

gg1 <- netVisual_bubble(cellchat, 
                        pairLR.use = gc.t.rdown.use, 
                        targets.use = "EC-GC", 
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "Down-regulated signaling of EC-GC in DKD")
gg1 <- gg1 + 
  geom_hline(yintercept=c(3.5, 4.5, 5.5, 21.5, 22.5, 23.5, 24.5,
                          28.5, 31.5, 32.5, 38.5, 42.5, 43.5, 44.5, 
                          45.5, 47.5, 63.5, 66.5, 69.5, 70.5, 71.5,
                          75.5, 77.5, 79.5, 81.5))
gg1 %>% insert_left(gene_grodown, width = 0.3)
  ##得到4.4 EC-GC的target下调受体气泡图B（logFC = 0.2）


