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



####EC-GC的上调配体气泡图（source）#####################
#EC-GC的上调配体气泡图（source）
gc.s.lup <- ligand.up %>%  
  filter(source == "EC-GC")
gc.s.lup <- gc.s.lup %>%
  filter(target %in% c("B", "T", "FIB", "MC", "POD", 
                       "EC-GC", "DC", "MAC"))


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
write.csv(summary_data, "data.csv")

p1 <- ggplot(gc.s.lup,
             aes(x = signway,
                 y = interaction_name_2)) + 
  geom_point(aes(size = pval, 
                 color = ligand.logFC)) +
  geom_point(shape = 1,
             aes(size = pval), 
             color="black")+
  scale_size(rang = c(1.5)) +
  labs(x=NULL,y=NULL,size="P value",
       title = "Up-regulated Ligand of EC-GC in DKD") +
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
  ##得到1.1 EC-GC的source上调配体气泡图A


gg1 <- netVisual_bubble(cellchat, 
                        pairLR.use = gc.s.lup.use, 
                        sources.use = "EC-GC",
                        targets.use = c("B", "T", "FIB", "MC", 
                                        "POD", "EC-GC", "DC", 
                                        "MAC"),
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "Up-regulated Ligand of EC-GC in DKD")
gg1 <- gg1 + 
  geom_hline(yintercept=c(2.5, 5.5, 6.5))
gg1 %>% insert_left(gene_group, width = 0.3)
  ##得到1.2 EC-GC的source上调配体气泡图B



####EC-GC的上调配体气泡图（target）#####################
#EC-GC的上调配体气泡图（target）
gc.t.lup <- ligand.up %>%  
  filter(target == "EC-GC")
gc.t.lup <- gc.t.lup %>%
  filter(source %in% c("B", "T", "FIB", "MC", "POD", 
                       "EC-GC", "DC", "MAC"))
gc.t.lup <- gc.t.lup %>%  
  arrange(substr(pathway_name, 1, 1)) 

unique(gc.t.lup$pathway_name)
gc.t.lup.use = gc.t.lup[, "interaction_name", drop = F]
gc.t.lup$interaction_name <- as.factor(gc.t.lup$interaction_name)
gc.t.lup$interaction_name <- fct_inorder(gc.t.lup$interaction_name)
gc.t.lup$interaction_name_2 <- as.factor(gc.t.lup$interaction_name_2)
gc.t.lup$interaction_name_2 <- fct_inorder(gc.t.lup$interaction_name_2)
gc.t.lup$pathway_name <- as.factor(gc.t.lup$pathway_name)
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
write.csv(summary_data, "data.csv")

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
       title = "Up-regulated Ligand of EC-GC in DKD") +
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
  geom_hline(yintercept=c(1.5, 5.5, 6.5, 9.5, 24.5, 25.5, 26.5, 
                          28.5, 29.5, 30.5, 33.5, 36.5, 42.5, 
                          43.5, 45.5, 46.5, 48.5, 52.5, 56.5,
                          64.5, 65.5, 66.5, 68.5, 69.5, 70.5,
                          73.5, 76.5, 77.5, 78.5, 80.5))+
  theme(plot.margin = margin(t = 5, 
                             r = 5, 
                             b = 5, 
                             l = 5, 
                             unit = "pt"))
p1

labels_and_y <- list("APP" = 1, 
                     "CCL" = 3.5, 
                     "CXCL" = 6, 
                     "COMPLEMENT" = 8, 
                     "COLLAGEN" = 17, 
                     "CD45" = 25, 
                     "CD48" = 26, 
                     "CD86" = 27.5, 
                     "CDH" = 29, 
                     "CDH1" = 30, 
                     "EPHA" = 32, 
                     "EPHB" = 35, 
                     "FN1" = 39.5, 
                     "GAS" = 43, 
                     "IL1" = 44.5, 
                     "LT" = 46, 
                     "MIF" = 47.5, 
                     "MK" = 50.5, 
                     "MHC-I" = 54.5, 
                     "MHC-II" = 60.5, 
                     "NRG" = 65, 
                     "Netrin" = 66, 
                     "PDGF" = 67.5, 
                     "PSAP" = 69, 
                     "PECAM1" = 70, 
                     "SPP1" = 72, 
                     "TGFb" = 75, 
                     "TRAIL" = 77, 
                     "TENASCIN" = 78, 
                     "THY1" = 79.5, 
                     "VISTA" = 81)  
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
  ##得到1.3 EC-GC的target上调配体气泡图A

gg1 <- netVisual_bubble(cellchat, 
                        pairLR.use = gc.t.lup.use, 
                        targets.use = "EC-GC", 
                        sources.use = c("B", "T", "FIB", "MC", 
                                        "POD", "EC-GC", "DC", 
                                        "MAC"),
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "Up-regulated Ligand of EC-GC in DKD")
gg1 <- gg1 + 
  geom_hline(yintercept=c(1.5, 5.5, 6.5, 9.5, 24.5, 25.5, 26.5, 
                          28.5, 29.5, 30.5, 33.5, 36.5, 42.5, 
                          43.5, 45.5, 46.5, 48.5, 52.5, 56.5,
                          64.5, 65.5, 66.5, 68.5, 69.5, 70.5,
                          73.5, 76.5, 77.5, 78.5, 80.5))
gg1 %>% insert_left(gene_group, width = 0.3)
  ##得到1.4 EC-GC的target上调配体气泡图B



####EC-GC的上调受体气泡图（source）#####################
#EC-GC的上调受体气泡图（source）
gc.s.rup <- receptor.up %>%  
  filter(source == "EC-GC")
gc.s.rup <- gc.s.rup %>%
  filter(target %in% c("B", "T", "FIB", "MC", "POD", 
                       "EC-GC", "DC", "MAC"))

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
write.csv(summary_data, "data.csv")

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
       title = "Up-regulated Receptor of EC-GC in DKD") +
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
  geom_hline(yintercept=c(1.5, 4.5, 8.5, 11.5, 12.5, 13.5, 
                          15.5, 16.5, 18.5, 23.5, 25.5, 26.5, 
                          27.5, 28.5, 31.5, 32.5, 33.5, 34.5,
                          38.5, 39.5, 40.5, 47.5, 50.5, 51.5, 
                          53.5, 60.5, 62.5, 63.5, 64.5, 65.5,
                          67.5, 70.5))+
  theme(plot.margin = margin(t = 5, 
                             r = 5, 
                             b = 5, 
                             l = 5, 
                             unit = "pt"))
p1

labels_and_y <- list("APP" = 1, 
                     "ADGRE" = 3, 
                     "ADGRG" = 6.5, 
                     "BMP" = 10, 
                     "BAFF" = 12, 
                     "CX3C" = 13, 
                     "CSF" = 14.5, 
                     "CD40" = 16, 
                     "COMPLEMENT" = 17.5, 
                     "COLLAGEN" = 21, 
                     "CD23" = 24.5, 
                     "CD99" = 26, 
                     "CDH" = 27, 
                     "CNTN" = 28, 
                     "FN1" = 30, 
                     "GAS" = 32, 
                     "GP1BA" = 33, 
                     "IGF" = 34, 
                     "IL2" = 36.5, 
                     "IL4" = 39, 
                     "IGFBP" = 40, 
                     "ICAM" = 44, 
                     "LAMININ" = 49, 
                     "L1CAM" = 51, 
                     "MIF" = 52.5, 
                     "MHC-II" = 57, 
                     "NOTCH" = 61.5, 
                     "PECAM1" = 63, 
                     "SIRP" = 64, 
                     "THBS" = 65, 
                     "THY1" = 66.5, 
                     "VEGF" = 69, 
                     "VCAM" = 71.5) 
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
  ##得到2.1 EC-GC的source上调受体气泡图A


gg1 <- netVisual_bubble(cellchat, 
                        pairLR.use = gc.s.rup.use, 
                        sources.use = "EC-GC",
                        targets.use = c("B", "T", "FIB", "MC",
                                        "POD", "EC-GC", "DC",
                                        "MAC"),
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "Up-regulated Receptor of EC-GC in DKD")
gg1 <- gg1 + 
  geom_hline(yintercept=c(1.5, 4.5, 8.5, 11.5, 12.5, 13.5, 
                          15.5, 16.5, 18.5, 23.5, 25.5, 26.5, 
                          27.5, 28.5, 31.5, 32.5, 33.5, 34.5,
                          38.5, 39.5, 40.5, 47.5, 50.5, 51.5, 
                          53.5, 60.5, 62.5, 63.5, 64.5, 65.5,
                          67.5, 70.5))
gg1 %>% insert_left(gene_group, width = 0.3)
  ##得到2.2 EC-GC的source上调受体气泡图B



####EC-GC的上调受体气泡图（target）#####################
#EC-GC的上调受体气泡图（target）
gc.t.rup <- receptor.up %>%  
  filter(target == "EC-GC")
gc.t.rup <- gc.t.rup %>%
  filter(source %in% c("B", "T", "FIB", "MC", "POD", 
                       "EC-GC", "DC", "MAC"))

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
write.csv(summary_data, "data.csv")

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
       title = "Up-regulated Receptor of EC-GC in DKD") +
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
  geom_hline(yintercept=c(1.5, 2.5, 3.5))+
  theme(plot.margin = margin(t = 5, 
                             r = 5, 
                             b = 5, 
                             l = 5, 
                             unit = "pt"))
p1

labels_and_y <- list("CNTN" = 1, 
                     "NOTCH" = 2, 
                     "PECAM1" = 3, 
                     "VEGF" = 5)  
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
  ##得到2.3 EC-GC的target上调受体气泡图A

gg1 <- netVisual_bubble(cellchat, 
                        pairLR.use = gc.t.rup.use, 
                        targets.use = "EC-GC", 
                        sources.use = c("B", "T", "FIB", "MC",
                                        "POD", "EC-GC", "DC",
                                        "MAC"),
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "Up-regulated Receptor of EC-GC in DKD")
gg1 <- gg1 + 
  geom_hline(yintercept=c(1.5, 2.5, 3.5))
gg1 %>% insert_left(gene_group, width = 0.3)
  ##得到2.4 EC-GC的target上调受体气泡图B



####EC-GC的下调配体气泡图（source）#####################
#EC-GC的下调配体气泡图（source）
gc.s.ldown <- ligand.down %>%  
  filter(source == "EC-GC")
gc.s.ldown <- gc.s.ldown %>%
  filter(target %in% c("B", "T", "FIB", "MC", "POD", 
                       "EC-GC", "DC", "MAC"))

gc.s.ldown <- gc.s.ldown %>%  
  arrange(substr(pathway_name, 1, 1)) 
rows_to_move <- gc.s.ldown[11:12, ]  
remaining_rows <- gc.s.ldown[-(11:12), ]  
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
write.csv(summary_data, "data.csv")

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
       title = "Down-regulated Ligand of EC-GC in DKD") +
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
  geom_hline(yintercept=c(2.5, 3.5, 6.5, 7.5, 8.5, 20.5, 22.5, 
                          23.5, 24.5, 25.5, 26.5, 28.5, 36.5,
                          37.5, 41.5, 42.5, 43.5, 44.5))+
  theme(plot.margin = margin(t = 5, 
                             r = 5, 
                             b = 5, 
                             l = 5, 
                             unit = "pt"))
p1

labels_and_y <- list("APP" = 1.5, 
                     "ADGRE" = 3, 
                     "ADGRG" = 5, 
                     "BTLA" = 7, 
                     "CX3C" = 8, 
                     "COLLAGEN" = 14.5, 
                     "CD99" = 21.5, 
                     "CDH5" = 23, 
                     "ESAM" = 24, 
                     "GAS" = 25, 
                     "GRN" = 26, 
                     "ICAM" = 27.5, 
                     "LAMININ" = 32.5, 
                     "MHC-I" = 37, 
                     "MHC-II" = 39.5, 
                     "Netrin" = 42, 
                     "PDGF" = 43, 
                     "SEMA6" = 44, 
                     "VCAM" = 45) 
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
  ##得到3.1 EC-GC的source下调配体气泡图A


gg1 <- netVisual_bubble(cellchat, 
                        pairLR.use = gc.s.ldown.use, 
                        sources.use = "EC-GC",
                        targets.use = c("B", "T", "FIB", "MC", "POD", 
                                        "EC-GC", "DC", "MAC"),
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "Down-regulated Ligand of EC-GC in DKD")
gg1 <- gg1 + 
  geom_hline(yintercept=c(2.5, 3.5, 6.5, 7.5, 8.5, 20.5, 22.5, 
                          23.5, 24.5, 25.5, 26.5, 28.5, 36.5,
                          37.5, 41.5, 42.5, 43.5, 44.5))
gg1 %>% insert_left(gene_grodown, width = 0.3)
  ##得到3.2 EC-GC的source下调配体气泡图B



####EC-GC的下调配体气泡图（target）#####################
#EC-GC的下调配体气泡图（target）
gc.t.ldown <- ligand.down %>%  
  filter(target == "EC-GC")
gc.t.ldown <- gc.t.ldown %>%
  filter(source %in% c("B", "T", "FIB", "MC", "POD", 
                       "EC-GC", "DC", "MAC"))

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
write.csv(summary_data, "data.csv")

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
       title = "Down-regulated Ligand of EC-GC in DKD") +
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
  geom_hline(yintercept=c(2.5, 3.5, 5.5, 6.5, 12.5, 13.5, 
                          16.5, 25.5, 26.5, 27.5, 29.5, 30.5,
                          31.5, 33.5, 34.5, 38.5, 39.5, 41.5,
                          43.5, 44.5, 56.5, 57.5, 62.5, 69.5,
                          73.5, 76.5, 77.5, 78.5, 79.5, 80.5,
                          85.5))+
  theme(plot.margin = margin(t = 5, 
                             r = 5, 
                             b = 5, 
                             l = 5, 
                             unit = "pt"))
p1

labels_and_y <- list("ANGPT" = 1.5, 
                     "APP" = 3, 
                     "BMP" = 4.5, 
                     "BTLA" = 6, 
                     "CCL" = 9.5, 
                     "CX3C" = 13, 
                     "COMPLEMENT" = 15, 
                     "COLLAGEN" = 21, 
                     "CD46" = 26, 
                     "CD48" = 27, 
                     "CD86" = 28.5, 
                     "CD99" = 30, 
                     "CDH5" = 31, 
                     "DHEAS" = 32.5, 
                     "ESAM" = 34, 
                     "FGF" = 36.5, 
                     "GAS" = 39, 
                     "IL1" = 40.5, 
                     "ICAM" = 42.5, 
                     "KLK" = 44, 
                     "LAMININ" = 50.5, 
                     "MIF" = 57, 
                     "MHC-I" = 60, 
                     "MHC-II" = 66, 
                     "NT" = 71.5, 
                     "NOTCH" = 75, 
                     "Netrin" = 77, 
                     "PDGF" = 78, 
                     "PLAU" = 79, 
                     "SEMA6" = 80, 
                     "VEGF" = 83, 
                     "VISFATIN" = 86)  
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
  ##得到3.3 EC-GC的target下调配体气泡图A

gg1 <- netVisual_bubble(cellchat, 
                        pairLR.use = gc.t.ldown.use, 
                        targets.use = "EC-GC", 
                        sources.use = c("B", "T", "FIB", "MC", "POD", 
                                        "EC-GC", "DC", "MAC"),
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "Down-regulated Ligand of EC-GC in DKD")
gg1 <- gg1 + 
  geom_hline(yintercept=c(2.5, 3.5, 5.5, 6.5, 12.5, 13.5, 
                          16.5, 25.5, 26.5, 27.5, 29.5, 30.5,
                          31.5, 33.5, 34.5, 38.5, 39.5, 41.5,
                          43.5, 44.5, 56.5, 57.5, 62.5, 69.5,
                          73.5, 76.5, 77.5, 78.5, 79.5, 80.5,
                          85.5))
gg1 %>% insert_left(gene_grodown, width = 0.3)
  ##得到3.4 EC-GC的target下调配体气泡图B



####EC-GC的下调受体气泡图（source）#####################
#EC-GC的下调受体气泡图（source）
gc.s.rdown <- receptor.down %>%  
  filter(source == "EC-GC")
gc.s.rdown <- gc.s.rdown %>% 
  filter(target %in% c("B", "T", "FIB", "MC", "POD", 
                       "EC-GC", "DC", "MAC"))

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
write.csv(summary_data, "data.csv")

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
       title = "Down-regulated Receptor of EC-GC in DKD") +
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
  geom_hline(yintercept=c(2.5, 3.5, 5.5, 6.5, 12.5, 13.5, 
                          16.5, 25.5, 26.5, 27.5, 29.5, 30.5, 
                          31.5, 33.5, 34.5, 38.5, 39.5, 41.5, 
                          43.5, 44.5, 56.5, 57.5, 62.5, 69.5, 
                          73.5, 76.5, 77.5, 78.5, 79.5, 80.5, 
                          85.5))+
  theme(plot.margin = margin(t = 5, 
                             r = 5, 
                             b = 5, 
                             l = 5, 
                             unit = "pt"))
p1

labels_and_y <- list("ANGPT" = 2, 
                     "APP" = 4.5, 
                     "BAFF" = 6, 
                     "CXCL" = 7, 
                     "CD40" = 8, 
                     "COMPLEMENT" = 9.5, 
                     "CALCR" = 11, 
                     "COLLAGEN" = 18.5, 
                     "CD46" = 26, 
                     "CD99" = 27, 
                     "CDH5" = 28, 
                     "EPHA" = 30, 
                     "EPHB" = 34.5, 
                     "ESAM" = 38, 
                     "FGF" = 41, 
                     "FN1" = 46.5, 
                     "GH" = 50, 
                     "IGF" = 51, 
                     "IL4" = 52, 
                     "IL6" = 53, 
                     "IL1" = 54.5, 
                     "IFN-II" = 56, 
                     "IGFBP" = 57, 
                     "LIFR" = 58.5, 
                     "LAMININ" = 62.5, 
                     "MSTN" = 66, 
                     "MIF" = 67.5, 
                     "NT" = 70, 
                     "NOTCH" = 72, 
                     "Netrin" = 73, 
                     "OSM" = 74, 
                     "PDGF" = 76, 
                     "PARs" = 78, 
                     "PLAU" = 79, 
                     "TGFb" = 84, 
                     "THBS" = 89, 
                     "VTN" = 90.5) 
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
  ##得到4.1 EC-GC的source下调受体气泡图A


gg1 <- netVisual_bubble(cellchat, 
                        pairLR.use = gc.s.rdown.use, 
                        sources.use = "EC-GC",
                        targets.use = c("B", "T", "FIB", "MC",
                                        "POD", "EC-GC", "DC",
                                        "MAC"),
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "Down-regulated Receptor of EC-GC in DKD")
gg1 <- gg1 + 
  geom_hline(yintercept=c(2.5, 3.5, 5.5, 6.5, 12.5, 13.5, 
                          16.5, 25.5, 26.5, 27.5, 29.5, 30.5, 
                          31.5, 33.5, 34.5, 38.5, 39.5, 41.5, 
                          43.5, 44.5, 56.5, 57.5, 62.5, 69.5, 
                          73.5, 76.5, 77.5, 78.5, 79.5, 80.5, 
                          85.5))
gg1 %>% insert_left(gene_grodown, width = 0.3)
  ##得到4.2 EC-GC的source下调受体气泡图B


####EC-GC的下调受体气泡图（target）#####################
#EC-GC的下调受体气泡图（target）
gc.t.rdown <- receptor.down %>%  
  filter(target == "EC-GC")
gc.t.rdown <- gc.t.rdown %>% 
  filter(source %in% c("B", "T", "FIB", "MC", "POD", 
                       "EC-GC", "DC", "MAC"))

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
write.csv(summary_data, "data.csv")

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
       title = "Down-regulated Receptor of EC-GC in DKD") +
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
  geom_hline(yintercept=c(3.5, 4.5, 5.5, 17.5, 18.5, 19.5, 
                          20.5, 24.5, 27.5, 28.5, 33.5, 37.5,
                          38.5, 39.5, 41.5, 49.5, 51.5, 52.5,
                          53.5, 54.5, 55.5, 58.5, 60.5, 61.5,
                          62.5))+
  theme(plot.margin = margin(t = 5, 
                             r = 5, 
                             b = 5, 
                             l = 5, 
                             unit = "pt"))
p1

labels_and_y <- list("ANGPT" = 2, 
                     "CD40" = 4, 
                     "CALCR" = 5, 
                     "COLLAGEN" = 11.5, 
                     "CD99" = 18, 
                     "CDH1" = 19, 
                     "CDH5" = 20, 
                     "EPHA" = 22.5, 
                     "EPHB" = 26, 
                     "ESAM" = 28, 
                     "FGF" = 31, 
                     "FN1" = 35.5, 
                     "IL6" = 38, 
                     "IL1" = 39, 
                     "LIFR" = 40.5, 
                     "LAMININ" = 45.5, 
                     "MK" = 50.5, 
                     "NOTCH" = 52, 
                     "OSM" = 53, 
                     "PTN" = 54, 
                     "Prostaglandin" = 55, 
                     "SPP1" = 57, 
                     "SELL" = 59.5, 
                     "TENASCIN" = 61, 
                     "VISFATIN" = 62, 
                     "VTN" = 63.5)  
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
  ##得到4.3 EC-GC的target下调受体气泡图A

gg1 <- netVisual_bubble(cellchat, 
                        pairLR.use = gc.t.rdown.use, 
                        targets.use = "EC-GC", 
                        sources.use = c("B", "T", "FIB", "MC",
                                        "POD", "EC-GC", "DC",
                                        "MAC"),
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "Down-regulated Receptor of EC-GC in DKD")
gg1 <- gg1 + 
  geom_hline(yintercept=c(3.5, 4.5, 5.5, 17.5, 18.5, 19.5, 
                          20.5, 24.5, 27.5, 28.5, 33.5, 37.5,
                          38.5, 39.5, 41.5, 49.5, 51.5, 52.5,
                          53.5, 54.5, 55.5, 58.5, 60.5, 61.5,
                          62.5))
gg1 %>% insert_left(gene_grodown, width = 0.3)
  ##得到4.4 EC-GC的target下调受体气泡图B


