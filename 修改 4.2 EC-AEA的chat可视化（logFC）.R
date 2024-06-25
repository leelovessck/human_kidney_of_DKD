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



####EC-AEA的上调配体气泡图（source）#####################
#EC-AEA的上调配体气泡图（source）
aea.s.lup <- ligand.up %>%  
  filter(source == "EC-AEA")
aea.s.lup <- aea.s.lup %>%
  filter(target %in% c("B", "T", "FIB", "MC", "POD", 
                       "EC-AEA", "DC", "MAC", "VSMC"))


unique(aea.s.lup$pathway_name)
aea.s.lup.use = aea.s.lup[, "interaction_name", drop = F]

aea.s.lup$interaction_name <- as.factor(aea.s.lup$interaction_name)
aea.s.lup$interaction_name <- fct_inorder(aea.s.lup$interaction_name)
aea.s.lup$interaction_name_2 <- as.factor(aea.s.lup$interaction_name_2)
aea.s.lup$interaction_name_2 <- fct_inorder(aea.s.lup$interaction_name_2)
aea.s.lup$pathway_name <- as.factor(aea.s.lup$pathway_name)
aea.s.lup$pathway_name <- fct_inorder(aea.s.lup$pathway_name)
aea.s.lup <- aea.s.lup %>% 
  mutate(p="")
aea.s.lup$signway <- paste(aea.s.lup$source,
                          "  ->  ", 
                          aea.s.lup$target) 

table(unique(aea.s.lup$interaction_name_2))
summary_data <- aea.s.lup %>%  
  group_by(pathway_name) %>%
  summarise(interaction_count = n_distinct(interaction_name_2))
print(summary_data)
write.csv(summary_data, "data.csv")

p1 <- ggplot(aea.s.lup,
             aes(x = signway,
                 y = interaction_name_2)) + 
  geom_point(aes(size = pval, 
                 color = ligand.logFC)) +
  geom_point(shape = 1,
             aes(size = pval), 
             color="black")+
  scale_size(rang = c(1.5)) +
  labs(x=NULL,y=NULL,size="P value",
       title = "Up-regulated Ligand of EC-AEA in DKD") +
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
  geom_hline(yintercept=c(3.5, 4.5, 5.5, 11.5, 12.5))+
  theme(plot.margin = margin(t = 5, 
                             r = 5, 
                             b = 5, 
                             l = 5, 
                             unit = "pt"))
p1

labels_and_y <- list("MIF" = 2, 
                     "SEMA3" = 4, 
                     "PSAP" = 5, 
                     "EPHA" = 8.5, 
                     "MHC-II" = 12, 
                     "PECAM1" = 13) 
gene_group <- ggplot(aea.s.lup,
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
  ##得到1.1 EC-AEA的source上调配体气泡图A


gg1 <- netVisual_bubble(cellchat, 
                        pairLR.use = aea.s.lup.use, 
                        sources.use = "EC-AEA",
                        targets.use = c("B", "T", "FIB", "MC", 
                                        "POD", "EC-AEA", "DC", 
                                        "MAC", "VSMC"),
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "Up-regulated Ligand of EC-AEA in DKD")
gg1 <- gg1 + 
  geom_hline(yintercept=c(3.5, 4.5, 5.5, 11.5, 12.5))
gg1 %>% insert_left(gene_group, width = 0.3)
  ##得到1.2 EC-AEA的source上调配体气泡图B



####EC-AEA的上调配体气泡图（target）#####################
#EC-AEA的上调配体气泡图（target）
aea.t.lup <- ligand.up %>%  
  filter(target == "EC-AEA")
aea.t.lup <- aea.t.lup %>%
  filter(source %in% c("B", "T", "FIB", "MC", "POD", 
                       "EC-AEA", "DC", "MAC", "VSMC"))
aea.t.lup <- aea.t.lup %>%  
  arrange(substr(pathway_name, 1, 1)) 

unique(aea.t.lup$pathway_name)
aea.t.lup.use = aea.t.lup[, "interaction_name", drop = F]
aea.t.lup$interaction_name <- as.factor(aea.t.lup$interaction_name)
aea.t.lup$interaction_name <- fct_inorder(aea.t.lup$interaction_name)
aea.t.lup$interaction_name_2 <- as.factor(aea.t.lup$interaction_name_2)
aea.t.lup$interaction_name_2 <- fct_inorder(aea.t.lup$interaction_name_2)
aea.t.lup$pathway_name <- as.factor(aea.t.lup$pathway_name)
aea.t.lup$pathway_name <- fct_inorder(aea.t.lup$pathway_name)
aea.t.lup <- aea.t.lup %>% 
  mutate(p="")

aea.t.lup$signway <- paste(aea.t.lup$source,
                          "  ->  ", 
                          aea.t.lup$target) 

table(unique(aea.t.lup$interaction_name_2))
summary_data <- aea.t.lup %>%  
  group_by(pathway_name) %>%
  summarise(interaction_count = n_distinct(interaction_name_2))
print(summary_data)
write.csv(summary_data, "data.csv")

p1 <- ggplot(aea.t.lup,
             aes(x = signway,
                 y = interaction_name_2)) + 
  geom_point(aes(size = pval, 
                 color = ligand.logFC)) +
  geom_point(shape = 1,
             aes(size = pval), 
             color="black")+
  scale_size(rang = c(1.5,6)) +
  labs(x=NULL,y=NULL,size="P value",
       title = "Up-regulated Ligand of EC-AEA in DKD") +
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
  geom_hline(yintercept=c(1.5, 5.5, 6.5, 9.5, 27.5, 28.5, 29.5, 
                          31.5, 32.5, 33.5, 45.5, 51.5, 52.5, 54.5, 
                          55.5, 58.5, 63.5, 67.5, 75.5, 76.5, 77.5,
                          80.5, 82.5, 83.5, 84.5, 89.5, 90.5, 93.5,
                          94.5, 95.5, 97.5, 99.5, 101.5))+
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
                     "COLLAGEN" = 18.5, 
                     "CD45" = 28, 
                     "CD48" = 29, 
                     "CD86" = 30.5, 
                     "CDH" = 32, 
                     "CDH1" = 33, 
                     "EPHA" = 39.5, 
                     "FN1" = 48.5, 
                     "GAS" = 52, 
                     "IL1" = 53.5, 
                     "LT" = 55, 
                     "MIF" = 57, 
                     "MK" = 61, 
                     "MHC-I" = 65.5, 
                     "MHC-II" = 71.5, 
                     "NRG" = 76, 
                     "Netrin" = 77, 
                     "PDGF" = 79, 
                     "PTN" = 81.5, 
                     "PSAP" = 83, 
                     "PECAM1" = 84, 
                     "SPP1" = 87, 
                     "SEMA3" = 90, 
                     "TGFb" = 92, 
                     "TRAIL" = 94, 
                     "TENASCIN" = 95, 
                     "THBS" = 96.5, 
                     "THY1" = 98.5, 
                     "VCAM" = 100.5, 
                     "VISTA" = 102)  
gene_group <- ggplot(aea.t.lup,
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
  ##得到1.3 EC-AEA的target上调配体气泡图A

gg1 <- netVisual_bubble(cellchat, 
                        pairLR.use = aea.t.lup.use, 
                        targets.use = "EC-AEA", 
                        sources.use = c("B", "T", "FIB", "MC", 
                                        "POD", "EC-AEA", "DC", 
                                        "MAC", "VSMC"),
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "Up-regulated Ligand of EC-AEA in DKD")
gg1 <- gg1 + 
  geom_hline(yintercept=c(1.5, 5.5, 6.5, 9.5, 27.5, 28.5, 29.5, 
                          31.5, 32.5, 33.5, 45.5, 51.5, 52.5, 54.5, 
                          55.5, 58.5, 63.5, 67.5, 75.5, 76.5, 77.5,
                          80.5, 82.5, 83.5, 84.5, 89.5, 90.5, 93.5,
                          94.5, 95.5, 97.5, 99.5, 101.5))
gg1 %>% insert_left(gene_group, width = 0.3)
  ##得到1.4 EC-AEA的target上调配体气泡图B



####EC-AEA的上调受体气泡图（source）#####################
#EC-AEA的上调受体气泡图（source）
aea.s.rup <- receptor.up %>%  
  filter(source == "EC-AEA")
aea.s.rup <- aea.s.rup %>%
  filter(target %in% c("B", "T", "FIB", "MC", "POD", 
                       "EC-AEA", "DC", "MAC", "VSMC"))

aea.s.rup <- aea.s.rup %>%  
  arrange(substr(pathway_name, 1, 1)) 

unique(aea.s.rup$pathway_name)
aea.s.rup.use = aea.s.rup[, "interaction_name", drop = F]
aea.s.rup$interaction_name <- as.factor(aea.s.rup$interaction_name)
aea.s.rup$interaction_name <- fct_inorder(aea.s.rup$interaction_name)
aea.s.rup$interaction_name_2 <- as.factor(aea.s.rup$interaction_name_2)
aea.s.rup$interaction_name_2 <- fct_inorder(aea.s.rup$interaction_name_2)
aea.s.rup$pathway_name <- as.factor(aea.s.rup$pathway_name)
aea.s.rup$pathway_name <- fct_inorder(aea.s.rup$pathway_name)
aea.s.rup <- aea.s.rup %>% 
  mutate(p="")
aea.s.rup$signway <- paste(aea.s.rup$source,
                          "  ->  ", 
                          aea.s.rup$target) 

table(unique(aea.s.rup$interaction_name_2))
summary_data <- aea.s.rup %>%  
  group_by(pathway_name) %>%
  summarise(interaction_count = n_distinct(interaction_name_2))
print(summary_data)
write.csv(summary_data, "data.csv")

p1 <- ggplot(aea.s.rup,
             aes(x = signway,
                 y = interaction_name_2)) + 
  geom_point(aes(size = pval, 
                 color = receptor.logFC)) +
  geom_point(shape = 1,
             aes(size = pval), 
             color="black")+
  scale_size(rang = c(1.5,6)) +
  labs(x=NULL,y=NULL,size="P value",
       title = "Up-regulated Receptor of EC-AEA in DKD") +
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
  geom_hline(yintercept=c(1.5, 4.5, 8.5, 11.5, 12.5, 14.5, 15.5, 17.5, 
                          21.5, 23.5, 24.5, 25.5, 26.5, 30.5, 31.5, 32.5, 
                          33.5, 35.5, 39.5, 40.5, 41.5, 48.5, 49.5, 50.5, 
                          59.5, 60.5, 63.5, 67.5, 74.5, 75.5, 78.5, 79.5, 
                          80.5, 81.5, 82.5, 84.5, 88.5, 90.5))+
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
                     "CSF" = 13.5, 
                     "CD40" = 15, 
                     "COMPLEMENT" = 16.5, 
                     "COLLAGEN" = 19.5, 
                     "CD23" = 22.5, 
                     "CD99" = 24, 
                     "CDH" = 25, 
                     "EDN" = 26, 
                     "FN1" = 28.5, 
                     "GAS" = 31, 
                     "GP1BA" = 32, 
                     "HGF" = 33, 
                     "IGF" = 34.5, 
                     "IL2" = 37.5, 
                     "IL4" = 40, 
                     "IGFBP" = 41, 
                     "ICAM" = 45, 
                     "JAM" = 49, 
                     "LT" = 50, 
                     "LAMININ" = 55, 
                     "L1CAM" = 60, 
                     "MIF" = 62, 
                     "MK" = 65.5, 
                     "MHC-II" = 71, 
                     "NT" = 75, 
                     "NOTCH" = 77, 
                     "PECAM1" = 79, 
                     "SIRP" = 80, 
                     "TENASCIN" = 81, 
                     "THBS" = 82, 
                     "THY1" = 83.5, 
                     "VEGF" = 86.5, 
                     "VCAM" = 89.5, 
                     "WNT" = 91.5) 
gene_group <- ggplot(aea.s.rup,
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
  ##得到2.1 EC-AEA的source上调受体气泡图A


gg1 <- netVisual_bubble(cellchat, 
                        pairLR.use = aea.s.rup.use, 
                        sources.use = "EC-AEA",
                        targets.use = c("B", "T", "FIB", "MC",
                                        "POD", "EC-AEA", "DC",
                                        "MAC", "VSMC"),
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "Up-regulated Receptor of EC-AEA in DKD")
gg1 <- gg1 + 
  geom_hline(yintercept=c(1.5, 4.5, 8.5, 11.5, 12.5, 14.5, 15.5, 17.5, 
                          21.5, 23.5, 24.5, 25.5, 26.5, 30.5, 31.5, 32.5, 
                          33.5, 35.5, 39.5, 40.5, 41.5, 48.5, 49.5, 50.5, 
                          59.5, 60.5, 63.5, 67.5, 74.5, 75.5, 78.5, 79.5, 
                          80.5, 81.5, 82.5, 84.5, 88.5, 90.5))
gg1 %>% insert_left(gene_group, width = 0.3)
  ##得到2.2 EC-AEA的source上调受体气泡图B



####EC-AEA的上调受体气泡图（target）#####################
#EC-AEA的上调受体气泡图（target）
aea.t.rup <- receptor.up %>%  
  filter(target == "EC-AEA")
aea.t.rup <- aea.t.rup %>%
  filter(source %in% c("B", "T", "FIB", "MC", "POD", 
                       "EC-AEA", "DC", "MAC", "VSMC"))

aea.t.rup <- aea.t.rup %>%  
  arrange(substr(pathway_name, 1, 1)) 

unique(aea.t.rup$pathway_name)
aea.t.rup.use = aea.t.rup[, "interaction_name", drop = F]
aea.t.rup$interaction_name <- as.factor(aea.t.rup$interaction_name)
aea.t.rup$interaction_name <- fct_inorder(aea.t.rup$interaction_name)
aea.t.rup$interaction_name_2 <- as.factor(aea.t.rup$interaction_name_2)
aea.t.rup$interaction_name_2 <- fct_inorder(aea.t.rup$interaction_name_2)
aea.t.rup$pathway_name <- as.factor(aea.t.rup$pathway_name)
aea.t.rup$pathway_name <- fct_inorder(aea.t.rup$pathway_name)
aea.t.rup <- aea.t.rup %>% 
  mutate(p="")

aea.t.rup$signway <- paste(aea.t.rup$source,
                          "  ->  ", 
                          aea.t.rup$target) 

table(unique(aea.t.rup$interaction_name_2))
summary_data <- aea.t.rup %>%  
  group_by(pathway_name) %>%
  summarise(interaction_count = n_distinct(interaction_name_2))
print(summary_data)
write.csv(summary_data, "data.csv")

p1 <- ggplot(aea.t.rup,
             aes(x = signway,
                 y = interaction_name_2)) + 
  geom_point(aes(size = pval, 
                 color = receptor.logFC)) +
  geom_point(shape = 1,
             aes(size = pval), 
             color="black")+
  scale_size(rang = c(1.5,6)) +
  labs(x=NULL,y=NULL,size="P value",
       title = "Up-regulated Receptor of EC-AEA in DKD") +
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
  geom_hline(yintercept=c(1.5, 2.5, 6.5))+
  theme(plot.margin = margin(t = 5, 
                             r = 5, 
                             b = 5, 
                             l = 5, 
                             unit = "pt"))
p1

labels_and_y <- list("PECAM1" = 1, 
                     "THBS" = 2, 
                     "VEGF" = 4.5, 
                     "VCAM" = 7.5)  
gene_group <- ggplot(aea.t.rup,
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
  ##得到2.3 EC-AEA的target上调受体气泡图A

gg1 <- netVisual_bubble(cellchat, 
                        pairLR.use = aea.t.rup.use, 
                        targets.use = "EC-AEA", 
                        sources.use = c("B", "T", "FIB", "MC",
                                        "POD", "EC-AEA", "DC",
                                        "MAC", "VSMC"),
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "Up-regulated Receptor of EC-AEA in DKD")
gg1 <- gg1 + 
  geom_hline(yintercept=c(1.5, 2.5, 6.5))
gg1 %>% insert_left(gene_group, width = 0.3)
  ##得到2.4 EC-AEA的target上调受体气泡图B



####EC-AEA的下调配体气泡图（source）#####################
#EC-AEA的下调配体气泡图（source）
aea.s.ldown <- ligand.down %>%  
  filter(source == "EC-AEA")
aea.s.ldown <- aea.s.ldown %>%
  filter(target %in% c("B", "T", "FIB", "MC", "POD", 
                       "EC-AEA", "DC", "MAC", "VSMC"))

aea.s.ldown <- aea.s.ldown %>%  
  arrange(substr(pathway_name, 1, 1)) 
#rows_to_move <- aea.s.ldown[11:12, ]  
#remaining_rows <- aea.s.ldown[-(11:12), ]  
#aea.s.ldown <- rbind(rows_to_move, remaining_rows)  
 
unique(aea.s.ldown$pathway_name)
aea.s.ldown.use = aea.s.ldown[, "interaction_name", drop = F]
aea.s.ldown$interaction_name <- as.factor(aea.s.ldown$interaction_name)
aea.s.ldown$interaction_name <- fct_inorder(aea.s.ldown$interaction_name)
aea.s.ldown$interaction_name_2 <- as.factor(aea.s.ldown$interaction_name_2)
aea.s.ldown$interaction_name_2 <- fct_inorder(aea.s.ldown$interaction_name_2)
aea.s.ldown$pathway_name <- as.factor(aea.s.ldown$pathway_name)
aea.s.ldown$pathway_name <- fct_inorder(aea.s.ldown$pathway_name)
aea.s.ldown <- aea.s.ldown %>% 
  mutate(p="")
aea.s.ldown$signway <- paste(aea.s.ldown$source,
                          "  ->  ", 
                          aea.s.ldown$target) 

table(unique(aea.s.ldown$interaction_name_2))
summary_data <- aea.s.ldown %>%  
  group_by(pathway_name) %>%
  summarise(interaction_count = n_distinct(interaction_name_2))
print(summary_data)
write.csv(summary_data, "data.csv")

p1 <- ggplot(aea.s.ldown,
             aes(x = signway,
                 y = interaction_name_2)) + 
  geom_point(aes(size = pval, 
                 color = ligand.logFC)) +
  geom_point(shape = 1,
             aes(size = pval), 
             color="black")+
  scale_size(rang = c(1.5,6)) +
  labs(x=NULL,y=NULL,size="P value",
       title = "Down-regulated Ligand of EC-AEA in DKD") +
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
  geom_hline(yintercept=c(1.5, 3.5, 4.5, 5.5, 8.5, 10.5, 12.5, 
                          14.5, 18.5))+
  theme(plot.margin = margin(t = 5, 
                             r = 5, 
                             b = 5, 
                             l = 5, 
                             unit = "pt"))
p1

labels_and_y <- list("CD46" = 1, 
                     "CD99" = 2.5, 
                     "ESAM" = 4, 
                     "GAS" = 5, 
                     "IGF" = 7, 
                     "ICAM" = 9.5, 
                     "MHC-I" = 11.5, 
                     "MHC-II" = 13.5, 
                     "NOTCH" = 16.5, 
                     "VEGF" = 20) 
gene_grodown <- ggplot(aea.s.ldown,
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
  ##得到3.1 EC-AEA的source下调配体气泡图A

gg1 <- netVisual_bubble(cellchat, 
                        pairLR.use = aea.s.ldown.use, 
                        sources.use = "EC-AEA",
                        targets.use = c("B", "T", "FIB", "MC", "POD", 
                                        "EC-AEA", "DC", "MAC", "VSMC"),
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "Down-regulated Ligand of EC-AEA in DKD")
gg1 <- gg1 + 
  geom_hline(yintercept=c(1.5, 3.5, 4.5, 5.5, 8.5, 10.5, 12.5, 
                          14.5, 18.5))
gg1 %>% insert_left(gene_grodown, width = 0.3)
  ##得到3.2 EC-AEA的source下调配体气泡图B



####EC-AEA的下调配体气泡图（target）#####################
#EC-AEA的下调配体气泡图（target）
aea.t.ldown <- ligand.down %>%  
  filter(target == "EC-AEA")
aea.t.ldown <- aea.t.ldown %>%
  filter(source %in% c("B", "T", "FIB", "MC", "POD", 
                       "EC-AEA", "DC", "MAC", "VSMC"))

aea.t.ldown <- aea.t.ldown %>%  
  arrange(substr(pathway_name, 1, 1)) 

unique(aea.t.ldown$pathway_name)
aea.t.ldown.use = aea.t.ldown[, "interaction_name", drop = F]
aea.t.ldown$interaction_name <- as.factor(aea.t.ldown$interaction_name)
aea.t.ldown$interaction_name <- fct_inorder(aea.t.ldown$interaction_name)
aea.t.ldown$interaction_name_2 <- as.factor(aea.t.ldown$interaction_name_2)
aea.t.ldown$interaction_name_2 <- fct_inorder(aea.t.ldown$interaction_name_2)
aea.t.ldown$pathway_name <- as.factor(aea.t.ldown$pathway_name)
aea.t.ldown$pathway_name <- fct_inorder(aea.t.ldown$pathway_name)
aea.t.ldown <- aea.t.ldown %>% 
  mutate(p="")

aea.t.ldown$signway <- paste(aea.t.ldown$source,
                          "  ->  ", 
                          aea.t.ldown$target) 

table(unique(aea.t.ldown$interaction_name_2))
summary_data <- aea.t.ldown %>%  
  group_by(pathway_name) %>%
  summarise(interaction_count = n_distinct(interaction_name_2))
print(summary_data)
write.csv(summary_data, "data.csv")

p1 <- ggplot(aea.t.ldown,
             aes(x = signway,
                 y = interaction_name_2)) + 
  geom_point(aes(size = pval, 
                 color = ligand.logFC)) +
  geom_point(shape = 1,
             aes(size = pval), 
             color="black")+
  scale_size(rang = c(1.5,6)) +
  labs(x=NULL,y=NULL,size="P value",
       title = "Down-regulated Ligand of EC-AEA in DKD") +
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
  geom_hline(yintercept=c(2.5, 3.5, 6.5, 12.5, 15.5, 18.5, 19.5, 
                          20.5, 22.5, 23.5, 25.5, 26.5, 30.5, 31.5, 
                          34.5, 36.5, 38.5, 39.5, 51.5, 53.5, 58.5,
                          65.5, 69.5, 72.5, 73.5, 74.5, 75.5, 80.5,
                          81.5))+
  theme(plot.margin = margin(t = 5, 
                             r = 5, 
                             b = 5, 
                             l = 5, 
                             unit = "pt"))
p1

labels_and_y <- list("ANGPT" = 1.5, 
                     "APP" = 3, 
                     "BMP" = 5, 
                     "CCL" = 9.5, 
                     "COMPLEMENT" = 14, 
                     "COLLAGEN" = 17, 
                     "CD46" = 19, 
                     "CD48" = 20, 
                     "CD86" = 21.5, 
                     "CD99" = 23, 
                     "DHEAS" = 24.5, 
                     "ESAM" = 26, 
                     "FGF" = 28.5, 
                     "GAS" = 31, 
                     "IGF" = 33, 
                     "IL1" = 35.5, 
                     "ICAM" = 37.5, 
                     "KLK" = 39, 
                     "LAMININ" = 45.5, 
                     "MIF" = 52.5, 
                     "MHC-I" = 56, 
                     "MHC-II" = 62, 
                     "NT" = 67.5, 
                     "NOTCH" = 71, 
                     "PDGF" = 73, 
                     "PLAU" = 74, 
                     "THBS" = 75, 
                     "VEGF" = 78, 
                     "VISFATIN" = 81, 
                     "VCAM" = 82)  
gene_grodown <- ggplot(aea.t.ldown,
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
  ##得到3.3 EC-AEA的target下调配体气泡图A

gg1 <- netVisual_bubble(cellchat, 
                        pairLR.use = aea.t.ldown.use, 
                        targets.use = "EC-AEA", 
                        sources.use = c("B", "T", "FIB", "MC", "POD", 
                                        "EC-AEA", "DC", "MAC", "VSMC"),
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "Down-regulated Ligand of EC-AEA in DKD")
gg1 <- gg1 + 
  geom_hline(yintercept=c(2.5, 3.5, 6.5, 12.5, 15.5, 18.5, 19.5, 
                          20.5, 22.5, 23.5, 25.5, 26.5, 30.5, 31.5, 
                          34.5, 36.5, 38.5, 39.5, 51.5, 53.5, 58.5,
                          65.5, 69.5, 72.5, 73.5, 74.5, 75.5, 80.5,
                          81.5))
gg1 %>% insert_left(gene_grodown, width = 0.3)
  ##得到3.4 EC-AEA的target下调配体气泡图B



####EC-AEA的下调受体气泡图（source）#####################
#EC-AEA的下调受体气泡图（source）
aea.s.rdown <- receptor.down %>%  
  filter(source == "EC-AEA")
aea.s.rdown <- aea.s.rdown %>% 
  filter(target %in% c("B", "T", "FIB", "MC", "POD", 
                       "EC-AEA", "DC", "MAC", "VSMC"))

aea.s.rdown <- aea.s.rdown %>%  
  arrange(substr(pathway_name, 1, 1)) 

unique(aea.s.rdown$pathway_name)
aea.s.rdown.use = aea.s.rdown[, "interaction_name", drop = F]
aea.s.rdown$interaction_name <- as.factor(aea.s.rdown$interaction_name)
aea.s.rdown$interaction_name <- fct_inorder(aea.s.rdown$interaction_name)
aea.s.rdown$interaction_name_2 <- as.factor(aea.s.rdown$interaction_name_2)
aea.s.rdown$interaction_name_2 <- fct_inorder(aea.s.rdown$interaction_name_2)
aea.s.rdown$pathway_name <- as.factor(aea.s.rdown$pathway_name)
aea.s.rdown$pathway_name <- fct_inorder(aea.s.rdown$pathway_name)
aea.s.rdown <- aea.s.rdown %>% 
  mutate(p="")
aea.s.rdown$signway <- paste(aea.s.rdown$source,
                          "  ->  ", 
                          aea.s.rdown$target) 

table(unique(aea.s.rdown$interaction_name_2))
summary_data <- aea.s.rdown %>%  
  group_by(pathway_name) %>%
  summarise(interaction_count = n_distinct(interaction_name_2))
print(summary_data)
write.csv(summary_data, "data.csv")

p1 <- ggplot(aea.s.rdown,
             aes(x = signway,
                 y = interaction_name_2)) + 
  geom_point(aes(size = pval, 
                 color = receptor.logFC)) +
  geom_point(shape = 1,
             aes(size = pval), 
             color="black")+
  scale_size(rang = c(1.5,6)) +
  labs(x=NULL,y=NULL,size="P value",
       title = "Down-regulated Receptor of EC-AEA in DKD") +
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
  geom_hline(yintercept=c(1.5, 3.5, 4.5, 5.5, 6.5, 8.5, 21.5, 22.5, 
                          23.5, 28.5, 32.5, 33.5, 37.5, 44.5, 45.5, 
                          47.5, 48.5, 49.5, 50.5, 51.5, 52.5, 54.5,
                          55.5, 77.5, 79.5, 82.5, 85.5, 88.5, 91.5,
                          92.5, 93.5, 94.5, 103.5, 104.5, 105.5))+
  theme(plot.margin = margin(t = 5, 
                             r = 5, 
                             b = 5, 
                             l = 5, 
                             unit = "pt"))
p1

labels_and_y <- list("ANGPT" = 1, 
                     "APP" = 2.5, 
                     "BAFF" = 4, 
                     "CXCL" = 5, 
                     "CD40" = 6, 
                     "COMPLEMENT" = 7.5, 
                     "COLLAGEN" = 15, 
                     "CD46" = 22, 
                     "CD99" = 23, 
                     "EPHA" = 26, 
                     "EPHB" = 30.5, 
                     "ESAM" = 33, 
                     "FGF" = 35.5, 
                     "FN1" = 41, 
                     "GH" = 45, 
                     "IGF" = 46.5, 
                     "IL4" = 48, 
                     "IL6" = 49, 
                     "IL1" = 50, 
                     "IFN-II" = 51, 
                     "IGFBP" = 52, 
                     "JAM" = 53.5, 
                     "LT" = 55, 
                     "LAMININ" = 66.5, 
                     "MSTN" = 78.5, 
                     "MIF" = 81, 
                     "MK" = 84, 
                     "NT" = 87, 
                     "PDGF" = 90, 
                     "PARs" = 92, 
                     "PLAU" = 93, 
                     "SEMA3" = 94, 
                     "TGFb" = 99, 
                     "THBS" = 104, 
                     "TENASCIN" = 105, 
                     "VTN" = 106.5) 
gene_grodown <- ggplot(aea.s.rdown,
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
  ##得到4.1 EC-AEA的source下调受体气泡图A


gg1 <- netVisual_bubble(cellchat, 
                        pairLR.use = aea.s.rdown.use, 
                        sources.use = "EC-AEA",
                        targets.use = c("B", "T", "FIB", "MC",
                                        "POD", "EC-AEA", "DC",
                                        "MAC", "VSMC"),
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "Down-regulated Receptor of EC-AEA in DKD")
gg1 <- gg1 + 
  geom_hline(yintercept=c(1.5, 3.5, 4.5, 5.5, 6.5, 8.5, 21.5, 22.5, 
                          23.5, 28.5, 32.5, 33.5, 37.5, 44.5, 45.5, 
                          47.5, 48.5, 49.5, 50.5, 51.5, 52.5, 54.5,
                          55.5, 77.5, 79.5, 82.5, 85.5, 88.5, 91.5,
                          92.5, 93.5, 94.5, 103.5, 104.5, 105.5))
gg1 %>% insert_left(gene_grodown, width = 0.3)
  ##得到4.2 EC-AEA的source下调受体气泡图B


####EC-AEA的下调受体气泡图（target）#####################
#EC-AEA的下调受体气泡图（target）
aea.t.rdown <- receptor.down %>%  
  filter(target == "EC-AEA")
aea.t.rdown <- aea.t.rdown %>% 
  filter(source %in% c("B", "T", "FIB", "MC", "POD", 
                       "EC-AEA", "DC", "MAC", "VSMC"))

aea.t.rdown <- aea.t.rdown %>%  
  arrange(substr(pathway_name, 1, 1)) 

unique(aea.t.rdown$pathway_name)
aea.t.rdown.use = aea.t.rdown[, "interaction_name", drop = F]
aea.t.rdown$interaction_name <- as.factor(aea.t.rdown$interaction_name)
aea.t.rdown$interaction_name <- fct_inorder(aea.t.rdown$interaction_name)
aea.t.rdown$interaction_name_2 <- as.factor(aea.t.rdown$interaction_name_2)
aea.t.rdown$interaction_name_2 <- fct_inorder(aea.t.rdown$interaction_name_2)
aea.t.rdown$pathway_name <- as.factor(aea.t.rdown$pathway_name)
aea.t.rdown$pathway_name <- fct_inorder(aea.t.rdown$pathway_name)
aea.t.rdown <- aea.t.rdown %>% 
  mutate(p="")

aea.t.rdown$signway <- paste(aea.t.rdown$source,
                          "  ->  ", 
                          aea.t.rdown$target) 

table(unique(aea.t.rdown$interaction_name_2))
summary_data <- aea.t.rdown %>%  
  group_by(pathway_name) %>%
  summarise(interaction_count = n_distinct(interaction_name_2))
print(summary_data)
write.csv(summary_data, "data.csv")

p1 <- ggplot(aea.t.rdown,
             aes(x = signway,
                 y = interaction_name_2)) + 
  geom_point(aes(size = pval, 
                 color = receptor.logFC)) +
  geom_point(shape = 1,
             aes(size = pval), 
             color="black")+
  scale_size(rang = c(1.5,6)) +
  labs(x=NULL,y=NULL,size="P value",
       title = "Down-regulated Receptor of EC-AEA in DKD") +
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
  geom_hline(yintercept=c(1.5, 2.5, 6.5, 8.5, 9.5, 14.5, 15.5, 
                          16.5, 17.5))+
  theme(plot.margin = margin(t = 5, 
                             r = 5, 
                             b = 5, 
                             l = 5, 
                             unit = "pt"))
p1

labels_and_y <- list("CD46" = 1, 
                     "CD99" = 2, 
                     "EPHA" = 4.5, 
                     "EPHB" = 7.5, 
                     "ESAM" = 9, 
                     "FGF" = 12, 
                     "LT" = 15, 
                     "MK" = 16, 
                     "PTN" = 17, 
                     "TWEAK" = 18)  
gene_grodown <- ggplot(aea.t.rdown,
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
  ##得到4.3 EC-AEA的target下调受体气泡图A

gg1 <- netVisual_bubble(cellchat, 
                        pairLR.use = aea.t.rdown.use, 
                        targets.use = "EC-AEA", 
                        sources.use = c("B", "T", "FIB", "MC",
                                        "POD", "EC-AEA", "DC",
                                        "MAC", "VSMC"),
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = FALSE,
                        title.name = "Down-regulated Receptor of EC-AEA in DKD")
gg1 <- gg1 + 
  geom_hline(yintercept=c(1.5, 2.5, 6.5, 8.5, 9.5, 14.5, 15.5,
                          16.5, 17.5))
gg1 %>% insert_left(gene_grodown, width = 0.3)
  ##得到4.4 EC-AEA的target下调受体气泡图B


