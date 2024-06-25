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
if(!require(ggalluvial))install.packages("ggalluvial")
if(!require(reticulate))install.packages("reticulate")
}



####数据载入############################
#数据载入
con.chat <- readRDS("./cellchat的rds/control（自、旁分泌）.rds")
dkd.chat <- readRDS("./cellchat的rds/DKD（自、旁分泌）.rds")




####control的outgoing模式#####################
#control的outgoing模式
selectK(con.chat, pattern = "outgoing")
  ##得到1.1 control的outgoing模式数量图（旁、自分泌）

con.chat <- identifyCommunicationPatterns(con.chat, 
                                          pattern = "outgoing", 
                                          k = 4,
                                          width = 15,
                                          height = 20)
  ##得到1.2 control的outgoing模式热图-模式（旁、自分泌）

con.chat <- identifyCommunicationPatterns(con.chat, 
                                          pattern = "outgoing", 
                                          k = 4,
                                          width = 5,
                                          height = 5)
  ##得到1.3 control的outgoing模式热图-细胞（旁、自分泌）

netAnalysis_river(con.chat,
                  pattern = "outgoing")
  ##得到1.4 control的outgoing模式冲击图-模式（旁、自分泌）
  ##得到1.5 control的outgoing模式冲击图-细胞（旁、自分泌）

netAnalysis_dot(con.chat,
                pattern = "outgoing")
  ##得到1.6 control的outgoing模式气泡图（旁、自分泌）



####control的incoming模式#####################
#control的incoming模式
selectK(con.chat, pattern = "incoming")
  ##得到2.1 control的incoming模式数量图（旁、自分泌）

con.chat <- identifyCommunicationPatterns(con.chat, 
                                          pattern = "incoming", 
                                          k = 4,
                                          width = 15,
                                          height = 20)
  ##得到2.2 control的incoming模式热图-模式（旁、自分泌）

con.chat <- identifyCommunicationPatterns(con.chat, 
                                          pattern = "incoming", 
                                          k = 4,
                                          width = 5,
                                          height = 5)
  ##得到2.3 control的incoming模式热图-细胞（旁、自分泌）

netAnalysis_river(con.chat,
                  pattern = "incoming")
  ##得到2.4 control的incoming模式冲击图-模式（旁、自分泌）
  ##得到2.5 control的incoming模式冲击图-细胞（旁、自分泌）

netAnalysis_dot(con.chat,
                pattern = "incoming")
  ##得到2.6 control的incoming模式气泡图（旁、自分泌）



####control的信号组########################################
#control的信号组
##基于功能
use_python("D:/anaconda3/python.exe", required = TRUE)
py_install("umap-learn")
con.chat <- computeNetSimilarity(con.chat, type = "functional")
con.chat <- netEmbedding(con.chat, type = "functional")
con.chat <- netClustering(con.chat, type = "functional")
netVisual_embedding(con.chat, type = "functional", label.size = 3.5)+
  ggtitle("Group by Function")
  ##得到3.1 control的功能信号组（旁、自分泌）


##基于结构
con.chat <- computeNetSimilarity(con.chat, type = "structural")
con.chat <- netEmbedding(con.chat, type = "structural")
con.chat <- netClustering(con.chat, type = "structural")
netVisual_embedding(con.chat, type = "structural", label.size = 3.5)+
  ggtitle("Group by Structural")
  ##得到3.2 control的结构信号组（旁、自分泌）



####DKD的outgoing模式#####################
#DKD的outgoing模式
selectK(dkd.chat, pattern = "outgoing")
  ##得到4.1 DKD的outgoing模式数量图（旁、自分泌）

dkd.chat <- identifyCommunicationPatterns(dkd.chat, 
                                          pattern = "outgoing", 
                                          k = 3,
                                          width = 15,
                                          height = 20)
  ##得到4.2 DKD的outgoing模式热图-模式（旁、自分泌）

dkd.chat <- identifyCommunicationPatterns(dkd.chat, 
                                          pattern = "outgoing", 
                                          k = 3,
                                          width = 5,
                                          height = 5)
  ##得到4.3 DKD的outgoing模式热图-细胞（旁、自分泌）

netAnalysis_river(dkd.chat,
                  pattern = "outgoing")
  ##得到4.4 DKD的outgoing模式冲击图-模式（旁、自分泌）
  ##得到4.5 DKD的outgoing模式冲击图-模式（旁、自分泌）

netAnalysis_dot(dkd.chat,
                pattern = "outgoing")
  ##得到4.6 DKD的outgoing模式气泡图（旁、自分泌）



####DKD的incoming模式#####################
#DKD的incoming模式
selectK(dkd.chat, pattern = "incoming")
  ##得到5.1 DKD的incoming模式数量图（旁、自分泌）

dkd.chat <- identifyCommunicationPatterns(dkd.chat, 
                                          pattern = "incoming", 
                                          k = 3,
                                          width = 15,
                                          height = 20)
  ##得到5.2 DKD的incoming模式热图-模式（旁、自分泌）

dkd.chat <- identifyCommunicationPatterns(dkd.chat, 
                                          pattern = "incoming", 
                                          k = 3,
                                          width = 5,
                                          height = 5)
  ##得到5.3 DKD的incoming模式热图-细胞（旁、自分泌）

netAnalysis_river(dkd.chat,
                  pattern = "incoming")
  ##得到5.4 DKD的incoming模式冲击图-模式（旁、自分泌）
  ##得到5.5 DKD的incoming模式冲击图-细胞（旁、自分泌）

netAnalysis_dot(dkd.chat,
                pattern = "incoming")
  ##得到5.6 DKD的incoming模式气泡图（旁、自分泌）



####DKD的信号组########################################
#DKD的信号组
##基于功能
dkd.chat <- computeNetSimilarity(dkd.chat, type = "functional")
dkd.chat <- netEmbedding(dkd.chat, type = "functional")
dkd.chat <- netClustering(dkd.chat, type = "functional")
netVisual_embedding(dkd.chat, type = "functional", label.size = 3.5)+
  ggtitle("Group by Function")
  ##得到6.1 DKD的功能信号组（旁、自分泌）


##基于结构
dkd.chat <- computeNetSimilarity(dkd.chat, type = "structural")
dkd.chat <- netEmbedding(dkd.chat, type = "structural")
dkd.chat <- netClustering(dkd.chat, type = "structural")
netVisual_embedding(dkd.chat, type = "structural", label.size = 3.5)+
  ggtitle("Group by Structural")
  ##得到6.2 DKD的结构信号组（旁、自分泌）
