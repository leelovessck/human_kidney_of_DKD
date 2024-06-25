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
}



####数据准备############################
#数据准备
human_all <- readRDS("./cellchat的rds/全部细胞（seurat对象）.rds")



####创建control的celllchat对象#################################
#创建control的celllchat对象
control <- subset(human_all, sampletype == "LD")
table(control$usetype)
data.control <- GetAssayData(control, layer = "data")
meta.control <- control@meta.data
con.chat <- createCellChat(object = data.control,
                           meta = meta.control,
                           group.by = "usetype")
con.chat <- addMeta(con.chat, meta = meta.control)
con.chat <- setIdent(con.chat, ident.use = "usetype")  
groupSize <- as.numeric(table(con.chat@idents))



####载入全部数据库（control）################
#载入全部数据库（control）
CellChatDB <- CellChatDB.human  ##CellChatDB.mouse用于小鼠
showDatabaseCategory(CellChatDB)  ##展示数据库的基本组成
con.chat.all <- con.chat
con.chat.all@DB <- CellChatDB


##信号通路的预测
con.chat.all <- subsetData(con.chat.all,features = NULL)
con.chat.all <- identifyOverExpressedGenes(con.chat.all) 
con.chat.all <- identifyOverExpressedInteractions(con.chat.all) 
con.chat.all <- projectData(con.chat.all, PPI.human) 
con.chat.all <- computeCommunProb(con.chat.all, raw.use = F)
con.chat.all <- computeCommunProbPathway(con.chat.all)
con.chat.all <- aggregateNet(con.chat.all)


##保存预测结果
con.net.all <- subsetCommunication(con.chat.all)
write.csv(con.net.all,"1.1 control组CellChat结果（全部）.csv")
saveRDS(con.chat.all,"./cellchat的rds/control（全部）.rds")


##可视化
par(mfrow = c(1,2), xpd=TRUE)

netVisual_circle(con.chat.all@net$count,
                 vertex.weight = groupSize,
                 weight.scale = T,
                 label.edge= F, 
                 title.name = "Number of interactions")
netVisual_circle(con.chat.all@net$weight, 
                 vertex.weight = groupSize, 
                 weight.scale = T, 
                 label.edge= F, 
                 title.name = "Interaction weights/strength")
  ##得到1.1 control组全部通路图
  
mat <- con.chat.all@net$weight
par(mfrow = c(4,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, 
                 nrow = nrow(mat), 
                 ncol = ncol(mat),
                 dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, 
                   vertex.weight = groupSize, 
                   weight.scale = T, 
                   edge.weight.max = max(mat), 
                   title.name = rownames(mat)[i],
                   arrow.width = 0.1,
                   arrow.size = 0.02)
}
  ##得到1.2 control组各细胞全部通路图（权重）

mat <- con.chat.all@net$count
par(mfrow = c(4,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, 
                 nrow = nrow(mat), 
                 ncol = ncol(mat), 
                 dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, 
                   vertex.weight = groupSize, 
                   weight.scale = T, 
                   edge.weight.max = max(mat), 
                   title.name = rownames(mat)[i],
                   arrow.width = 0.1,
                   arrow.size = 0.02)
}
  ##得到1.3 control组各细胞全部通路图（数量）



####载入自、旁分泌（control）################
#载入自、旁分泌（control）
CellChatDB <- CellChatDB.human  ##CellChatDB.mouse用于小鼠
showDatabaseCategory(CellChatDB)  ##展示数据库的基本组成
CellChatDB.SS  <- subsetDB(CellChatDB, search = "Secreted Signaling")
con.chat.ss <- con.chat
con.chat.ss@DB <- CellChatDB.SS


##信号通路的预测
con.chat.ss <- subsetData(con.chat.ss,features = NULL)
con.chat.ss <- identifyOverExpressedGenes(con.chat.ss) 
con.chat.ss <- identifyOverExpressedInteractions(con.chat.ss) 
con.chat.ss <- projectData(con.chat.ss, PPI.human) 
con.chat.ss <- computeCommunProb(con.chat.ss, raw.use = F)
con.chat.ss <- computeCommunProbPathway(con.chat.ss)
con.chat.ss <- aggregateNet(con.chat.ss)


##保存预测结果
con.net.ss <- subsetCommunication(con.chat.ss)
write.csv(con.net.ss,"1.2 control组CellChat结果（自、旁分泌）.csv")
saveRDS(con.chat.ss,"./cellchat的rds/control（自、旁分泌）.rds")


##可视化
par(mfrow = c(1,2), xpd=TRUE)

netVisual_circle(con.chat.ss@net$count,
                 vertex.weight = groupSize,
                 weight.scale = T,
                 label.edge= F, 
                 title.name = "Number of interactions")
netVisual_circle(con.chat.ss@net$weight, 
                 vertex.weight = groupSize, 
                 weight.scale = T, 
                 label.edge= F, 
                 title.name = "Interaction weights/strength")
  ##得到2.1 control组自、旁分泌通路图

mat <- con.chat.ss@net$weight
par(mfrow = c(4,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, 
                 nrow = nrow(mat), 
                 ncol = ncol(mat),
                 dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, 
                   vertex.weight = groupSize, 
                   weight.scale = T, 
                   edge.weight.max = max(mat), 
                   title.name = rownames(mat)[i],
                   arrow.width = 0.1,
                   arrow.size = 0.02)
}
  ##得到2.2 control组各细胞自、旁分泌通路图（权重）

mat <- con.chat.ss@net$count
par(mfrow = c(4,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, 
                 nrow = nrow(mat), 
                 ncol = ncol(mat), 
                 dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, 
                   vertex.weight = groupSize, 
                   weight.scale = T, 
                   edge.weight.max = max(mat), 
                   title.name = rownames(mat)[i],
                   arrow.width = 0.1,
                   arrow.size = 0.02)
}
  ##得到2.3 control组各细胞自、旁分泌通路图（数量）



####载入细胞外基质（control）################
#载入细胞外基质（control）
CellChatDB <- CellChatDB.human  ##CellChatDB.mouse用于小鼠
showDatabaseCategory(CellChatDB)  ##展示数据库的基本组成
CellChatDB.ECM  <- subsetDB(CellChatDB, search = "ECM-Receptor")
con.chat.ecm <- con.chat
con.chat.ecm@DB <- CellChatDB.ECM


##信号通路的预测
con.chat.ecm <- subsetData(con.chat.ecm,features = NULL)
con.chat.ecm <- identifyOverExpressedGenes(con.chat.ecm) 
con.chat.ecm <- identifyOverExpressedInteractions(con.chat.ecm) 
con.chat.ecm <- projectData(con.chat.ecm, PPI.human) 
con.chat.ecm <- computeCommunProb(con.chat.ecm, raw.use = F)
con.chat.ecm <- computeCommunProbPathway(con.chat.ecm)
con.chat.ecm <- aggregateNet(con.chat.ecm)


##保存预测结果
con.net.ecm <- subsetCommunication(con.chat.ecm)
write.csv(con.net.ecm,"1.3 control组CellChat结果（细胞外基质）.csv")
saveRDS(con.chat.ecm,"./cellchat的rds/control（细胞外基质）.rds")


##可视化
par(mfrow = c(1,2), xpd=TRUE)

netVisual_circle(con.chat.ecm@net$count,
                 vertex.weight = groupSize,
                 weight.scale = T,
                 label.edge= F, 
                 title.name = "Number of interactions")
netVisual_circle(con.chat.ecm@net$weight, 
                 vertex.weight = groupSize, 
                 weight.scale = T, 
                 label.edge= F, 
                 title.name = "Interaction weights/strength")
  ##得到3.1 control组细胞外基质通路图

mat <- con.chat.ecm@net$weight
par(mfrow = c(4,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, 
                 nrow = nrow(mat), 
                 ncol = ncol(mat),
                 dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, 
                   vertex.weight = groupSize, 
                   weight.scale = T, 
                   edge.weight.max = max(mat), 
                   title.name = rownames(mat)[i],
                   arrow.width = 0.1,
                   arrow.size = 0.02)
}
  ##得到3.2 control组各细胞细胞外基质通路图（权重）

mat <- con.chat.ecm@net$count
par(mfrow = c(4,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, 
                 nrow = nrow(mat), 
                 ncol = ncol(mat), 
                 dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, 
                   vertex.weight = groupSize, 
                   weight.scale = T, 
                   edge.weight.max = max(mat), 
                   title.name = rownames(mat)[i],
                   arrow.width = 0.1,
                   arrow.size = 0.02)
}
  ##得到3.3 control组各细胞细胞外基质通路图（数量）



####载入细胞连接（control）################
#载入细胞连接（control）
CellChatDB <- CellChatDB.human  ##CellChatDB.mouse用于小鼠
showDatabaseCategory(CellChatDB)  ##展示数据库的基本组成
CellChatDB.CC  <- subsetDB(CellChatDB, search = "Cell-Cell Contact")
con.chat.cc <- con.chat
con.chat.cc@DB <- CellChatDB.CC


##信号通路的预测
con.chat.cc <- subsetData(con.chat.cc,features = NULL)
con.chat.cc <- identifyOverExpressedGenes(con.chat.cc) 
con.chat.cc <- identifyOverExpressedInteractions(con.chat.cc) 
con.chat.cc <- projectData(con.chat.cc, PPI.human) 
con.chat.cc <- computeCommunProb(con.chat.cc, raw.use = F)
con.chat.cc <- computeCommunProbPathway(con.chat.cc)
con.chat.cc <- aggregateNet(con.chat.cc)


##保存预测结果
con.net.cc <- subsetCommunication(con.chat.cc)
write.csv(con.net.cc,"1.4 control组CellChat结果（细胞连接）.csv")
saveRDS(con.chat.cc,"./cellchat的rds/control（细胞连接）.rds")


##可视化
par(mfrow = c(1,2), xpd=TRUE)

netVisual_circle(con.chat.cc@net$count,
                 vertex.weight = groupSize,
                 weight.scale = T,
                 label.edge= F, 
                 title.name = "Number of interactions")
netVisual_circle(con.chat.cc@net$weight, 
                 vertex.weight = groupSize, 
                 weight.scale = T, 
                 label.edge= F, 
                 title.name = "Interaction weights/strength")
  ##得到4.1 control组细胞连接通路图

mat <- con.chat.cc@net$weight
par(mfrow = c(4,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, 
                 nrow = nrow(mat), 
                 ncol = ncol(mat),
                 dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, 
                   vertex.weight = groupSize, 
                   weight.scale = T, 
                   edge.weight.max = max(mat), 
                   title.name = rownames(mat)[i],
                   arrow.width = 0.1,
                   arrow.size = 0.02)
}
  ##得到4.2 control组各细胞细胞连接通路图（权重）

mat <- con.chat.cc@net$count
par(mfrow = c(4,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, 
                 nrow = nrow(mat), 
                 ncol = ncol(mat), 
                 dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, 
                   vertex.weight = groupSize, 
                   weight.scale = T, 
                   edge.weight.max = max(mat), 
                   title.name = rownames(mat)[i],
                   arrow.width = 0.1,
                   arrow.size = 0.02)
}
  ##得到4.3 control组各细胞细胞连接通路图（数量）



####载入非蛋白（control）################
#载入非蛋白（control）
CellChatDB <- CellChatDB.human  ##CellChatDB.mouse用于小鼠
showDatabaseCategory(CellChatDB)  ##展示数据库的基本组成
CellChatDB.NP  <- subsetDB(CellChatDB, search = "Non-protein Signaling")
con.chat.np <- con.chat
con.chat.np@DB <- CellChatDB.NP


##信号通路的预测
con.chat.np <- subsetData(con.chat.np,features = NULL)
con.chat.np <- identifyOverExpressedGenes(con.chat.np) 
con.chat.np <- identifyOverExpressedInteractions(con.chat.np) 
con.chat.np <- projectData(con.chat.np, PPI.human) 
con.chat.np <- computeCommunProb(con.chat.np, raw.use = F)
con.chat.np <- computeCommunProbPathway(con.chat.np)
con.chat.np <- aggregateNet(con.chat.np)


##保存预测结果
con.net.np <- subsetCommunication(con.chat.np)
write.csv(con.net.np,"1.5 control组CellChat结果（非蛋白通路）.csv")
saveRDS(con.chat.np,"./cellchat的rds/control（非蛋白通路）.rds")


##可视化
par(mfrow = c(1,2), xpd=TRUE)

netVisual_circle(con.chat.np@net$count,
                 vertex.weight = groupSize,
                 weight.scale = T,
                 label.edge= F, 
                 title.name = "Number of interactions")
netVisual_circle(con.chat.np@net$weight, 
                 vertex.weight = groupSize, 
                 weight.scale = T, 
                 label.edge= F, 
                 title.name = "Interaction weights/strength")
  ##得到5.1 control组非蛋白通路图

mat <- con.chat.np@net$weight
par(mfrow = c(4,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, 
                 nrow = nrow(mat), 
                 ncol = ncol(mat),
                 dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, 
                   vertex.weight = groupSize, 
                   weight.scale = T, 
                   edge.weight.max = max(mat), 
                   title.name = rownames(mat)[i],
                   arrow.width = 0.1,
                   arrow.size = 0.02)
}
  ##得到5.2 control组各细胞非蛋白通路图（权重）

mat <- con.chat.np@net$count
par(mfrow = c(4,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, 
                 nrow = nrow(mat), 
                 ncol = ncol(mat), 
                 dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, 
                   vertex.weight = groupSize, 
                   weight.scale = T, 
                   edge.weight.max = max(mat), 
                   title.name = rownames(mat)[i],
                   arrow.width = 0.1,
                   arrow.size = 0.02)
}
  ##得到5.3 control组各细胞非蛋白通路图（数量）



####创建DKD的celllchat对象#################################
#创建DKD的celllchat对象
DKD <- subset(human_all, sampletype == "DKD")
table(DKD$usetype)
data.DKD <- GetAssayData(DKD, layer = "data")
meta.DKD <- DKD@meta.data
dkd.chat <- createCellChat(object = data.DKD,
                           meta = meta.DKD,
                           group.by = "usetype")
dkd.chat <- addMeta(dkd.chat, meta = meta.DKD)
dkd.chat <- setIdent(dkd.chat, ident.use = "usetype")  
groupSize <- as.numeric(table(dkd.chat@idents))



####载入全部数据库（DKD）################
#载入全部数据库（DKD）
CellChatDB <- CellChatDB.human  ##CellChatDB.mouse用于小鼠
showDatabaseCategory(CellChatDB)  ##展示数据库的基本组成
dkd.chat.all <- dkd.chat
dkd.chat.all@DB <- CellChatDB


##信号通路的预测
dkd.chat.all <- subsetData(dkd.chat.all,features = NULL)
dkd.chat.all <- identifyOverExpressedGenes(dkd.chat.all) 
dkd.chat.all <- identifyOverExpressedInteractions(dkd.chat.all) 
dkd.chat.all <- projectData(dkd.chat.all, PPI.human) 
dkd.chat.all <- computeCommunProb(dkd.chat.all, raw.use = F)
dkd.chat.all <- computeCommunProbPathway(dkd.chat.all)
dkd.chat.all <- aggregateNet(dkd.chat.all)


##保存预测结果
dkd.net.all <- subsetCommunication(dkd.chat.all)
write.csv(dkd.net.all,"2.1 DKD组CellChat结果（全部）.csv")
saveRDS(dkd.chat.all,"./cellchat的rds/DKD（全部）.rds")


##可视化
par(mfrow = c(1,2), xpd=TRUE)

netVisual_circle(dkd.chat.all@net$count,
                 vertex.weight = groupSize,
                 weight.scale = T,
                 label.edge= F, 
                 title.name = "Number of interactions")
netVisual_circle(dkd.chat.all@net$weight, 
                 vertex.weight = groupSize, 
                 weight.scale = T, 
                 label.edge= F, 
                 title.name = "Interaction weights/strength")
  ##得到6.1 DKD组全部通路图

mat <- dkd.chat.all@net$weight
par(mfrow = c(4,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, 
                 nrow = nrow(mat), 
                 ncol = ncol(mat),
                 dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, 
                   vertex.weight = groupSize, 
                   weight.scale = T, 
                   edge.weight.max = max(mat), 
                   title.name = rownames(mat)[i],
                   arrow.width = 0.1,
                   arrow.size = 0.02)
}
  ##得到6.2 DKD组各细胞全部通路图（权重）

mat <- dkd.chat.all@net$count
par(mfrow = c(4,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, 
                 nrow = nrow(mat), 
                 ncol = ncol(mat), 
                 dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, 
                   vertex.weight = groupSize, 
                   weight.scale = T, 
                   edge.weight.max = max(mat), 
                   title.name = rownames(mat)[i],
                   arrow.width = 0.1,
                   arrow.size = 0.02)
}
  ##得到6.3 DKD组各细胞全部通路图（数量）



####载入自、旁分泌（DKD）################
#载入自、旁分泌（DKD）
CellChatDB <- CellChatDB.human  ##CellChatDB.mouse用于小鼠
showDatabaseCategory(CellChatDB)  ##展示数据库的基本组成
CellChatDB.SS  <- subsetDB(CellChatDB, search = "Secreted Signaling")
dkd.chat.ss <- dkd.chat
dkd.chat.ss@DB <- CellChatDB.SS


##信号通路的预测
dkd.chat.ss <- subsetData(dkd.chat.ss,features = NULL)
dkd.chat.ss <- identifyOverExpressedGenes(dkd.chat.ss) 
dkd.chat.ss <- identifyOverExpressedInteractions(dkd.chat.ss) 
dkd.chat.ss <- projectData(dkd.chat.ss, PPI.human) 
dkd.chat.ss <- computeCommunProb(dkd.chat.ss, raw.use = F)
dkd.chat.ss <- computeCommunProbPathway(dkd.chat.ss)
dkd.chat.ss <- aggregateNet(dkd.chat.ss)


##保存预测结果
dkd.net.ss <- subsetCommunication(dkd.chat.ss)
write.csv(dkd.net.ss,"2.2 DKD组CellChat结果（自、旁分泌）.csv")
saveRDS(dkd.chat.ss,"./cellchat的rds/DKD（自、旁分泌）.rds")


##可视化
par(mfrow = c(1,2), xpd=TRUE)

netVisual_circle(dkd.chat.ss@net$count,
                 vertex.weight = groupSize,
                 weight.scale = T,
                 label.edge= F, 
                 title.name = "Number of interactions")
netVisual_circle(dkd.chat.ss@net$weight, 
                 vertex.weight = groupSize, 
                 weight.scale = T, 
                 label.edge= F, 
                 title.name = "Interaction weights/strength")
  ##得到7.1 DKD组自、旁分泌通路图

mat <- dkd.chat.ss@net$weight
par(mfrow = c(4,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, 
                 nrow = nrow(mat), 
                 ncol = ncol(mat),
                 dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, 
                   vertex.weight = groupSize, 
                   weight.scale = T, 
                   edge.weight.max = max(mat), 
                   title.name = rownames(mat)[i],
                   arrow.width = 0.1,
                   arrow.size = 0.02)
}
  ##得到7.2 DKD组各细胞自、旁分泌通路图（权重）

mat <- dkd.chat.ss@net$count
par(mfrow = c(4,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, 
                 nrow = nrow(mat), 
                 ncol = ncol(mat), 
                 dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, 
                   vertex.weight = groupSize, 
                   weight.scale = T, 
                   edge.weight.max = max(mat), 
                   title.name = rownames(mat)[i],
                   arrow.width = 0.1,
                   arrow.size = 0.02)
}
  ##得到7.3 DKD组各细胞自、旁分泌通路图（数量）



####载入细胞外基质（DKD）################
#载入细胞外基质（DKD）
CellChatDB <- CellChatDB.human  ##CellChatDB.mouse用于小鼠
showDatabaseCategory(CellChatDB)  ##展示数据库的基本组成
CellChatDB.ECM  <- subsetDB(CellChatDB, search = "ECM-Receptor")
dkd.chat.ecm <- dkd.chat
dkd.chat.ecm@DB <- CellChatDB.ECM


##信号通路的预测
dkd.chat.ecm <- subsetData(dkd.chat.ecm,features = NULL)
dkd.chat.ecm <- identifyOverExpressedGenes(dkd.chat.ecm) 
dkd.chat.ecm <- identifyOverExpressedInteractions(dkd.chat.ecm) 
dkd.chat.ecm <- projectData(dkd.chat.ecm, PPI.human) 
dkd.chat.ecm <- computeCommunProb(dkd.chat.ecm, raw.use = F)
dkd.chat.ecm <- computeCommunProbPathway(dkd.chat.ecm)
dkd.chat.ecm <- aggregateNet(dkd.chat.ecm)


##保存预测结果
dkd.net.ecm <- subsetCommunication(dkd.chat.ecm)
write.csv(dkd.net.ecm,"2.3 DKD组CellChat结果（细胞外基质）.csv")
saveRDS(dkd.chat.ecm,"./cellchat的rds/DKD（细胞外基质）.rds")


##可视化
par(mfrow = c(1,2), xpd=TRUE)

netVisual_circle(dkd.chat.ecm@net$count,
                 vertex.weight = groupSize,
                 weight.scale = T,
                 label.edge= F, 
                 title.name = "Number of interactions")
netVisual_circle(dkd.chat.ecm@net$weight, 
                 vertex.weight = groupSize, 
                 weight.scale = T, 
                 label.edge= F, 
                 title.name = "Interaction weights/strength")
  ##得到8.1 DKD组细胞外基质通路图

mat <- dkd.chat.ecm@net$weight
par(mfrow = c(4,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, 
                 nrow = nrow(mat), 
                 ncol = ncol(mat),
                 dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, 
                   vertex.weight = groupSize, 
                   weight.scale = T, 
                   edge.weight.max = max(mat), 
                   title.name = rownames(mat)[i],
                   arrow.width = 0.1,
                   arrow.size = 0.02)
}
  ##得到8.2 DKD组各细胞细胞外基质通路图（权重）

mat <- dkd.chat.ecm@net$count
par(mfrow = c(4,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, 
                 nrow = nrow(mat), 
                 ncol = ncol(mat), 
                 dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, 
                   vertex.weight = groupSize, 
                   weight.scale = T, 
                   edge.weight.max = max(mat), 
                   title.name = rownames(mat)[i],
                   arrow.width = 0.1,
                   arrow.size = 0.02)
}
  ##得到8.3 DKD组各细胞细胞外基质通路图（数量）



####载入细胞连接（DKD）################
#载入细胞连接（DKD）
CellChatDB <- CellChatDB.human  ##CellChatDB.mouse用于小鼠
showDatabaseCategory(CellChatDB)  ##展示数据库的基本组成
CellChatDB.CC  <- subsetDB(CellChatDB, search = "Cell-Cell Contact")
dkd.chat.cc <- dkd.chat
dkd.chat.cc@DB <- CellChatDB.CC


##信号通路的预测
dkd.chat.cc <- subsetData(dkd.chat.cc,features = NULL)
dkd.chat.cc <- identifyOverExpressedGenes(dkd.chat.cc) 
dkd.chat.cc <- identifyOverExpressedInteractions(dkd.chat.cc) 
dkd.chat.cc <- projectData(dkd.chat.cc, PPI.human) 
dkd.chat.cc <- computeCommunProb(dkd.chat.cc, raw.use = F)
dkd.chat.cc <- computeCommunProbPathway(dkd.chat.cc)
dkd.chat.cc <- aggregateNet(dkd.chat.cc)


##保存预测结果
dkd.net.cc <- subsetCommunication(dkd.chat.cc)
write.csv(dkd.net.cc,"2.4 DKD组CellChat结果（细胞连接）.csv")
saveRDS(dkd.chat.cc,"./cellchat的rds/DKD（细胞连接）.rds")


##可视化
par(mfrow = c(1,2), xpd=TRUE)

netVisual_circle(dkd.chat.cc@net$count,
                 vertex.weight = groupSize,
                 weight.scale = T,
                 label.edge= F, 
                 title.name = "Number of interactions")
netVisual_circle(dkd.chat.cc@net$weight, 
                 vertex.weight = groupSize, 
                 weight.scale = T, 
                 label.edge= F, 
                 title.name = "Interaction weights/strength")
  ##得到9.1 DKD组细胞连接通路图

mat <- dkd.chat.cc@net$weight
par(mfrow = c(4,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, 
                 nrow = nrow(mat), 
                 ncol = ncol(mat),
                 dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, 
                   vertex.weight = groupSize, 
                   weight.scale = T, 
                   edge.weight.max = max(mat), 
                   title.name = rownames(mat)[i],
                   arrow.width = 0.1,
                   arrow.size = 0.02)
}
  ##得到9.2 DKD组各细胞细胞连接通路图（权重）

mat <- dkd.chat.cc@net$count
par(mfrow = c(4,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, 
                 nrow = nrow(mat), 
                 ncol = ncol(mat), 
                 dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, 
                   vertex.weight = groupSize, 
                   weight.scale = T, 
                   edge.weight.max = max(mat), 
                   title.name = rownames(mat)[i],
                   arrow.width = 0.1,
                   arrow.size = 0.02)
}
  ##得到9.3 DKD组各细胞细胞连接通路图（数量）



####载入非蛋白（DKD）################
#载入非蛋白（DKD）
CellChatDB <- CellChatDB.human  ##CellChatDB.mouse用于小鼠
showDatabaseCategory(CellChatDB)  ##展示数据库的基本组成
CellChatDB.NP  <- subsetDB(CellChatDB, search = "Non-protein Signaling")
dkd.chat.np <- dkd.chat
dkd.chat.np@DB <- CellChatDB.NP


##信号通路的预测
dkd.chat.np <- subsetData(dkd.chat.np,features = NULL)
dkd.chat.np <- identifyOverExpressedGenes(dkd.chat.np) 
dkd.chat.np <- identifyOverExpressedInteractions(dkd.chat.np) 
dkd.chat.np <- projectData(dkd.chat.np, PPI.human) 
dkd.chat.np <- computeCommunProb(dkd.chat.np, raw.use = F)
dkd.chat.np <- computeCommunProbPathway(dkd.chat.np)
dkd.chat.np <- aggregateNet(dkd.chat.np)


##保存预测结果
dkd.net.np <- subsetCommunication(dkd.chat.np)
write.csv(dkd.net.np,"2.5 DKD组CellChat结果（非蛋白通路）.csv")
saveRDS(dkd.chat.np,"./cellchat的rds/DKD（非蛋白通路）.rds")


##可视化
par(mfrow = c(1,2), xpd=TRUE)

netVisual_circle(dkd.chat.np@net$count,
                 vertex.weight = groupSize,
                 weight.scale = T,
                 label.edge= F, 
                 title.name = "Number of interactions")
netVisual_circle(dkd.chat.np@net$weight, 
                 vertex.weight = groupSize, 
                 weight.scale = T, 
                 label.edge= F, 
                 title.name = "Interaction weights/strength")
  ##得到10.1 DKD组非蛋白通路图

mat <- dkd.chat.np@net$weight
par(mfrow = c(4,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, 
                 nrow = nrow(mat), 
                 ncol = ncol(mat),
                 dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, 
                   vertex.weight = groupSize, 
                   weight.scale = T, 
                   edge.weight.max = max(mat), 
                   title.name = rownames(mat)[i],
                   arrow.width = 0.1,
                   arrow.size = 0.02)
}
  ##得到10.2 DKD组各细胞非蛋白通路图（权重）

mat <- dkd.chat.np@net$count
par(mfrow = c(4,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, 
                 nrow = nrow(mat), 
                 ncol = ncol(mat), 
                 dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, 
                   vertex.weight = groupSize, 
                   weight.scale = T, 
                   edge.weight.max = max(mat), 
                   title.name = rownames(mat)[i],
                   arrow.width = 0.1,
                   arrow.size = 0.02)
}
  ##得到10.3 DKD组各细胞非蛋白通路图（数量）
