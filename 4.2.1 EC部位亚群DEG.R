rm(list = ls())
gc()



#打开必要的package
{
  if(!require(Seurat))install.packages("Seurat")
  if(!require(EnhancedVolcano))BiocManager::install("EnhancedVolcano")
  if(!require(presto))devtools::install_github("immunogenomics/presto")
  if(!require(tidyverse))install.packages("tidyverse")
  if(!require(ggrepel))install.packages("ggrepel")
}



setwd("D:/2023.10/调试/糖尿病肾病/1.人DKD分析")
getwd()
human_all <- readRDS("./data source/merge_rds/全部样本（自定义注释）.rds")



####全部EC############################
#全部EC
EC <- subset(human_all, subclass.l1 == "EC")
table(EC$subclass.l1)
table(EC$usetype)
EC$usetype <- droplevels(EC$usetype,
                         exclude = setdiff(
                           levels(EC$usetype),
                           unique(EC$usetype)))
EC_DEG <- FindMarkers(EC, 
                      min.pct = 0.10, 
                      logfc.threshold = 0.10,
                      group.by = "sampletype",
                      ident.1 = "DKD",
                      ident.2 = "LD")
  ##这里ident.1是实验组，ident.2是对照组
write.csv(EC_DEG, "全部EC的DEG.csv")

EnhancedVolcano(EC_DEG,
                lab = rownames(EC_DEG),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 10e-5,
                FCcutoff = 1,
                pointSize = 1.0,
                labSize = 6.0,
                title = "全部EC的DEG（2倍）")
  ##得到1.全部EC的DEG（2倍）

EnhancedVolcano(EC_DEG,
                lab = rownames(EC_DEG),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 10e-5,
                FCcutoff = 2,
                pointSize = 1.0,
                labSize = 6.0,
                title = "全部EC的DEG（4倍）")
  ##得到2.全部EC的DEG（4倍）



####EC-PTC############################
#EC-PTC
EC_PTC <- subset(EC, usetype == "EC-PTC")
table(EC_PTC$usetype)
EC_PTC$usetype <- droplevels(EC_PTC$usetype,
                         exclude = setdiff(
                           levels(EC_PTC$usetype),
                           unique(EC_PTC$usetype)))
EC_PTC_DEG <- FindMarkers(EC_PTC, 
                          min.pct = 0.10, 
                          logfc.threshold = 0.10,
                          group.by = "sampletype",
                          ident.1 = "DKD",
                          ident.2 = "LD")
  ##这里ident.1是实验组，ident.2是对照组
write.csv(EC_PTC_DEG, "EC-PTC的DEG.csv")

EnhancedVolcano(EC_PTC_DEG,
                lab = rownames(EC_PTC_DEG),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 10e-5,
                FCcutoff = 1,
                pointSize = 1.0,
                labSize = 6.0,
                title = "EC-PTC的DEG（2倍）")
  ##得到3.EC-PTC的DEG（2倍）

EnhancedVolcano(EC_PTC_DEG,
                lab = rownames(EC_PTC_DEG),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 10e-5,
                FCcutoff = 2,
                pointSize = 1.0,
                labSize = 6.0,
                title = "EC-PTC的DEG（4倍）")
  ##得到4.EC-PTC的DEG（4倍）



####EC-GC############################
#EC-GC
EC_GC <- subset(EC, usetype == "EC-GC")
table(EC_GC$usetype)
EC_GC$usetype <- droplevels(EC_GC$usetype,
                             exclude = setdiff(
                               levels(EC_GC$usetype),
                               unique(EC_GC$usetype)))
EC_GC_DEG <- FindMarkers(EC_GC, 
                          min.pct = 0.10, 
                          logfc.threshold = 0.10,
                          group.by = "sampletype",
                          ident.1 = "DKD",
                          ident.2 = "LD")
##这里ident.1是实验组，ident.2是对照组
write.csv(EC_GC_DEG, "EC-GC的DEG.csv")

EnhancedVolcano(EC_GC_DEG,
                lab = rownames(EC_GC_DEG),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 10e-5,
                FCcutoff = 1,
                pointSize = 1.0,
                labSize = 6.0,
                title = "EC-GC的DEG（2倍）")
  ##得到5.EC-GC的DEG（2倍）

EnhancedVolcano(EC_GC_DEG,
                lab = rownames(EC_GC_DEG),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 10e-5,
                FCcutoff = 2,
                pointSize = 1.0,
                labSize = 6.0,
                title = "EC-GC的DEG（4倍）")
  ##得到6.EC-GC的DEG（4倍）



####EC-AEA############################
#EC-AEA
EC_AEA <- subset(EC, usetype == "EC-AEA")
table(EC_AEA$usetype)
EC_AEA$usetype <- droplevels(EC_AEA$usetype,
                            exclude = setdiff(
                              levels(EC_AEA$usetype),
                              unique(EC_AEA$usetype)))
EC_AEA_DEG <- FindMarkers(EC_AEA, 
                         min.pct = 0.10, 
                         logfc.threshold = 0.10,
                         group.by = "sampletype",
                         ident.1 = "DKD",
                         ident.2 = "LD")
##这里ident.1是实验组，ident.2是对照组
write.csv(EC_AEA_DEG, "EC-AEA的DEG.csv")

EnhancedVolcano(EC_AEA_DEG,
                lab = rownames(EC_AEA_DEG),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 10e-5,
                FCcutoff = 1,
                pointSize = 1.0,
                labSize = 6.0,
                title = "EC-AEA的DEG（2倍）")
  ##得到7.EC-AEA的DEG（2倍）

EnhancedVolcano(EC_AEA_DEG,
                lab = rownames(EC_AEA_DEG),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 10e-5,
                FCcutoff = 2,
                pointSize = 1.0,
                labSize = 6.0,
                title = "EC-AEA的DEG（4倍）")
  ##得到8.EC-AEA的DEG（4倍）



########################
#可视化
##数据准备
EC_DEG$celltype <- "all EC"
EC_PTC_DEG$celltype <- "EC-PTC"
EC_AEA_DEG$celltype <- "EC-AEA"
EC_GC_DEG$celltype <- "EC-GC"
diff_cell<-rbind(EC_DEG, EC_PTC_DEG, EC_AEA_DEG, EC_GC_DEG)
diff_cell$GENE <- rownames(diff_cell)
head(diff_cell)


##显著性
diff_cell$label <- ifelse(diff_cell$p_val_adj<0.01,"adjust P-val<0.01","adjust P-val>=0.01")
top10_EC <- filter(diff_cell,celltype=="all EC") %>% 
  distinct(GENE,.keep_all = T) %>% 
  top_n(10,abs(avg_log2FC))
top10_EC_PTC <- filter(diff_cell,celltype=="EC-PTC") %>% 
  distinct(GENE,.keep_all = T) %>% 
  top_n(10,abs(avg_log2FC))
top10_EC_AEA <- filter(diff_cell,celltype=="EC-AEA") %>% 
  distinct(GENE,.keep_all = T) %>% 
  top_n(10,abs(avg_log2FC))
top10_EC_GC <- filter(diff_cell,celltype=="EC-GC") %>% 
  distinct(GENE,.keep_all = T) %>% 
  top_n(10,abs(avg_log2FC))

top10 <- rbind(top10_EC, top10_EC_PTC, top10_EC_AEA, top10_EC_GC)
diff_cell$size <- case_when(!(diff_cell$GENE %in% top10$GENE)~ 1,
                            diff_cell$GENE %in% top10$GENE ~ 2)


##提取非Top10的基因表格；
dt <- filter(diff_cell,size==1)
head(dt)


##绘图
p <- ggplot()+
  geom_jitter(data = dt,
              aes(x = celltype, y = avg_log2FC, color = label),
              size = 0.85,
              width =0.4)+
  geom_jitter(data = top20,
              aes(x = celltype, y = avg_log2FC, color = label),
              size = 1,
              width =0.4)
p

dfbar<-data.frame(x=c("all EC", "EC-AEA","EC-GC","EC-PTC"),
                  y=c(12, 8, 12, 11))
dfbar1<-data.frame(x=c("all EC", "EC-AEA","EC-GC","EC-PTC"),
                   y=c(-15, -13, -13, -15))
p1 <- ggplot()+
  geom_col(data = dfbar,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)+
  geom_col(data = dfbar1,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)
p1

p2 <- ggplot()+
  geom_col(data = dfbar,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)+
  geom_col(data = dfbar1,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)+
  geom_jitter(data = dt,
              aes(x = celltype, y = avg_log2FC, color = label),
              size = 0.85,
              width =0.4)+
  geom_jitter(data = top10,
              aes(x = celltype, y = avg_log2FC, color = label),
              size = 1,
              width =0.4)
dfcol<-data.frame(x=c("all EC", "EC-AEA","EC-GC","EC-PTC"),
                  y=0,
                  label=c("all EC", "EC-AEA","EC-GC","EC-PTC"))
p2

mycol <- c("#E64B357F","#00A0877F","#34887F","#F39B7F7F")
p3 <- p2 + geom_tile(data = dfcol,
                     aes(x=x,y=y),
                     height=0.8,
                     color = "black",
                     fill = mycol,
                     alpha = 0.6,
                     show.legend = F)
p3
p4 <- p2 + geom_tile(data = dfcol,
                     aes(x=x,y=y),
                     height=0.8,
                     color = "black",
                     fill = mycol,
                     alpha = 0.6,
                     show.legend = F)+
  geom_text_repel(
    data=top10,
    aes(x=celltype,y=avg_log2FC,label=GENE),
    size =3,
    arrow = arrow(length = unit(0.008, "npc"),
                  type = "open", ends = "last"), 
    max.overlaps = 100
  )
p4
p5 <- p4+
  labs(x="celltype",y="avg_log2FC")+
  geom_text(data=dfcol,
            aes(x=x,y=y,label=label),
            size =3,
            color ="white")
p5
light_red <- rgb(250, 100, 100, maxColorValue = 255)
p6 <- p4 +
  scale_color_manual(name=NULL,
                     values = c(light_red,"grey"))+
  labs(x="celltype",y="avg_log2FC")+
  geom_text(data=dfcol,
            aes(x=x,y=y,label=label),
            size =4,
            color ="white")
p6
  ##得到9.内皮细胞DEG汇总
write.csv(diff_cell, "EC亚群DEG.csv")
