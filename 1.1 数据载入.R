rm(list = ls())
gc()



#打开必要的package
{
if(!require(multtest))BiocManager::install("multtest")
if(!require(Seurat))install.packages("Seurat")
if(!require(dplyr))install.packages("dplyr")
if(!require(tidyverse))install.packages("tidyverse")
}



setwd("D:/2023.10/调试/糖尿病肾病/1.人DKD分析")
getwd()
gc()



####载入数据#############################
#载入数据
setwd("./data source/raw_rds")
rds_files <- list()  
file_patterns <- c("CON", "DKD")  
numbers <- 1:6  
suffix <- "n.rds"  
for (pattern in file_patterns) {  
  for (num in numbers) {  
    filename <- paste0(pattern, num, suffix)  
    file_path <- file.path(filename) 
    
    if (file.exists(file_path)) {  
      var_name <- paste0(pattern, num)  
      assign(var_name, readRDS(file_path), envir = .GlobalEnv)  
    } else {  
      warning(paste("文件", filename, "不存在，请检查路径和文件名。"))  
    }  
  }  
}  



####写入病例信息#############################
#写入病例信息
columns_to_remove <- c("SpecimenID", "ClusterNumber", "ClusterClass", 
                       "LibraryID", "age", "seurat_clusters", "gender",
                       "state", "tissuetype", "celltype")  


##删除多余的meta.data  
for (i in 1:6) {  
  obj_name <- paste0("CON", i)  
 
  if (exists(obj_name, envir = .GlobalEnv) && is(get(obj_name, envir = .GlobalEnv), "Seurat")) {  
 
    seurat_obj <- get(obj_name, envir = .GlobalEnv)  
  
    seurat_obj@meta.data <- seurat_obj@meta.data[, !(names(seurat_obj@meta.data) %in% columns_to_remove), drop = FALSE]  
  
    assign(obj_name, seurat_obj, envir = .GlobalEnv)  
  } else {  
    warning(paste0("对象 '", obj_name, "' 不存在或不是Seurat对象。"))  
  }  
}  
for (i in 1:6) {  
  obj_name <- paste0("DKD", i)  
  
  if (exists(obj_name, envir = .GlobalEnv) && is(get(obj_name, envir = .GlobalEnv), "Seurat")) {  
 
    seurat_obj <- get(obj_name, envir = .GlobalEnv)  

    seurat_obj@meta.data <- seurat_obj@meta.data[, !(names(seurat_obj@meta.data) %in% columns_to_remove), drop = FALSE]  
 
    assign(obj_name, seurat_obj, envir = .GlobalEnv)  
  } else {  
    warning(paste0("对象 '", obj_name, "' 不存在或不是Seurat对象。"))  
  }  
}  


##修改orig.ident 
for (i in 1:6) {  
  con_obj_name <- paste0("CON", i)  
  
  if (exists(con_obj_name, envir = .GlobalEnv) && is(get(con_obj_name, envir = .GlobalEnv), "Seurat")) {  
    seurat_obj <- get(con_obj_name, envir = .GlobalEnv)  
    seurat_obj@meta.data$orig.ident <- con_obj_name  
    assign(con_obj_name, seurat_obj, envir = .GlobalEnv)  
  } else {  
    warning(paste0("对象 '", con_obj_name, "' 不存在或不是Seurat对象。"))  
  }  
}  
for (i in 1:6) {
  dkd_obj_name <- paste0("DKD", i)  

  if (exists(dkd_obj_name, envir = .GlobalEnv) && is(get(dkd_obj_name, envir = .GlobalEnv), "Seurat")) {  
    seurat_obj <- get(dkd_obj_name, envir = .GlobalEnv)  
    seurat_obj@meta.data$orig.ident <- dkd_obj_name  
    assign(dkd_obj_name, seurat_obj, envir = .GlobalEnv)  
  } else {  
    warning(paste0("对象 '", dkd_obj_name, "' 不存在或不是Seurat对象。"))  
  }  
}  


##增加性别信息
for (i in 1:3) {  
  obj_name <- paste0("CON", i)  
  obj_name_obj <- get(obj_name)  
  obj_name_obj@meta.data$sex <- "M" 
  assign(obj_name, obj_name_obj)  
}  
for (i in 1:3) {  
  obj_name <- paste0("DKD", i)  
  obj_name_obj <- get(obj_name) 
  obj_name_obj@meta.data$sex <- "M"    
  assign(obj_name, obj_name_obj)  
}
for (i in 4:6) {  
  obj_name <- paste0("CON", i)  
  obj_name_obj <- get(obj_name)  
  obj_name_obj@meta.data$sex <- "F" 
  assign(obj_name, obj_name_obj)  
}  
for (i in 4:6) {  
  obj_name <- paste0("DKD", i)  
  obj_name_obj <- get(obj_name) 
  obj_name_obj@meta.data$sex <- "F"    
  assign(obj_name, obj_name_obj)  
}



####合并对象################
#合并对象
seurat_list <- list(CON1, CON2, CON3, CON4, CON5, CON6, 
                    DKD1, DKD2, DKD3, DKD4, DKD5, DKD6)  
merged_human <- Reduce(function(x, y) {
  merge(x, y, add.cell.ids = list(names(x), names(y)))  
}, seurat_list)  
table(merged_human$orig.ident)

setwd("D:/2023.10/调试/糖尿病肾病/1.人DKD分析/data source/merge_rds")
saveRDS(merged_human, "全部样本（原始）.rds")
