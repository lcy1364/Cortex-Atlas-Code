library(purrr)
library(magrittr)
library(tidyverse)
library(Seurat)
library(harmony)
library(ape)
library(uwot)
library(ggtree)
library(treeio)
# library(future)
setwd("~/cortex/SnRNA/3_mergingDatasets/") 
set.seed(123)

qsFiles <- list.files(".","merge.qs",full.names = T)

sampledCell <- function(x){
  set.seed(1)
  x[,sample(1:ncol(x),size = ceiling(ncol(x)/3),replace = F)]
}


seu <- parallel::mclapply(qsFiles,mc.cores = length(qsFiles),FUN = function(x){
  tmp = qs::qread(x)
  tmp = tmp[,tmp$region != "MTG"]
  tmp = tmp[, (tmp$source != "lein_10x_layer5_only"| is.na(tmp$source)) & tmp$region != "MTG"] 
  tmp %>% sampledCell %>% DietSeurat()
})  %>% {merge(.[[1]],.[-c(1)])}

Meta <- read.csv("./SnRNA_Meta.csv",row.names = 1)

seu <- readRDS("../SnRNA_seurat.RDS")

# colnames(seu)[!colnames(seu) %in% rownames(Meta)]

table(colnames(seu) %in% rownames(Meta))

seu@meta.data <- Meta[colnames(seu),]

seu %<>% FindVariableFeatures(nfeatures = 1000) %>% ScaleData(split.by = "batch") %>% RunPCA(verbose = T) 

seu %<>% RunHarmony(group.by.vars = "batch",ncores = 10,dims.use = 1:40) #%>% 

seu %<>% RunUMAP(reduction = "harmony",reduction.name = "harmony_umap",reduction.key = "harmonyumap_",dims = 1:25) 

{DimPlot(seu,group.by = "subclass",split.by = "batch",label = T) + theme(legend.position = "none")} %>% ggsave(filename = "4_merge_seu.png",width = 10,height = 5, dpi = 150)

seu %>% saveRDS("../SnRNA_seurat.RDS", compress = "gzip")

#
x = qsFiles[[1]]
seu <- parallel::mclapply(qsFiles,mc.cores = length(qsFiles),FUN = function(x){
  tmp = qs::qread(x)
  tmp = tmp[,tmp$region != "MTG"]
  tmp = tmp[, (tmp$source != "lein_10x_layer5_only"| is.na(tmp$source)) & tmp$region != "MTG"] 
  tmp@meta.data <- Meta[colnames(tmp),]
  tmp %>% qs::qsave(x)
})
