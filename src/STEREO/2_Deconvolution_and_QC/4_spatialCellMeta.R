library(tidyverse)
library(qs)
library(parallel)
library(magrittr)
library(RANN)
library(ggridges)
library(dendextend)
library(ggpubr)
devtools::load_all("~/spacexr-master/")
setwd("~/cortex/STEREO/2_Deconvolution_and_QC/")

chipList <- read.delim("~/cortex/STEREO/cortexMeta.txt") %>% {setNames(nm = .$chip,.$region)}

layer_anno <- list.files("../3_domainClustering/","_[2,1]00.csv",recursive = T,full.names = T)  %>% keep(~str_detect(.,str_c(names(chipList),collapse = "|")))
fpath <- layer_anno[[1]]

rctdRes <- list.files("../2_Deconvolution_and_QC/","qs",full.names = T) %>% keep(~str_detect(.,str_c(names(chipList),collapse = "|")))


anno <- mclapply(layer_anno,mc.cores = length(layer_anno),FUN = function(fpath){
  tmp <- read.csv(fpath)
  bin_size <-  basename(fpath) %>% str_extract("[1,2]00.*.csv") %>% str_extract("[1,2]00") %>% as.numeric()
  tmp <- tmp[colnames(tmp) %>% keep(~ str_detect(.,str_c("x","y","mainLayer","chip","domain",sep = "|")))]
  tmp <- tmp[colnames(tmp)[c(1,2,ncol(tmp))]]
  colnames(tmp) <- c("x","y","domain")
  tmp$bin_area <- (bin_size*0.5/1000)^2 # mm2
  tmp$chip <- basename(fpath) %>% str_extract("^[A-Z0-9]*")
  tmp
  }) %>% purrr::reduce(rbind)

anno$domain %<>% str_to_upper()
anno$domain %<>% str_remove(regex("[A,B,C]$",ignore_case = T))

region_to_depth = c(0,0.5,NA,1,2,3,4,5,6,7)
names(region_to_depth) <- c("ARACHNOID", "ARACHNOID L1", "EMPTY", "L1", "L2", "L3", "L4", 
                            "L5", "L6", "WM")
anno$depth <- region_to_depth[anno$domain]
anno %<>% group_by(chip) %>% group_split()
names(anno) <- purrr::map(anno,~.$chip[[1]]) %>% unlist


x = anno[[1]]
annoArea_total <- mclapply(anno,mc.cores = length(anno),FUN = function(x){
  x %>% group_by(chip) %>% summarise(area = sum(bin_area)) 
}) %>% purrr::reduce(rbind)

x = anno[[1]]
annoArea_layer <- mclapply(anno,mc.cores = length(anno),FUN = function(x){
  x %>% group_by(chip,domain) %>% summarise(area = sum(bin_area)) 
}) %>% purrr::reduce(rbind)

x = rctdRes[[1]]

spatialCellMeta <- mclapply(rctdRes,mc.cores = min(length(rctdRes),100),mc.preschedule = F,FUN = function(x){
  meta = basename(x) %>% str_split("_") %>% unlist
  ct = meta[[2]]
  chip = meta[[1]]
  seu <- qread(x)
  res <- cbind(seu@results$results_df,seu@spatialRNA@coords)
  res$chip <- chip
  res$region <- chipList[chip]
  tmp <- RANN::nn2(anno[[chip]][c("x","y")],query = res[c("x","y")],k = 1)
  res$domain <- anno[[chip]]$domain[tmp$nn.idx]
  res$depth <- anno[[chip]]$depth[tmp$nn.idx]
  res$subclass <- ct
  res$cluster <- res$first_type
  res
})# %>% purrr::reduce(rbind)

spatialCellMeta %>% keep(~!is.data.frame(.))
spatialCellMeta %<>% purrr::reduce(rbind)

spatialCellMeta <- spatialCellMeta[c("x", "y", "domain", "subclass", "cluster", "chip", "region", 
  "depth")]


OPC_MICRO = read_csv("OPC_MICRO.csv")
spatialCellMeta <- spatialCellMeta[colnames(OPC_MICRO)]
spatialCellMeta %<>% rbind(OPC_MICRO)

dend <- readRDS("../../../SnRNA/dend.RDS") 

spatialCellMeta$cluster %<>% str_replace("L2_3","L2/3") %>% str_replace("L5_6","L5/6") 

spatialCellMeta$cluster %>% unique %>% keep(~ !. %in% labels(dend)) %>% unique
# character(0)
labels(dend) %>% unique %>% keep(~ !. %in% unique(spatialCellMeta$cluster)) 
# character(0)

spatialCellMeta$cluster %<>% factor(levels = dend %>% labels)

write_csv(spatialCellMeta,"spatialCellMeta.csv")

