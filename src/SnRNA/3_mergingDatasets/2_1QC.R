library(Seurat)
library(tidyverse)
library(parallel)
library(magrittr)
library(ggtree)
library(ape)
library(patchwork)
library(scrattch.hicat)
setwd("~/cortex/SnRNA/3_mergingDatasets//")

seu <- qs::qread("~/cortex/SnRNA/3_mergingDatasets/SnRNA_seurat.qs") 
# seu <- qs::qread("~/cortex/SnRNA/1_SnRNA_preprocessing/SnRNA_seurat.qs")

seu_list <- seu %>% SplitObject(split.by = "library_prep")

x = seu_list[[1]]
seu_list <- mclapply(seu_list,mc.cores = length(seu_list),FUN = function(x){
  x %<>% NormalizeData(scale.factor = 10^6) %>% FindVariableFeatures() %>% ScaleData(features = rownames(x)) %>% RunPCA() %>% RunUMAP(dims = 1:40)
})

tmp <- seu_list[[1]]

genes <- read.delim("../2_codePreprocessingExternalData/gene.txt",header = F)

mclapply(seu_list,mc.cores = length(seu_list),FUN = function(tmp){
  plot1 <- DimPlot(tmp,group.by = "subclass",label = T,repel = T,shuffle = T,cols = subclass_color <-c(AST = "#665C47", ENDO = "#604B47", ET = "#CEC823", CHANDELIER = "#E25691", 
                                                                        `L2-L3 IT LINC00507` = "#07D8D8", `L3-L4 IT RORB` = "#09B2B2", 
                                                                        `L4-L5 IT RORB` = "#69B199", `L6 CAR3` = "#898325", `L6 CT` = "#2D8CB8", 
                                                                        `L6 IT` = "#19E3BE", L6b = "#7944AA", LAMP5 = "#D96C07", MICRO = "#513577", 
                                                                        NDNF = "#f2798d", NP = "#93C43B", OLIGO = "#5E8A79", OPC = "#2E3E39", 
                                                                        PAX6 = "#E96DF2", PVALB = "#FF2D4E", SST = "#F2B003", VIP = "#9C28CC", 
                                                                        VLMC = "#697255")) + theme(legend.position = "none") + coord_fixed()
  
  
  aver_exp <- AverageExpression(tmp,group.by = "subclass",features = FindVariableFeatures(tmp,nfeatures = 3000) %>% VariableFeatures()) 
  
  tree <- aver_exp[[1]] %>% pvclust::pvclust() 
  
  plot2 <- ggtree(tree) + geom_tiplab() 
  
  
  tmp@misc$plot <- plot2
  
  tmp2 <- tmp[rownames(tmp) %in% genes,] # filter genes 
  aver_exp  <- AverageExpression(tmp2,group.by = "subclass",features = FindVariableFeatures(tmp2,nfeatures = 3000) %>% VariableFeatures()) 
  tree <- aver_exp[[1]] %>% pvclust::pvclust() 
  
  plot3 <- ggtree(tree) + geom_tiplab()
  
  cl.mat <- tmp2@assays$RNA@data[FindVariableFeatures(tmp2,nfeatures = 3000) %>% VariableFeatures(),]
  
  cl.med <- do.call("cbind", tapply(names(tmp2$subclass), tmp2$subclass %>% setNames(NULL), function(x) {
    matrixStats::rowMedians(as.matrix(cl.mat[, x]))
  }))
  
  cl.tree <- cl.med %>% pvclust::pvclust(parallel = 5L) 
  
  plot4 <- ggtree(cl.tree) + geom_tiplab()
  

  ggsave(filename = str_c("2_QC_",tmp$orig.ident[[1]],".png"),  wrap_plots(list(plot1,plot2, plot3, plot4),nrow = 1),width = 18, height = 8,limitsize = F)
  
  })
