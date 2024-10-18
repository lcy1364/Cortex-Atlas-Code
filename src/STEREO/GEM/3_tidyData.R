library(tidyverse)
library(qs)
library(parallel)
library(BPCells)
library(SeuratObject)
library(SeuratDisk)
library(magrittr)
library(Matrix)
library(RANN)
devtools::load_all("~/seurat/")
setwd("~/cortex/STEREO/GEM/")

cortexMeta <- read.delim("../cortex") %>% column_to_rownames("chip")

qsFiles <- list.files("bin200/","*qs",full.names = T) %>% keep(~str_detect(.,str_c(rownames(cortexMeta),collapse = "|")))

qsFiles %>% basename %>% str_extract("^[0-9A-Z]*") %>% unique %>% length

geneKept <- read_csv("../ensemble93gtf_rmXY.csv")


#read annotation

annoFiles <- list.files("../3_domainProcess/annotation/","csv",full.names = T)

names(annoFiles) <- str_extract(annoFiles %>% basename,"^[0-9A-Z]*")

# [1] 44
foo = qsFiles[[31]]

seu_list <- mclapply(qsFiles, mc.cores = length(qsFiles),FUN = function(foo){
  
  chipID <-  foo %>% basename() %>%  str_extract("^[0-9A-Z]*")
  
  seu <- foo %>% qread %>% DietSeurat()
  
  seu@images <- list()
  
  seu <- RenameCells(seu, new.names = paste0(chipID,"_", Cells(seu)))
  
  writePath = str_c('./bin200/',chipID)
  
  if(!dir.exists(writePath)){
    write_matrix_dir(mat = GetAssayData(seu,slot = "counts"), dir = writePath)
    counts.mat <- open_matrix_dir(dir = writePath)
    } else {
      counts.mat <- GetAssayData(seu,slot = "counts")
    }
  
  seu <- CreateSeuratObject(counts = counts.mat,assay = "Spatial", meta.data = seu@meta.data[c("x","y")])
  
  cell_coords <- seu@meta.data[c('x', 'y')]
  cell_coords['cells'] <- row.names(cell_coords)
  
  seu@images[[str_c(cortexMeta[chipID,"region"],"_",chipID)]] <- new(Class = 'SlideSeq', coordinates = cell_coords)
  seu@images[[1]]@key <- str_c("image",str_to_lower(chipID),"_")
  seu@images[[1]]@assay <- "Spatial"
  
  seu$chip <- chipID
  seu$donor <- cortexMeta[chipID,"donor"]
  seu$region <- cortexMeta[chipID,"region"]

  anno <- read_csv(annoFiles[chipID])
  
  anno <- if ("region" %in% colnames(anno)) {
    cbind(anno[c("x","y")],anno["region"])
  }else{
    cbind(anno[c("x","y")],anno[last(colnames(anno) %>% keep(~str_starts(.,"domain")))])
  }
  
  colnames(anno) <- c("x","y","domain")
  
  anno %<>% filter(domain != "NVU")
  
  res <- RANN::nn2(anno[c("x","y")],seu@meta.data[c("x","y")],k = 1)
  
  seu@meta.data$domain <- anno$domain[res$nn.idx]

  seu@meta.data$domain %<>% str_to_title
  
  seu <- seu[,seu$domain != "Empty"]
  
  seu %>% qs::qsave(file.path("bin200",str_c(chipID,".qs")))
  
  seu
}) 

merge_seu <- merge(seu_list[[1]],seu_list[c(-1)])

merge_seu %<>% NormalizeData(scale.factor = 10^6) 

merge_seu <- merge_seu[rownames(merge_seu) %in% geneKept$gene_name,]

merge_seu %<>% JoinLayers()
merge_seu %<>% FindVariableFeatures 
merge_seu %<>% ScaleData(verbose = T)

merge_seu@assays[[1]]@layers$scale.data %<>% as.matrix()

merge_seu %<>% RunPCA(verbose = T)

merge_seu %<>% harmony::RunHarmony(group.by.vars = "chip",theta = 20)

merge_seu %<>% RunUMAP(dims = 1:30,verbose = T,reduction = "harmony")

DimPlot(merge_seu,split.by = "region",group.by = "domain",label = T)

merge_seu$majorDomain <- strtrim(merge_seu$domain,width = 2) %>% str_to_upper()

merge_seu$majorDomain[merge_seu$majorDomain == "AR"] <- "ARACHNOID"

merge_seu %>% qs::qsave("~/cortex/STEREO/st_domain_seu_44slides.qs")
