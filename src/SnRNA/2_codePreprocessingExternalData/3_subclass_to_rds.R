library(parallel)
library(magrittr)
library(tidyverse)
library(rtracklayer)
library(anndataR)
library(Seurat)
setwd("~/cortex/SnRNA/2_codePreprocessingExternalData")

gtf <- readRDS("../1_SnRNA_preprocessing/geneSym_to_geneID.RDS")

files <- list.files(".","_subclass.h5ad")  

subclasses <- files %>% str_remove("_subclass.h5ad") %>% str_to_upper()
subclass = "ASTRO"

# trans back to counts

subclass_data <- mclapply(subclasses,mc.cores = length(unique(subclasses)),function(subclass){
  
  adata <- files[subclasses == subclass] %>% read_h5ad(to = "InMemoryAnnData")  
  
  data <- adata$X
  
  data@Dimnames[[1]] <- adata$obs_names
  data@Dimnames[[2]] <- adata$var_names
  
  counts <- data

  counts@x <- (exp(counts@x) - 1)/1000000

  nCount <-  apply(counts,MARGIN = 1, FUN = function(x) {
    (1 / (x[x>0] %>% min) ) %>% round() })

  counts <- counts * nCount
  counts@x %>% table
  counts@x %<>% round

  # tmp <- LogNormalize(counts,scale.factor = 10^6) 
  
  seu <- CreateSeuratObject(counts = t(counts),meta.data = adata$obs)
  
  seu <- NormalizeData(seu,scale.factor = 10^6) # logCPM
  
  colnames(seu@meta.data) %<>% str_to_lower()
  seu$library_prep <- str_remove(colnames(seu),  pattern = "[ATGC]+-[0-9]{0,3}") %>% str_remove("-[0-9]*$")
  seu$subclass <- seu$crossarea_subclass %>% str_to_upper()
  seu$cluster <- seu$crossarea_cluster %>% str_to_upper()
  seu$donor <- seu$donor_id
  # seu <- seu[rownames(seu) %in% gtf,]
  
  
  seu$subclass[str_detect(seu$subclass,regex("MICRO",ignore_case = T))] = "MICRO"
  seu$subclass[str_detect(seu$subclass,regex("SNCG",ignore_case = T))] = "NDNF"
  seu$subclass[str_detect(seu$subclass,regex("LAMP5",ignore_case = T))] = "LAMP5" # -1
  seu$subclass[str_detect(seu$subclass,regex("SST",ignore_case = T))] = "SST" # -1
  
  seu$subclass[str_detect(seu$subclass,regex("L2/3 IT",ignore_case = T))] = "L2-L3 IT LINC00507" #  
  seu$subclass[str_detect(seu$subclass,regex("L4 IT",ignore_case = T))] = "L3-L4 IT RORB"#  
  seu$subclass[str_detect(seu$subclass,regex("L5 IT",ignore_case = T))] = "L4-L5 IT RORB"
  seu$subclass[str_detect(seu$subclass,regex("ET",ignore_case = T))] =  "ET"
  seu$subclass[str_detect(seu$subclass,regex("NP",ignore_case = T))] =  "NP"
  seu$subclass[str_detect(seu$subclass,regex("CAR3",ignore_case = T))] =  "L6 CAR3"
  seu$subclass[str_detect(seu$subclass,regex("L6b",ignore_case = T))] =  "L6b"
  seu$subclass[str_detect(seu$subclass,regex("ASTRO",ignore_case = T))] =  "AST"
  
  
  seu$cluster[str_detect(seu$cluster,regex("SST CHODL_1",ignore_case = T))] = "SST_38" # -1
  seu$cluster[str_detect(seu$cluster,"LHX6")] = "LAMP5_8" # -1
  
  seu$cluster <- str_c(seu$subclass,str_extract(seu$cluster,"_[0-9]*$"))  #%>% unique
  seu
})


# merge some of them


subclasses[str_detect(subclasses,regex("LAMP5",ignore_case = T))] = "LAMP5" # -1
subclasses[str_detect(subclasses,regex("SST",ignore_case = T))] = "SST" # -1

sampled_data <- mclapply(unique(subclasses),mc.cores = length(unique(subclasses)),function(subclass){
  
  seu <- subclass_data[subclasses == subclass] %>% purrr::reduce(merge)

  qs::qsave(seu,str_c(str_replace(seu$subclass[[1]],"/","_"),"_subclass.qs"),nthreads = 10,preset = "fast")
  
  # sample 10% for annotating datasets globally
  
  seu[,sample(colnames(seu),ceiling(length(colnames(seu))*0.2),replace = F)]
})  


seu_edlein <- merge(sampled_data[[1]],sampled_data[-c(1)])

seu_edlein$batch = "edlein"

seu_edlein %>% qs::qsave("SnRNA_0.2DownSampled.qs",nthreads = 20,preset = "fast")


