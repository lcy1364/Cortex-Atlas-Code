library(Seurat)
library(parallel)
setwd("~/cortex/SnRNA/3_mergingDatasets/")
seu_us <- qs::qread("SnRNA_seurat.qs")

# filter by each sample

seu_list <- seu_us %>% SplitObject(split.by = "library_prep")

mclapply(seu_list,mc.cores = length(seu_list),FUN = function(x){
  x %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunUMAP()
  
  
})
