library(Seurat)
library(tidyverse)

setwd("/media/desk16/luomeng/data/cortex/SnRNA/2_codePreprocessingExternalData")
files <- list.files(".","subclass.qs")

All_gene <- parallel::mclapply(files,mc.cores = length(files),FUN = function(x){
  qs::qread(x) %>% rownames
})



All_gene %>% purrr::reduce(c) %>% unique %>% cat(file = "gene.txt",sep = "\n")
