library(tidyverse)
setwd("~/cortex/SnRNA/3_mergingDatasets/")
library(data.table)
library(magrittr)
resFiles <- list.files("./","merge.qs")

x = resFiles[[1]]

Meta <- parallel::mclapply(resFiles,mc.cores = length(resFiles),FUN = function(x){
  seu <- qs::qread(x)
  seu@meta.data$sex[str_detect(seu@meta.data$donor,"^S[0-9]{4}")] <- "male" 
  seu@meta.data$sex[seu@meta.data$donor == "S0406"] <- "female" 
  seu@meta.data %<>% filter((source != "lein_10x_layer5_only"| is.na(source)) & region != "MTG") %>% select(where(~ !any(is.na(.))))
  
  seu@meta.data$library_prep[seu@meta.data$batch == "us"] <- seu@meta.data$orig.ident[seu@meta.data$batch == "us"]
  seu@meta.data
})


sharedCols <- Meta %>% purrr::map(colnames) %>% purrr::reduce(intersect)
Meta %<>% purrr::reduce(~rbind(.x[sharedCols],.y[sharedCols]))



Meta[Meta$region == "A1","region"] <- "STG"
Meta[Meta$region == "ANG","region"] <- "AG"

Meta$subclass %<>% str_to_upper()

Meta[Meta$subclass %in% c("ET", "L2-L3 IT LINC00507", "L3-L4 IT RORB", "L4-L5 IT RORB", 
                          "L6 CAR3", "L6 CT", "L6 IT", "L6b","L6B", "NP"),"class"] <- "excitatory"

Meta[Meta$subclass %in% c("AST", "ENDO", "MICRO", "OLIGO", "OPC", "VLMC"),"class"] <- "non-neuronal"

Meta[Meta$subclass %in% c("LAMP5", "NDNF", "PAX6", "PVALB", "SST", "VIP","CHANDELIER"),"class"] <- "inhibitory"



Meta$subclass[str_detect(Meta$subclass,"MICRO")] = "MICRO"
Meta$subclass[str_detect(Meta$subclass,"SNCG")] = "NDNF"
Meta$subclass[str_detect(Meta$subclass,"LAMP5")] = "LAMP5" # -1
Meta$subclass[str_detect(Meta$subclass,"SST")] = "SST" # -1

Meta$subclass[str_detect(Meta$cluster,"CHANDELIER")] = "CHANDELIER" # 

Meta$cluster %<>% str_to_upper()
Meta$cluster[str_detect(Meta$cluster,"SST CHODL_1")] = "SST_38" # -1

Meta %>% write.csv("SnRNA_Meta.csv")
