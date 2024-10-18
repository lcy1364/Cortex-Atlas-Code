library(qs)
library(tidyverse)
library(magrittr)
setwd("~/cortex/STEREO/3_domainClustering/")

Slides <- read.delim("../cortex")


annoFiles <- list.files("~/cortex/STEREO/3_domainClustering/","*r.csv",full.names = T) # %>% str_subset(selectedSlides$chip %>% str_c(collapse = "|")) 

bin200SeuFiles <- list.files(".","200.qs",full.names = T) 

# check whether two files can map or not 

problems_file <- list()
i = 1
for (i in 1:nrow(selectedSlides)) {
  chip <- selectedSlides$chip[[i]]
  region <-  selectedSlides$region[[i]]
  # print(chip)
  # timestamp()
  seu <- bin200SeuFiles %>% str_subset(chip) %>% qread()
  
  seu@meta.data$dissection <- region
  seu@meta.data %<>% rownames_to_column("cellID") %>% unite(col = "locID",c(x,y),remove = F) %>% column_to_rownames("cellID")
  
  anno <- annoFiles %>% str_subset(chip) %>% read.csv(row.names = 1) %>% unite(col = "locID",c(x,y),remove = F) 
  seu@meta.data %<>% left_join(anno %>% select(locID,region,mainLayer),by = c("locID" = "locID")) %>% column_to_rownames("cellId")
  
  notMatched <- sum(!(anno$locID %in% seu$locID))
  if(notMatched) {
    print(str_c("i=",i," ",chip," ",notMatched,"/", nrow(anno), " ", round(100*notMatched/nrow(anno),2),"%" ))
  } 
  qsave(seu,bin200SeuFiles %>% str_subset(chip) %>% basename(),nthreads = 10)
}
qsFiles =  list.files(".","200.qs",full.names = T) 


datalist <- parallel::mclapply(qsFiles, mc.cores = 10, function(qsFile) {
  # for (qsFile in qsFiles) {
  print(timestamp())
  print(qsFile)
  
  tmp <- qread(qsFile,nthreads = 20)
  chipID  <- qsFile %>% basename %>% str_extract("^[A-Z0-9]*")
  tmp$chipID <- chipID
  
  tmp %<>% PercentageFeatureSet(pattern = "^MT-",col.name = "mitoRatio") %>%
    NormalizeData() %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData() # %>%
  #SCTransform(vars.to.regress = c("mitoRatio"),vst.flavor = "v2",return.only.var.genes = F)
  
  
  tmp %T>% qsave(qsFile,nthreads = 20)
})

# datalist <- parallel::mclapply(qsFiles, mc.cores = length(qsFiles),qread) 


seurat_merge <- merge(datalist[[1]],datalist[-c(1)],add.cell.ids = str_extract(basename(qsFiles),"[A-Z0-9]*")) 
VariableFeatures(seurat_merge) <- SelectIntegrationFeatures(object.list = datalist,nfeatures = 5000)

seurat_merge <- seurat_merge[,!is.na(seurat_merge$mainLayer)]
seurat_merge %>% qsave("Seu_merge_slides.qs",nthreads = 50)

