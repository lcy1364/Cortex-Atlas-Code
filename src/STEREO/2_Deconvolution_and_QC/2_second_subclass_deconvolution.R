library(qs)
library(tidyverse)
library(sf)
library(magrittr)
library(Seurat)
devtools::load_all("~/spacexr-master/")

setwd("~/cortex/STEREO/2_Deconvolution_and_QC/")

STEREO <- list.files("./","*.rds",full.names = T)

selectedFiles <- read.delim("./") 

refFiles <- list.files("~/DATA/NeoCortex_EdLein_Sten/EdLein/","qs",full.names = T) 

seu <- readRDS(file)
# ct = "MICRO"
cts = selectedFiles$subclass %>% unique %>% sort

for (ct in cts) {

  
  print(ct)
  datasets <- selectedFiles %>% filter(subclass ==ct) %$% subclass2 %>% c %>% purrr::map(~qs::qread(str_c("/home/luomeng/DATA/NeoCortex_EdLein_Sten/EdLein/",.)))
    
    datasets %<>% purrr::reduce(merge)
    
 
    ctCount <- datasets$CrossArea_cluster %>% table()
    ctKept <- names(ctCount)[ctCount >= 25] # 少于25个细胞的类没办法解
    
    datasets <- datasets[,datasets$CrossArea_cluster %in% ctKept]
    datasets@meta.data$CrossArea_cluster %<>% as.character()
    datasets@meta.data$CrossArea_cluster %<>% str_replace_all("/","_")
    datasets@meta.data$CrossArea_cluster %<>% as.character() %>% factor()
    
    
    # datasets@assays$RNA@counts[]
    reference_tmp <- Reference(counts = datasets@assays$RNA@counts,cell_types = datasets@meta.data$CrossArea_cluster %>% setNames(datasets@meta.data %>% rownames()))
    qsave(reference_tmp,str_c(ct,"_ref.qs"))
 
}
# })

referenceFiles <- list.files(".","_ref.qs",full.names = T)

referenceDat <- purrr::map(referenceFiles,~qread(.))
names(referenceDat) <- basename(referenceFiles) %>% str_remove("_ref.qs")

gtf <- read.csv("../../ensemble93gtf_rmXY.csv")

name_id <- gtf$gene_id %>% setNames(gtf$gene_name)

x = STEREO
parallel::mclapply(X = STEREO1,mc.cores = 44,mc.preschedule = F,FUN = function(x){
  seu <- readRDS(x)
  chipName <- x %>% basename() %>% str_remove_all("_cellbin.rds")

  seu_list <- seu %>% SplitObject(split.by = "cellType")
  ct = "ENDO"
  seu_list <- seu_list[!names(seu_list) %in% c("OPC","MICRO")] 
  for (ct in names(seu_list)) {
    # spatial stuff
    ST <- seu_list[[ct]]
    ST <- ST[rowSums(ST@assays$Spatial@counts) >0,]
    countSP <- ST@assays$Spatial@counts
    countSP@Dimnames[[1]] <- name_id[countSP@Dimnames[[1]]]
    countSP <- countSP[!is.na(rownames(countSP)),]

    puck <- SpatialRNA(ST@meta.data[c("x", "y")],
                       countSP)

    ref_tmp <- referenceDat[[ct]]
    ref_tmp@counts <- ref_tmp@counts[rownames(ref_tmp@counts) %>% keep(~. %in% rownames(countSP)),]
    
    # RCTD
    myRCTD <- create.RCTD(puck, referenceDat[[ct]], max_cores = 5)

    qsave(myRCTD, str_c(chipName,"_",ct,"_RCTD.qs"), nthreads = 1)

  }

})
