library(Seurat)
library(magrittr)
library(tidyverse)
library(tcltk)
library(harmony)
library(rtracklayer)
library(harmony)
library(MASS)
library(patchwork)
library(FNN)
setwd("~/cortex/SnRNA/3_mergingDatasets/")
devtools::load_all("~/ClusterGVis-main/")
setwd("~/cortex/SnRNA/3_mergingDatasets/")

# initial assignment ------

# external
seu_edlein <- qs::qread("~/cortex/SnRNA/2_codePreprocessingExternalData/SnRNA_0.2DownSampled.qs")  # counts and data are both lognormalization
seu_edlein$batch = "edlein"
# gtf <- import("../1_SnRNA_preprocessing/gencode.v32.primary_assembly.annotation.gtf") %>% as.data.frame() %>% mutate(gene_id = str_remove(gene_id,"\\..*")) #%>% select(seqnames,gene_name,gene_id) %>% mutate(gene_id = str_remove(gene_id,"\\..*")) %>% distinct() %>% {setNames(.$gene_id,nm = .$gene_name)}

# ours
seu_us <- qs::qread("~/cortex/SnRNA/1_SnRNA_preprocessing/SnRNA_seurat.qs")

seu_us$library_prep = seu_us$orig.ident
seu_us$batch = "us"
seu_us$subclass_old <- seu_us$subclass

sharedGenes <- intersect(rownames(seu_us),rownames(seu_edlein))

timestamp()
seu <- merge(seu_us[sharedGenes,],seu_edlein[sharedGenes,]) 

seu %<>% NormalizeData(scale.factor = 10^6) %>% FindVariableFeatures() %>% ScaleData(split.by = "batch") %>% RunPCA() 

seu %<>% RunHarmony(group.by.vars = "batch",lambda = 2) %>% RunUMAP(reduction = "harmony",dims = 1:30,reduction.name = "umap_har")
timestamp()

{DimPlot(seu,reduction = "umap_har",split.by = "batch",group.by = "subclass",shuffle = T,label = T,label.box = T) + coord_fixed() + theme(legend.position = "none")} %>% ggsave(filename = "1_subclass_mapping.png",width = 20,height = 13,limitsize = F)

seu %>% qs::qsave("1_merged_seu.qs",preset= "fast")

# assign cells 

cellID <- FetchData(seu,c("UMAP_1","UMAP_2","batch","subclass","class"))

library(RANN)

edlein <- cellID[cellID$batch=="edlein",]
us <-  cellID[cellID$batch=="us",]

library(patchwork)
library(ggrastr)


res <- RANN::nn2(edlein[c("UMAP_1","UMAP_2")],query = us[c("UMAP_1","UMAP_2")],k = 20)
x = 1

us$subclass <- parallel::mclapply(1:nrow(res$nn.idx),mc.cores = 50,FUN = function(x){
  edlein[res$nn.idx[x,],"subclass"]  %>% table %>% which.max %>% names 
  }) %>% unlist


{ggplot(cellID,mapping = aes(UMAP_1,UMAP_2,color = subclass)) + geom_point_rast() + facet_wrap(~batch) + coord_fixed() + theme(legend.position = "none")} %>% ggsave(filename = "1_subclass_mapping.png",width = 20,height = 13,limitsize = F)

# {ggplot(cellID,mapping = aes(UMAP_1,UMAP_2,color = subclass)) + geom_point_rast()} %>% ggsave("1_subclass_mapping.png",width = 20,height = 10,limi)


seu_us$subclass <- us[colnames(seu_us),"subclass"]
seu_us$class <- us[colnames(seu_us),"class"]

# fiter by nFeatures based on major types

seu_us <- seu_us[, seu_us$nFeature_RNA >=1000|str_detect(seu_us$class,regex("non",ignore_case = T))] # filter by class

seu_us %>% qs::qsave("~/cortex/SnRNA/3_mergingDatasets/SnRNA_seurat.qs",preset = "fast")

# further assignment ------

# define function to find the most frequent knn items around our data
most_frequent_elements <- function(row) {
  freq_table <- table(row)
  most_common_element <- names(freq_table)[which.max(freq_table)]
  return(most_common_element)
}


files <- list.files("../2_codePreprocessingExternalData/","_subclass.qs",full.names = T) 
files %<>% setNames(nm = files %>% basename %>% str_remove_all("_subclass.qs") %>% str_to_upper()) 

seu_us_list <- seu_us %>% SplitObject("subclass") 

names(seu_us_list) %<>% str_to_upper() %>% str_replace("/","_")

setequal(unique(files %>% names) %>% sort,unique(seu_us_list %>% names) %>% sort) #%>% keep(~.!="IT")

subclasses <- intersect(unique(files %>% names),unique(seu_us_list %>% names))

x = subclasses[[1]]
i = 1

seu_list <- parallel::mclapply(1:length(subclasses),mc.cores = length(subclasses),FUN = function(i){
  x = subclasses[[i]]
  edlein <- files[names(files) == x] %>% purrr::map(~qs::qread(.)) %>% purrr::reduce(merge)
  edlein <- edlein[,edlein$region != "MTG"]
  edlein$batch = "edlein"
  
  edlein <- edlein[rowSums(edlein) >0,] #
  
  us <- seu_us_list[names(seu_us_list) == x] %>% purrr::reduce(merge)
  us <- us[rowSums(us) >0,] 
  us %<>% NormalizeData(scale.factor = 10^6)
  us$batch = "us"
  
  sharedGenes <- intersect(rownames(us),rownames(edlein))
  us <- us[sharedGenes,]
  edlein <- edlein[sharedGenes,]
  
  merged_subclass <- merge(us,edlein)
  
  merged_subclass %<>% ScaleData() %>% FindVariableFeatures() %>% RunPCA() %>% RunHarmony(group.by.vars = "batch",dims.use = 1:30,ncores = 5,lambda = 2) %>% 
    RunUMAP(reduction = "harmony",reduction.name = "harmony_umap",reduction.key = "harmonyumap_",dims = 1:30) 
  
  # assign cells 
  
  umapData <- FetchData(merged_subclass,c("harmonyumap_1","harmonyumap_2","batch","cluster")) #%>% group_by(batch) %>% group_split() 
  edlein_umap <- umapData[umapData$batch=="edlein",] 
  us_umap <- umapData[umapData$batch=="us",] 
  res <- RANN::nn2(edlein_umap[c("harmonyumap_1","harmonyumap_2")],query = us_umap[c("harmonyumap_1","harmonyumap_2")],k = 20)
  
  most_frequent <- parallel::mclapply(1:nrow(res$nn.idx), function(i) most_frequent_elements(edlein_umap[res$nn.idx[i, ],"cluster"]), mc.cores = 5)
  
  merged_subclass@meta.data[colnames(us),"cluster"] = unlist(most_frequent)
  
  Idents(merged_subclass) <- "cluster"
  
  DimPlot(merged_subclass,reduction = "harmony_umap",split.by = "batch",group.by = "cluster",shuffle = T,label = T,label.box = T) %>% ggsave(filename = str_c("2_",x,".png"),width = 10,height = 6,limitsize = F)
  
  
  merged_subclass@active.ident <- fct_relevel(merged_subclass@active.ident, sort(levels(merged_subclass@active.ident)))
  
  merged_subclass %>% qs::qsave(str_c("2_",x,"_merge.qs"),nthreads = 10,preset = "fast") 
  
  merged_subclass
})

names(seu_list) <- purrr::map(seu_list,~.$subclass[[1]]) %>% unlist 

x = seu_list[[1]]

parallel::mclapply(seu_list, mc.cores = length(seu_list),FUN = function(x){
  
  edlein <- x[,x$batch == "edlein"]
  
  edlein <- edlein[,edlein@meta.data %>% rownames_to_column("cellID") %>% group_by(cluster) %>%
                     group_modify(~ {
                       if (nrow(.x) < 1000) {
                         return(.x)
                       } else {
                         return(sample_n(.x, 1000))
                       }
                     }) %$% cellID]
  
  
  Idents(edlein) <- edlein$cluster %>% factor()
  edlein@active.ident <- fct_relevel(edlein@active.ident, sort(levels(edlein@active.ident)))
  celltypes <- edlein@active.ident %>% levels
  
  marker <-
    mclapply(celltypes,mc.cores = 5,FUN = function(ct){
      markers <- FindMarkers(edlein,ident.1 = ct,ident.2 = NULL,only.pos = T)
      markers$cluster = ct
      markers$gene = rownames(markers)
      markers
    }) %>% purrr::reduce(rbind)
  
  marker %<>% group_by(cluster) %>% slice_max(order_by = avg_log2FC,n = 200) %>% arrange(cluster,desc(avg_log2FC))
  
  marker$subclass <- x$subclass[[1]]
  x@misc$marker <- marker
  
  Idents(x) <- str_c(x$cluster,x$batch,sep = "_") %>% factor()
  
  plotdata <- prepareDataFromscRNA(object = x,
                                   diffData = marker,
                                   showAverage = TRUE)
  
  edleinClusterNum = edlein$cluster %>% unique %>% length
  pdf(str_c("2_",x$subclass[[1]],"_markers.pdf"),height = edleinClusterNum * 1,width = edleinClusterNum * 0.5 + 2,onefile = F)
  visCluster(object = plotdata,
             plot.type = "heatmap",
             column_names_rot = 45,
             cluster_rows = F,
             cluster.order = 1:length(plotdata$wide.res$cluster %>% unique)
  )
  dev.off()
  
  marker %>% qs::qsave(str_c("2_",x$subclass[[1]],"_marker.qs"),nthreads = 10, preset = "fast")
  
})







