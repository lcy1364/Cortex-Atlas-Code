library(qs)
library(tidyverse)
library(Seurat)
library(magrittr)
library(treeio)
library(ggtree)
library(ape)
library(scRNAtoolVis)
setwd("~/cortex/fig4/")

merge_seu <- qs::qread("merge_seu.qs",nthreads = 10)
geneId_Name <-read_csv("../SnRNA/1_SnRNA_preprocessing/gene_kept.csv") %>% {setNames(object = .$gene_uni,nm = .$gene_id)}

names(geneId_Name) %<>% str_extract("ENSG[0-9]*")

sst <- merge_seu[,str_detect(merge_seu$subclass,"SST")]
sst %<>% DietSeurat(counts = F) %>% ScaleData() %>% FindVariableFeatures() %>% RunPCA()  %>% RunHarmony(group.by.vars = "batch",ncores = 10) 
  
sst %<>%  RunUMAP(reduction = "harmony",reduction.name = "umap",reduction.key = "harmonyumap_",dims = 1:20) 

spatialCellMeta <- read_csv("../STEREO/spatialCellMeta.csv")

cluster_depth <- spatialCellMeta %>% filter(subclass == "SST") %>% group_by(cluster) %>% summarise(depth = mean(depth)) %>% {setNames(.$depth,nm = .$cluster)}

sst@meta.data$depth <- cluster_depth[sst@meta.data$cluster]
fig4c_1 <-  featureCornerAxes(object = sst,reduction = 'umap',
                              groupFacet = "NULL",
                              relLength = 0.5,
                              relDist = 0.2,
                              features = "depth")

fig4c_1 <- FeaturePlot(sst,features = "depth") + scale_color_gradientn(colours = MetBrewer::MetPalettes$Greek[[1]] %>% rev,breaks = 0:7,labels = c("ARACHNOID", "L1", "L2", "L3", "L4", "L5", "L6", "WM"),limits = c(1,6)) + cowplot::theme_map() + coord_fixed()
ggsave("fig4c_1.pdf",fig4c_1,width = 5,height = 5)
fig4c_2 <- DimPlot(object = sst,reduction = 'umap',group.by = "cluster",label = T, size = 0.1,label.size = 2,repel = T) + coord_fixed() + scale_color_manual(values = MetBrewer::met.brewer("Benedictus",n = 38)) + cowplot::theme_map()
fig4c_2$layers <- fig4c_2$layers[2]
ggsave("fig4c_2.pdf",fig4c_2,width = 5,height = 5)

sst$depth_cluster = ifelse(sst$cluster %in% c("SST_16", "SST_36", "SST_4", "SST_23", "SST_3", "SST_13", "SST_28", 
                                              "SST_37", "SST_11", "SST_6", "SST_8", "SST_20", "SST_29", "SST_18", 
                                              "SST_22", "SST_27", "SST_15", "SST_9", "SST_1"),yes = "deep",no = "shallow")

qsave(sst,"sstRNA.qs")

