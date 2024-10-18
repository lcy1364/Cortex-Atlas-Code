library(tidyverse)
setwd("~/cortex/fig4/")
library(data.table)
library(magrittr)
library(Seurat)

subclass_color <-c(AST = "#665C47", ENDO = "#604B47", ET = "#CEC823", CHANDELIER = "#E25691", 
                   `L2-L3 IT LINC00507` = "#07D8D8", `L3-L4 IT RORB` = "#09B2B2", 
                   `L4-L5 IT RORB` = "#69B199", `L6 CAR3` = "#898325", `L6 CT` = "#2D8CB8", 
                   `L6 IT` = "#19E3BE", L6B = "#7944AA", L6b = "#7944AA",LAMP5 = "#D96C07", MICRO = "#513577", 
                   NDNF = "#f2798d", NP = "#93C43B", OLIGO = "#5E8A79", OPC = "#2E3E39", 
                   PAX6 = "#E96DF2", PVALB = "#FF2D4E", SST = "#F2B003", VIP = "#9C28CC", 
                   VLMC = "#697255")

merge_seu <- readRDS("../SnRNA/SnRNA_seurat.RDS")

INH <- merge_seu[,str_detect(merge_seu$class,"inhibitory")]

INH %<>% DietSeurat(counts = F) %>% ScaleData() %>% FindVariableFeatures() %>% RunPCA()  %>% RunHarmony(group.by.vars = "batch",ncores = 10) %>% RunUMAP(reduction = "harmony",reduction.name = "harmony_umap",reduction.key = "harmonyumap_",dims = 1:40) 


P_INT_unlabeled <- DimPlot(INH,group.by = "subclass",,pt.size = 0.1,reduction = "harmony_umap",shuffle = T) + scale_fill_manual(values = subclass_color) + scale_color_manual(values = subclass_color)  + cowplot::theme_map() + theme(legend.position = "none",title = element_text()) + coord_fixed()
ggsave("INH_harmonyed_withoutLabel.pdf",P_INT_unlabeled,height = 8, width = 8)
ggsave("INH_harmonyed_withoutLabel.png",P_INT_unlabeled,height = 8, width = 8)

P_INT_labeled <- DimPlot(INH,group.by = "subclass",label.box = T,repel = T,label = T,pt.size = 0.1,shuffle = T,reduction = "harmony_umap") + scale_fill_manual(values = subclass_color,na.value = "#F14315") + scale_color_manual(values = subclass_color,na.value = "#F14315")  + cowplot::theme_map() + theme(legend.position = "none",title = element_text()) + coord_fixed()
ggsave("INH_harmonyed_withLabel.pdf",P_INT_labeled,height = 8, width = 8)
ggsave("INH_harmonyed_withLabel.png",P_INT_labeled,height = 8, width = 8)

