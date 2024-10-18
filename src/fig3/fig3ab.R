library(tidyverse)
library(magrittr)
library(Seurat)
library(harmony)
setwd("~/cortex/fig3")

subclass_color <-c(AST = "#665C47", ENDO = "#604B47", ET = "#CEC823", CHANDELIER = "#E25691", 
                   `L2-L3 IT LINC00507` = "#07D8D8", `L3-L4 IT RORB` = "#09B2B2", 
                   `L4-L5 IT RORB` = "#69B199", `L6 CAR3` = "#898325", `L6 CT` = "#2D8CB8", 
                   `L6 IT` = "#19E3BE", L6B = "#7944AA", L6b = "#7944AA",LAMP5 = "#D96C07", MICRO = "#513577", 
                   NDNF = "#f2798d", NP = "#93C43B", OLIGO = "#5E8A79", OPC = "#2E3E39", 
                   PAX6 = "#E96DF2", PVALB = "#FF2D4E", SST = "#F2B003", VIP = "#9C28CC", 
                   VLMC = "#697255")

region_color <- c(FPPFC = "#3F4587", DLPFC = "#8562AA", VLPFC = "#EC8561", M1 = "#B97CB5", 
                  S1 = "#D43046", S1E = "#F0592B", PoCG = "#ED4A96", SPL = "#593C97", 
                  SMG = "#A54486", AG = "#FBDE13", V1 = "#299FAE", ITG = "#75CCE3", 
                  STG = "#0C6939", ACC = "#0D9547")


merge_seu <- readRDS("../SnRNA/SnRNA_seurat.RDS")

EXC <- merge_seu[,str_detect(merge_seu$class, "exc")]

EXC %<>% FindVariableFeatures() %>% ScaleData()  %>% RunPCA() 
EXC %<>%  RunUMAP(dims = 1:20)

# EXC$subclass %<>% str_to_upper()

P_EXC_labeled <- DimPlot(EXC,group.by = "subclass",reduction = "harmony_umap",label.box = T,repel = T,label = T,pt.size = 0.5) + scale_fill_manual(values = subclass_color,na.value = "#F14315") + scale_color_manual(values = subclass_color,na.value = "#F14315")  + cowplot::theme_map() + theme(legend.position = "none",title = element_text()) + coord_fixed()
ggsave("fig3a_EXC_harmonyed_withLabel.pdf",P_EXC_labeled,height = 10, width = 10,dpi = 500)

P_EXC_region_labeled <- DimPlot(EXC,group.by = "region",reduction = "harmony_umap",shuffle = T, pt.size = 1) + scale_fill_manual(values = region_color,na.value = "#F14315") + scale_color_manual(values = region_color,na.value = "#F14315")  + cowplot::theme_map() + theme(legend.position = "none",title = element_text()) + coord_fixed()
ggsave("fig3b_EXC_harmonyed_withLabel_region.pdf",P_EXC_region_labeled,height = 10, width = 10,dpi = 500)
