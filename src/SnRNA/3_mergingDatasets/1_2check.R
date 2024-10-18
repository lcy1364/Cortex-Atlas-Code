library(Seurat)
library(tidyverse)
setwd("~/cortex/SnRNA/3_mergingDatasets")
subclass_color <-c(AST = "#665C47", ENDO = "#604B47", ET = "#CEC823", CHANDELIER = "#E25691", 
                   `L2-L3 IT LINC00507` = "#07D8D8", `L3-L4 IT RORB` = "#09B2B2", 
                   `L4-L5 IT RORB` = "#69B199", `L6 CAR3` = "#898325", `L6 CT` = "#2D8CB8", 
                   `L6 IT` = "#19E3BE", L6b = "#7944AA", LAMP5 = "#D96C07", MICRO = "#513577", 
                   NDNF = "#f2798d", NP = "#93C43B", OLIGO = "#5E8A79", OPC = "#2E3E39", 
                   PAX6 = "#E96DF2", PVALB = "#FF2D4E", SST = "#F2B003", VIP = "#9C28CC", 
                   VLMC = "#697255")

seu_us_org <- qs::qread("~/luomeng/DATA/data/SnRNA/3_5_seu.rna.rmXY.qs")
seu_us <- qs::qread("~/cortex/SnRNA/3_mergingDatasets/SnRNA_seurat.qs")
seu_us_2 <- qs::qread("~/cortex/SnRNA/1_SnRNA_preprocessing/SnRNA_seurat.qs")

seu_us_org$new_subclass <- seu_us@meta.data[colnames(seu_us_org),"subclass"]

DimPlot(seu_us_org,group.by = c("subclass","new_subclass"),label = T) + scale_color_manual(values = subclass_color)

{DimPlot(seu_us_org,group.by = c("new_subclass"),label = T) + facet_wrap(~new_subclass) + scale_color_manual(values = subclass_color) + theme(legend.position = "none")} %>% ggsave(filename = "1_new_subclass_mapping.png",width = 15,height = 15,dpi=150)
