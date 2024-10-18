library(tidyverse)
setwd("~/cortex/fig1/")
library(data.table)
library(magrittr)
library(Seurat)
library(harmony)

merge_seu <- readRDS("../SnRNA/SnRNA_seurat.RDS")

subclass_color <-c(AST = "#665C47", ENDO = "#604B47", ET = "#CEC823", CHANDELIER = "#E25691", 
                   `L2-L3 IT LINC00507` = "#07D8D8", `L3-L4 IT RORB` = "#09B2B2", 
                   `L4-L5 IT RORB` = "#69B199", `L6 CAR3` = "#898325", `L6 CT` = "#2D8CB8", 
                   `L6 IT` = "#19E3BE", L6B = "#7944AA", L6b = "#7944AA",LAMP5 = "#D96C07", MICRO = "#513577", 
                   NDNF = "#f2798d", NP = "#93C43B", OLIGO = "#5E8A79", OPC = "#2E3E39", 
                   PAX6 = "#E96DF2", PVALB = "#FF2D4E", SST = "#F2B003", VIP = "#9C28CC", 
                   VLMC = "#697255")

P_global_label <- DimPlot(merge_seu,group.by = "subclass",reduction = "harmony_umap",repel = T,split.by = "batch",label.box = T,label = T,pt.size = 2) + scale_fill_manual(values = subclass_color) + scale_color_manual(values = subclass_color)  + cowplot::theme_map() + theme(legend.position = "none",title = element_text()) + coord_fixed()
ggsave("1_b_merge_harmonyed_batch.pdf",P_global_label,height = 5, width = 10)
ggsave("1_b_merge_harmonyed_batch.png",P_global_label,height = 5, width = 10)

P_global_unlabel <- DimPlot(merge_seu,group.by = "subclass",reduction = "harmony_umap",split.by = "batch",pt.size = 2) + scale_fill_manual(values = subclass_color) + scale_color_manual(values = subclass_color)  + cowplot::theme_map() + theme(legend.position = "none",title = element_text()) + coord_fixed()
ggsave("fig1b_merge_harmonyed_batch_withoutLabel.pdf",P_global_unlabel,height = 5, width = 10)
ggsave("fig1b_merge_harmonyed_batch_withoutLabel.png",P_global_unlabel,height = 5, width = 10)

P_global_whole <- DimPlot(merge_seu,group.by = "subclass",reduction = "harmony_umap",label.box = T,repel = T,label = T,pt.size = 2) + scale_fill_manual(values = subclass_color) + scale_color_manual(values = subclass_color)  + cowplot::theme_map() + theme(legend.position = "none",title = element_text()) + coord_fixed()

ggsave("fig1b_merge_harmonyed.pdf",P_global_whole,height = 5, width = 5)
ggsave("fig1b_merge_harmonyed.png",P_global_whole,height = 5, width = 5)

P_global_whole_unlabel <- DimPlot(merge_seu,group.by = "subclass",reduction = "harmony_umap",pt.size = 2) + scale_fill_manual(values = subclass_color) + scale_color_manual(values = subclass_color)  + cowplot::theme_map() + theme(legend.position = "none",title = element_text()) + coord_fixed()

ggsave("fig1b_merge_harmonyed_withoutLabel.pdf",P_global_whole_unlabel,height = 5, width = 5)
ggsave("fig1b_merge_harmonyed_withoutLabel.png",P_global_whole_unlabel,height = 5, width = 5)


