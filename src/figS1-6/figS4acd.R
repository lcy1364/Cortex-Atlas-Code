library(tidyverse)
library(data.table)
library(magrittr)
library(Seurat)
library(harmony)
setwd("~/cortex/figS1-6/")

merge_seu <- readRDS("../SnRNA/SnRNA_seurat.RDS")
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

NON <- merge_seu[,merge_seu$class == "non-neuronal"]

NON %<>% ScaleData() %>% FindVariableFeatures() %>% RunPCA() %>% RunHarmony(group.by.vars = "batch",ncores = 10) %>% RunUMAP(reduction = "harmony",reduction.name = "harmony_umap",reduction.key = "harmonyumap_",dims = 1:40) 

NON$subclass %<>% str_to_upper()

# figS4a 
P_labeled <- DimPlot(NON,group.by = "subclass",reduction = "harmony_umap",label.box = T,repel = T,label.size = 10,label = T,pt.size = 0.1) + scale_fill_manual(values = subclass_color,na.value = "#F14315") + scale_color_manual(values = subclass_color,na.value = "#F14315")  + cowplot::theme_map() + theme(legend.position = "none",title = element_text()) + coord_fixed()
ggsave("figS4a_NON_harmonyed_withLabel.pdf",P_labeled,height = 10, width = 10)
ggsave("figS4a_NON_harmonyed_withLabel.png",P_labeled,height = 10, width = 10)

# figS4c

Meta <- read.csv("../SnRNA/3_mergingDatasets/SnRNA_Meta.csv")
{ggplot(Meta %>% filter(str_detect(class,"non")), aes(x = region, fill = subclass)) + geom_bar(position = "fill") + scale_fill_manual(values = subclass_color) + cowplot::theme_cowplot() + theme(axis.text.x = element_text(angle = 45,hjust = 1, vjust = 1)) } %>% 
  ggsave(filename = "figS4c.pdf",height = 5,width = 6)


# figS4d 
P_rg_labeled <- DimPlot(NON,group.by = "region",reduction = "harmony_umap", repel = T,pt.size = 0.1) + scale_fill_manual(values = region_color,na.value = "#F14315") + scale_color_manual(values = region_color,na.value = "#F14315")  + cowplot::theme_map() + theme(legend.position = "none",title = element_text()) + coord_fixed()
ggsave("figS4d_NON_harmonyed_withLabel_region.pdf",P_rg_labeled,height = 10, width = 10)
ggsave("figS4d_NON_harmonyed_withLabel.png",P_rg_labeled,height = 10, width = 10)
