library(magrittr)
library(qs)
library(Seurat)
library(stringr)
library(tidyverse)
setwd("~/cortex/figS1-6/")
merge_seu <- readRDS("../SnRNA/SnRNA_seurat.RDS")

opcRNA <- merge_seu[,merge_seu$subclass == "OPC"]
geneName_id <- read.csv("../SnRNA/1_SnRNA_preprocessing/gene_kept.csv") %>% {setNames(.$gene_id,.$gene_name)}
setwd("~/cortex/figS1-6/")


PlotData = DotPlot(opcRNA,features = geneName_id[c("FOXP2","THEMIS","PCDH15")],group.by = "region") %$% data

{ggplot(PlotData,aes(x = feature.groups, y = id, size =  avg.exp.scaled,  fill =  avg.exp.scaled)) + 
    geom_point(shape = 22)  + 
    scale_size_continuous(name = "average expression",range = c(0,15)) + 
    MetBrewer::scale_fill_met_c(name = "Benedictus",direction = -1) + 
    cowplot::theme_minimal_grid()  + 
    guides(fill = guide_legend("expression"), size = guide_legend("expression"))  + 
    theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))} %>% 
  ggsave(filename = "figS4f.pdf",width = 5, height = 8)
