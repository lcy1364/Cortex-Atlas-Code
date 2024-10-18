# fig6d -----------
library(tidyverse)
library(Seurat)
library(magrittr)

seu <- qs::qread("merge_seu.qs")

SpatialDimPlot(seu[,seu$region == "V1"],group.by = "domain.fine") %>% ggsave(filename = "fig6d.pdf")

# for fig6e refer to the previous cellbin plot in fig1b please
