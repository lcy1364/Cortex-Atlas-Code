library(qs)
library(parallel)
library(magrittr)
library(tidyverse)
library(Seurat)
library(org.Hs.eg.db)
library(rrvgo)
setwd("~/data/STEREO/AnalysisPlot/")
seu <- qread("./Seu_merge.qs")

genes <- read.csv("../../ensemble93gtf_rmXY.csv")
seu <- seu[rownames(seu) %in% genes$gene_name,]

colorPallete <- c(ggsci::pal_aaas()(10),
                  ggsci::pal_jama()(10),
                  ggsci::pal_npg()(10),
                  ggsci::pal_lancet()(10),
                  ggsci::pal_frontiers()(10),
                  ggsci::pal_nejm()(10)) %>% unique %>% str_subset(negate = T, "1B1919")



# fig6i -------
cl = "L4a"
# for (cl in markGenes$domain %>% unique) {
mg <- c("EGR1","ENC1","HOPX")

tmp <- SpatialFeaturePlot(seu_merge,images = c("image_B01012C5"),stroke = 0,slot = "scale.data",pt.size.factor = 3,features = mg,alpha = c(0,1),combine = F)  %>% purrr::map(~. + coord_flip() + scale_y_reverse() )
tmp <- patchwork::wrap_plots(tmp, nrow = 1,byrow = T) 
ggsave(str_c("L4abc_",cl,"_image_B01012C5_markers.png"),tmp, height = 21/8,width =  92/8/2,dpi = 150,limitsize = FALSE )

tmp <- SpatialFeaturePlot(seu_merge,images = c("image_B02222E1"),stroke = 0,slot = "scale.data",pt.size.factor = 3,features = mg,alpha = c(0,1),combine = F)  %>% purrr::map(~. + coord_flip() + scale_y_reverse() )
tmp <- patchwork::wrap_plots(tmp, nrow = 1,byrow = T) 
ggsave(str_c("L4abc_",cl,"_image_B02222E1_markers.png"),tmp, height = 21/8,width =  92/8/2,dpi = 150,limitsize = FALSE )

cl = "L4b"
# for (cl in markGenes$domain %>% unique) {
mg <- c("NEFH","SYT2","SMYD2")

tmp <- SpatialFeaturePlot(seu_merge,images = c("image_B01012C5"),stroke = 0,slot = "scale.data",pt.size.factor = 3,features = mg,alpha = c(0,1),combine = F)  %>% purrr::map(~. + coord_flip() + scale_y_reverse() )
tmp <- patchwork::wrap_plots(tmp, nrow = 1,byrow = T) 
ggsave(str_c("L4abc_",cl,"_image_B01012C5_markers.png"),tmp, height = 21/8,width =  92/8/2,dpi = 150,limitsize = FALSE )

tmp <- SpatialFeaturePlot(seu_merge,images = c("image_B02222E1"),stroke = 0,slot = "scale.data",pt.size.factor = 3,features = mg,alpha = c(0,1),combine = F)  %>% purrr::map(~. + coord_flip() + scale_y_reverse() )
tmp <- patchwork::wrap_plots(tmp, nrow = 1,byrow = T) 
ggsave(str_c("L4abc_",cl,"_image_B02222E1_markers.png"),tmp, height = 21/8,width =  92/8/2,dpi = 150,limitsize = FALSE )

cl = "L4c"
 
mg <- c("EGR1","ENC1","HOPX","NEFH","SYT2","CD74","RAB3B","CNTN5")


tmp <- SpatialFeaturePlot(seu_merge,images = c("image_B01012C5"),stroke = 0,slot = "scale.data",pt.size.factor = 3,features = mg,alpha = c(0,1),combine = F)  %>% purrr::map(~. + coord_flip() + scale_y_reverse() )
tmp <- patchwork::wrap_plots(tmp, nrow = 1,byrow = T) 
ggsave(str_c("L4abc_",cl,"_image_B01012C5_markers.png"),tmp, height = 21/8,width =  92/8/2,dpi = 150,limitsize = FALSE )

tmp <- SpatialFeaturePlot(seu_merge,images = c("image_B02222E1"),stroke = 0,slot = "scale.data",pt.size.factor = 3,features = mg,alpha = c(0,1),combine = F)  %>% purrr::map(~. + coord_flip() + scale_y_reverse() )
tmp <- patchwork::wrap_plots(tmp, nrow = 1,byrow = T) 
ggsave(str_c("L4abc_",cl,"_image_B02222E1_markers.png"),tmp, height = 21/8,width =  92/8/2,dpi = 150,limitsize = FALSE )

 

 

