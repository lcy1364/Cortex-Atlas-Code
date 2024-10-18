library(Seurat)
library(data.table)
library(ggplot2)
library(plotly)
library(ggforce)
library(rtracklayer)
library(ggtree)
library(tidyverse)
library(magrittr)
library(patchwork)
library(qs)
library(parallel)
library(ggtern)
library(sf)
setwd("~/data/STEREO/AnalysisPlot/")

chipID = chipIDs[[1]]

gtf <- import("../Homo_sapiens.GRCh38.93.gtf") %>% as.data.frame() 

gtf %<>% dplyr::select(seqnames,gene_name,gene_name) %>% distinct()
geneID_name <- gtf$gene_name %>% setNames(gtf$gene_id)

outputDir <- "figS4e/"

dir.create(outputDir)

class_markers <- c(
  "FOXP2", "THEMIS", "PCDH15"
) 

 shapeQsMain <- list.files("./1_mainLayer/","_shape.qs",full.names = T)

seu <- qread("merge_slides.qs") %>% SplitObject(split.by = "chipID")

 
mclapply(chipIDs, mc.cores = length(chipIDs), function(chipID){
   
   print(chipID)
   
   binSize <- 200
  
  expression <-  seu$chipID

  points = expression@meta.data[c("x","y")] %>% as.matrix %>% st_multipoint() %>% st_geometry() %>% st_cast("POINT")

  expression %<>% NormalizeData(features = rownames(expression))
 
  t2 <- list()
  
  for (markerG in class_markers) {
    
    markerData <- FetchData(expression, c("x","y",markerG))

    print(markerG)
    t2[[markerG]] <- ggplot(markerData %>% filter(get(markerG) > 0) %>%  top_frac(n = 1) %>% mutate(exp = scales::oob_squish(scale(get(markerG)), range(-3,3)))) + 
      geom_point(size = 0.5, mapping = aes_string("x","y",color = "exp", alpha = "exp")) + 
      scale_color_gradientn(colors = RColorBrewer::brewer.pal(11,"Spectral") %>% rev, limits = c(-3,3)) + 
      geom_sf(data = shapeSF, linewidth = 0.1, mapping = aes(geometry = geometry ), fill = NA, color = shapeSF$color) + cowplot::theme_map() + labs(title = markerG) + theme(plot.title = element_text(size = 20))
    
  } 
 
  ggsave(file.path(outputDir,str_c("figS4e_",chipID,"_","_markers.pdf")),patchwork::wrap_plots(t2,nrow = 1),height = 5,width = 16, dpi = 100,limitsize = FALSE)

  
} )
