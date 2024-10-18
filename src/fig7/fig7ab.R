library(qs)
library(magrittr)
library(MASS)
library(tidyverse)
library(sf)
library(Seurat)
library(patchwork)
setwd("/home/luomeng/data/STEREO/AnalysisPlot//")
seu <- qs::qread("merge_slides_seu.qs")


# fig7a ---------

SpatialDimPlot(seu,group.by = "domaine.fine") %>% ggsave("fig7a.pdf",height = 10,width = 20,limitsize = F)

# fig7b  ---------
genesList = c("MBP","THEMIS","CTGF","SEMA3E","MGST1")
chipIDs <- cortexMeta  %$% chip
seu_list <- seu %>% SplitObject(split.by = "chipID")
mclapply(chipIDs, mc.cores = length(chipIDs), function(chipID){
  # for (chipID in chipIDs) {
  print(chipID)
  expression <- qread(str_subset(bins,chipID) %>% str_subset("200.qs"))

  markerGPlot = list()
  for (i in 1:length(genesList)) {
    markerG = genesList[[i]]
    # maererGID = genesList[i,"gene"]
    # type = genesList[i,"group"]
    print(markerG)
    if (!markerG %in% rownames(expression)) {
      print(str_c(markerG," not exist"))
      next
    }
    tmp <- ggplot(FetchData(seu_list$chipID, c("x","y",markerG)) %>% filter(get(markerG) > 0) %>%  top_frac(n = 1) %>% mutate(exp = scales::oob_squish(scale(get(markerG)), range(-3,3)))) + geom_point(size = 0.5, mapping = aes_string("x","y",color = "exp", alpha = "exp")) + 
      scale_color_gradientn(colors = RColorBrewer::brewer.pal(11,"Spectral") %>% rev, limits = c(-3,3)) + 
      cowplot::theme_nothing() + theme(plot.title = element_text(size = 20)) +    theme(
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_rect(fill='transparent')
      )
    ggsave(file.path(outputDir,str_c("./",markerG,"_",chipID,"_",cortex[chipID],".png")), tmp,height = 3,width = 3 , dpi = 100,limitsize = FALSE)
  } 
} )


