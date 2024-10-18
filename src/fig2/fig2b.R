library(Seurat)
library(tidyverse)
library(stringr)
setwd("~/cortex/fig2/")

domainColor <- c(
  ARACHNOID = "#8a3b35",  
  L1 = "#8ba28e",         
  L2 = "#9ec87e",         
  L3 = "#669c68",         
  L4 = "#67b8bb",         
  L5 = "#5687ac",         
  L6 = "#5d5d8d",         
  WM = "#b5b3bb"          
)

gene_list <-
  c(
    "MBP",
    "SNAP25",
    "CTGF",
    "THEMIS",
    "PCP4",
    "RORB",
    "COL5A2",
    "CUX2",
    "LAMP5",
    "PCDH8",
    "LINC00507",
    "RELN"
  )

chipIDs <- cortexMeta %$% ChipID
t2 <- list()
bins <- list.files("../STEREO/GEM/bin200/", "200.qs", full.names = T) # bin200 seurat obejct
mclapply(chipIDs, mc.cores = length(chipIDs), function(chipID) {
  print(chipID)
  expression <- qread(bins %>% keep(~ str_detect(., chipID)))
  for (markerG in gene_list) {
    print(markerG)
    t2[[markerG]] <-
      ggplot(
        FetchData(expression, c("x", "y", markerG)) %>% filter(get(markerG) > 0) %>%  top_frac(n = 1) %>% mutate(exp = scales::oob_squish(scale(get(
          markerG
        )), range(-3, 3)))
      ) + geom_point(size = 0.5,
                     mapping = aes_string("x", "y", color = "exp", alpha = "exp")) +
      scale_color_gradientn(colors = RColorBrewer::brewer.pal(11, "Spectral") %>% rev,
                            limits = c(-3, 3)) + cowplot::theme_nothing() + labs(title = markerG) + theme(plot.title = element_text(size = 20))
    
  }
  
  ggsave(
    file.path(
      outputDir,
      str_c(chipID, "_", "_mainLayer_markers.pdf")
    ),
    patchwork::wrap_plots(t2, nrow = 1),
    height = 4,
    width = 4 * 10,
    dpi = 100,
    limitsize = FALSE
  )
  ggsave(
    file.path(
      outputDir,
      str_c(chipID, "_", "_mainLayer_markers.png")
    ),
    patchwork::wrap_plots(t2, nrow = 1),
    height = 4,
    width = 4 * 10,
    dpi = 100,
    limitsize = FALSE
  )
})
