library(Seurat)
library(tidyverse)
library(stringr)
setwd("~/cortex/fig4/")


gene_list <-
  c(
    "PDZD2",
    "GNAL",
    "GRIA4",
    "CALB1",
    "DCC",
    "TRHDE",
    "GRIN3A")

chipIDs <- cortexMeta %$% ChipID
t2 <- list()
bins <- list.files("../GEF/bins/", "200.qs", full.names = T) # bin200 seurat obejct
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
      str_c(chipID, ".pdf")
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
      str_c(chipID, ".png")
    ),
    patchwork::wrap_plots(t2, nrow = 1),
    height = 4,
    width = 4 * 10,
    dpi = 100,
    limitsize = FALSE
  )
})
