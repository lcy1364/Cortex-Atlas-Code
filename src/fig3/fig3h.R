library(Seurat)
library(tidyverse)
spatialMeta <- read.csv("spatialMeta.meta")
 
{ggplot(spatialMeta %>% filter(str_detect(subclass,"IT"),region == "ACC"),aes(x,y,color = subclass)) + geom_point() + facet_grid(~subclass)} %>% ggsave(filename = "fig3h.pdf")

