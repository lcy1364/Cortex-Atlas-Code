library(Seurat)
library(tidyverse)
set("~/cortex/figS1-6/")

spatialMeta <- read.csv("spatialMeta.csv")

{ggplot(spatialMeta,aes(x,y,color = subclass)) + geom_point() + facet_grid(~subclass)} %>% ggsave(filename = "figS1a.pdf")
