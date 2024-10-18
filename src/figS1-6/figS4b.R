library(Seurat)
library(tidyverse)
spatialMeta <- read.csv("spatialMeta.meta")

{ggplot(spatialMeta %>% filter(chipID == "SS200001075BR" & class == "NON"),aes(x,y,color = subclass)) + geom_point()} %>% ggsave(filename = "figS4b.pdf")

