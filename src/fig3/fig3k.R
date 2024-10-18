library(Seurat)
library(tidyverse)
spatialMeta <- read.csv("spatialMeta.meta")

#fig3k
{ggplot(spatialMeta %>% filter(str_detect(subclass,"ET"),region %in% c("S1","S1E","PoCG","ACC","ITG","AG","FPPFC")),aes(x,y,color = subclass)) + geom_density()  + geom_point() + facet_grid(~subclass)} %>% ggsave(filename = "fig3k.pdf")

