library(Seurat)
library(tidyverse)
setwd("~/cortex/fig3/")
spatialMeta <- read.csv("../STEREO/spatialCellMeta.csv")


#fig3d 
{ggplot(spatialMeta %>% filter(str_detect(subclass,"IT")&str_detect(subclass,"CAR3",negate = T)),aes(x,y,color = subclass)) + geom_point() + facet_grid(~subclass+chip)} %>% ggsave(filename = "fig3d.pdf")

