library(Seurat)
library(tidyverse)
library(magrittr)
setwd("~/cortex/fig3/")
spatialMeta <- read.csv("../STEREO/spatialCellMeta.csv")

#fig3c left 
{ggplot(spatialMeta %>% filter(chip == "SS200001075BR"),aes(x,y,color = subclass)) + geom_point()} %>% ggsave(filename = "fig3c.pdf")

#fig3c right
{ggplot(spatialMeta %>% filter(chip == "SS200001075BR", str_detect(subclass,"IT")&str_detect(subclass,"CAR3",negate = T)),aes(x,y,color = subclass)) + geom_point() + facet_grid(~subclass)} %>% ggsave(filename = "fig3c.pdf")

