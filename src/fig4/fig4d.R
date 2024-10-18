library(qs)
library(tidyverse)
library(Seurat)
library(magrittr)
library(ggwordcloud)
setwd("~/cortex/fig4/")

spatialMeta <- read_csv("spatialMeta.csv")
spatialMeta %<>% filter(str_detect(subclass,"SST"))
spatialMeta$deep_shallow <- ifelse(sst$cluster %in% c("SST_16", "SST_36", "SST_4", "SST_23", "SST_3", "SST_13", "SST_28", 
                                                      "SST_37", "SST_11", "SST_6", "SST_8", "SST_20", "SST_29", "SST_18", 
                                                      "SST_22", "SST_27", "SST_15", "SST_9", "SST_1"),yes = "deep",no = "shallow")

ggplot(spatialMeta,mapping = aes(x,y,color = deep_shallow)) + geom_point() + facet_wrap(~chipID,deep_shallow)



