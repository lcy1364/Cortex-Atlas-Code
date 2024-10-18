library(qs)
library(magrittr)
library(MASS)
library(tidyverse)
library(sf)
library(patchwork)
setwd("/home/luomeng/data/STEREO/AnalysisPlot/16_Plot/")

spatialCellMeta <- read_csv("spatialCellMeta.csv")

subclass_color <- source("./glossary.R")

pdf("fig4b.pdf",height = 5,width = 20)
ggplot(spatialCellMeta %>% filter(subclass %in%  c("NDNF","LAMP5","PAX6","VIP")) , aes( x,y, color = subclass )) + facet_wrap(~subclass,nrow = 1) 
dev.off()
