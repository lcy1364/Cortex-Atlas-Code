library(tidyverse)
library(qs)
library(parallel)
library(magrittr)
library(RANN)
library(ggridges)
library(dendextend)
library(ggpubr)
devtools::load_all("~/spacexr-master/")
setwd("~/cortex/fig1/")

spatialCellMeta <- read_csv("../STEREO/spatialCellMeta.csv")

subclass_color <-c(AST = "#665C47", ENDO = "#604B47", ET = "#CEC823", CHANDELIER = "#E25691", 
                   `L2-L3 IT LINC00507` = "#07D8D8", `L3-L4 IT RORB` = "#09B2B2", 
                   `L4-L5 IT RORB` = "#69B199", `L6 CAR3` = "#898325", `L6 CT` = "#2D8CB8", 
                   `L6 IT` = "#19E3BE", L6b = "#7944AA", LAMP5 = "#D96C07", MICRO = "#513577", 
                   NDNF = "#f2798d", NP = "#93C43B", OLIGO = "#5E8A79", OPC = "#2E3E39", 
                   PAX6 = "#E96DF2", PVALB = "#FF2D4E", SST = "#F2B003", VIP = "#9C28CC", 
                   VLMC = "#697255")

region_color <- c(FPPFC = "#3F4587", DLPFC = "#8562AA", VLPFC = "#EC8561", M1 = "#B97CB5", 
                  S1 = "#D43046", S1E = "#F0592B", PoCG = "#ED4A96", SPL = "#593C97", 
                  SMG = "#A54486", AG = "#FBDE13", V1 = "#299FAE", ITG = "#75CCE3", 
                  STG = "#0C6939", ACC = "#0D9547")

#fig1c
spatialDist <- spatialCellMeta %>% group_by(chip) %>% group_split()
foo = spatialDist[[1]]
mclapply(spatialDist,mc.cores = length(spatialDist),FUN = function(foo){
  ggsave(str_c(foo$region[[1]],"_",foo$chip[[1]],".png"),height = 500,width = 500,dpi = 0.5,limitsize = F,
         ggplot(foo,aes(x,y,color = subclass)) + geom_point(size = 0.01) + cowplot::theme_nothing() + scale_color_manual(values = subclass_color) + coord_fixed())
})
