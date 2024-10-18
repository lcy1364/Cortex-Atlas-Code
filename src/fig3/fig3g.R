library(tidyverse)
library(EnhancedVolcano)
library(magrittr)
library(ggprism)
setwd("~/cortex/fig3/")

resFiles <- list.files("./fig3g/","neuronsubclass_results.csv",full.names = T,recursive = T)

region_color <- c("#2F4587", "#8562AA", "#EC8561", "#B97CB5", "#D43046", "#F0592B", 
                  "#ED4A96", "#593C97", "#A54486", "#FBDE13", "#299FAE", "#75CCE3", 
                  "#0C6939", "#0D9547") %>% setNames(c("FPPFC", "DLPFC", "VLPFC", "M1", "S1", "S1E", "PoCG", "SPL",
                                                       "SMG", "AG", "V1", "ITG", "STG", "ACC"))

foo <- resFiles[[1]]

tmp <- read_csv(foo)

res <- lapply(resFiles,FUN = function(foo){
  res <- read_csv(foo) %>% filter(str_detect(Covariate,"Treatment") ) %>% select(`Final Parameter`,`log2-fold change`,`Cell Type`) %>% group_by(`Cell Type`) %>% summarise(FinalParamer = -mean(`Final Parameter`),Log2FC = -mean(`log2-fold change`))
  res$region <- foo  %>% dirname %>% basename()
  res
}) %>% purrr::reduce(rbind)

res$label = str_c(res$`Cell Type`, res$region)

res %<>% filter(`Cell Type` %in% c( "ET", "L2-L3 IT LINC00507", "L3-L4 IT RORB","L4-L5 IT RORB", "L6 CAR3", "L6 CT", "L6 IT", "L6B", "NP"))

pdf("fig3g.pdf",height = 6, width = 7)
ggplot(res,aes(FinalParamer,abs(Log2FC),color = region,label = label)) + geom_text(data = res %>% filter(region == "ACC")) + geom_point() + theme_prism() + scale_color_manual(values = region_color) + coord_fixed()
dev.off()
