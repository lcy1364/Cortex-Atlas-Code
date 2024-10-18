library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggsignif)
setwd("~/cortex/figS1-6/")

spatialCellMeta <- read_csv("../STEREO/spatialCellMeta.csv")

ast_distribution <- spatialCellMeta %>% filter(subclass == "AST")

countByRegion <- ast_distribution %>% filter(str_detect(cluster,"5|4")  )%>% group_by(chip,region,cluster) %>% summarise(count = n()) 

pdf("figS4i-j.pdf",width = 6,height = 3)
ggplot(countByRegion ,aes( x = fct_reorder(region,count,.desc = T), y = count, color = region)) + geom_boxplot() + facet_wrap(~cluster) + 
  cowplot::theme_cowplot() + theme(legend.position = "none") +  ggpubr::rotate_x_text(angle = 45) + ggpubr::labs_pubr() + ggpubr::theme_cleveland() +
  geom_hline(yintercept = mean(countByRegion$count), linetype = 2) + xlab("count") + 
  stat_compare_means(label = "p.signif", method = "wilcox.test",
                     ref.group = ".all.", hide.ns = TRUE) + scale_color_manual(values = region_color)   
dev.off()
# draw acc ast

drawChips <- spatialCellMeta %>% select(region,chip) %>% distinct() %>% filter(region %in%  c("ACC","S1")) %$% chip #,collapse  = "|")

foo = drawChips[[1]]
mclapply(drawChips,mc.cores = length(drawChips),FUN = function(foo){

  ggsave(str_c("ast_",chipList[foo],"_",foo,"_astro_4",".pdf"),height = 5,width = 5,
         ggplot() + geom_point(data = spatialCellMeta %>% dplyr::filter(chip == foo,subclass == "AST",str_detect(cluster,"4") ) %>% mutate(x = x/2, y = y/2),aes(x,y),size = 0.5, color = "grey20") +
           xlab("micrometer") + ylab("micrometer") + labs(title = chipList[foo]) + ggpubr::labs_pubr() + scale_x_continuous(breaks = seq(0,12000,3000)) + scale_y_continuous(breaks = seq(0,12000,3000))  + scale_color_manual(values = setNames(shape$color,shape$layer))  + cowplot::theme_minimal_grid() + theme(legend.position = "none")
  )
  
  ggsave(str_c("ast_",chipList[foo],"_",foo,"_astro_5",".pdf"),height = 5,width = 5,
         ggplot() + geom_point(data = spatialCellMeta %>% dplyr::filter(chip == foo,subclass == "AST",str_detect(cluster,"5") ) %>% mutate(x = x/2, y = y/2),aes(x,y),size = 0.5, color = "grey20") +
           xlab("micrometer") + ylab("micrometer") + labs(title = chipList[foo]) + ggpubr::labs_pubr() + coord_sf(xlim = c(0,12000),ylim = c(0,12000)) + scale_x_continuous(breaks = seq(0,12000,3000)) + scale_y_continuous(breaks = seq(0,12000,3000)) +
           scale_color_manual(values = setNames(shape$color,shape$layer))  + cowplot::theme_minimal_grid() + theme(legend.position = "none")
  )
})
