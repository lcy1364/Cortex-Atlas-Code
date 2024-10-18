library(tidyverse)
library(qs)
library(parallel)
library(magrittr)
library(RANN)
library(ggridges)
library(dendextend)
library(ggpubr)

# fig i
ast_distribution <- spatialCellMeta %>% filter(subclass == "AST")

countByRegion <- ast_distribution %>% filter(str_detect(cluster,"5|4")  )%>% group_by(chip,region,cluster) %>% summarise(count = n()) 
pdf("ast_count_boxplot.pdf",width = 6,height = 3)
ggplot(countByRegion ,aes( x = fct_reorder(region,count,.desc = T), y = count, color = region)) + geom_boxplot() + facet_wrap(~cluster) + 
  cowplot::theme_cowplot() + theme(legend.position = "none") +  rotate_x_text(angle = 45) + ggpubr::labs_pubr() + ggpubr::theme_cleveland() +
  geom_hline(yintercept = mean(countByRegion$count), linetype = 2) + xlab("count") + 
  stat_compare_means(label = "p.signif", method = "wilcox.test",
                     ref.group = ".all.", hide.ns = TRUE) + scale_color_manual(values = region_color)   
dev.off()

# fig j
densityByRegion <- ast_distribution %>% filter(str_detect(cluster,"5|4")) %>% left_join(annoArea_total,by = c("chip" = "chip")) %>% group_by(chip,region,cluster,area) %>% summarise(count = n())  %>% mutate(density = count / area)
# densitystatByRegion <- densityByRegion %>% group_by(region,cluster) %>% summarise(sd = sd(density),mean = mean(density)) 
densityByRegion$region %<>% fct_relevel("ACC",after = 0)

pdf("ast_density_boxplot.pdf",width = 6,height = 3)

ggplot(densityByRegion ,aes( x = fct_reorder(region,density,.desc = T), y = density, color = region)) + geom_boxplot() + facet_wrap(~cluster) + 
  cowplot::theme_cowplot() + theme(legend.position = "none") +  rotate_x_text(angle = 45)+  ggpubr::labs_pubr() + ggpubr::theme_cleveland()+
  geom_hline(yintercept = mean(densityByRegion$density), linetype = 2)+  xlab("density") + 
  stat_compare_means(label = "p.signif", method = "wilcox.test",
                     ref.group = ".all.", hide.ns = TRUE) + scale_color_manual(values = region_color)     # Pairwise comparison against all
dev.off()


# fig l
drawChips <- spatialCellMeta %>% select(region,chip) %>% distinct() %>% filter(region %in%  c("ACC","S1")) %$% chip #,collapse  = "|")
sf <- list.files("~/DATA/data/STEREO/AnalysisPlot/1_mainLayer/","qs",full.names = T) %>% keep(~str_detect(., str_c(drawChips,collapse ="|")))

# fig m
foo = drawChips[[1]]
mclapply(drawChips,mc.cores = length(drawChips),FUN = function(foo){
  shape <- sf %>% keep(~str_detect(.,foo)) %>% qread()
  shape$geometry = shape$geometry/2400*26460/2  
  
  ggsave(str_c("ast_",chipList[foo],"_",foo,"_astro_4",".pdf"),height = 5,width = 5,
         ggplot() + geom_point(data = spatialCellMeta %>% dplyr::filter(chip == foo,subclass == "AST",str_detect(cluster,"4") ) %>% mutate(x = x/2, y = y/2),aes(x,y),size = 0.5, color = "grey20") +
           xlab("micrometer") + ylab("micrometer") + labs(title = chipList[foo]) + ggpubr::labs_pubr() + coord_sf(xlim = c(0,12000),ylim = c(0,12000)) + scale_x_continuous(breaks = seq(0,12000,3000)) + scale_y_continuous(breaks = seq(0,12000,3000)) +
           geom_sf(data = shape %>% filter(layer %in% c("WM", "L3", "L6", "L5", "ARACHNOID", "L2", "L1", "L4")),mapping = aes(geometry = geometry,color = layer),fill = NA)  + scale_color_manual(values = setNames(shape$color,shape$layer))  + cowplot::theme_minimal_grid() + theme(legend.position = "none")
  )
  
  ggsave(str_c("ast_",chipList[foo],"_",foo,"_astro_5",".pdf"),height = 5,width = 5,
         ggplot() + geom_point(data = spatialCellMeta %>% dplyr::filter(chip == foo,subclass == "AST",str_detect(cluster,"5") ) %>% mutate(x = x/2, y = y/2),aes(x,y),size = 0.5, color = "grey20") +
           xlab("micrometer") + ylab("micrometer") + labs(title = chipList[foo]) + ggpubr::labs_pubr() + coord_sf(xlim = c(0,12000),ylim = c(0,12000)) + scale_x_continuous(breaks = seq(0,12000,3000)) + scale_y_continuous(breaks = seq(0,12000,3000)) +
           geom_sf(data = shape %>% filter(layer %in% c("WM", "L3", "L6", "L5", "ARACHNOID", "L2", "L1", "L4")),mapping = aes(geometry = geometry,color = layer),fill = NA)  + scale_color_manual(values = setNames(shape$color,shape$layer))  + cowplot::theme_minimal_grid() + theme(legend.position = "none")
  )
  
  
})