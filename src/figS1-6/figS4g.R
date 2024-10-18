library(qs)
library(sf)
setwd("/media/desk16/luomeng/data/STEREO/AnalysisPlot/")
qsFiles <- list.files("~/data/STEREO/AnalysisPlot/","qs$",full.names = T) %>% str_subset("cut")


subclass_color <- c(AST = "#665C47", ENDO = "#604B47", ET = "#CEC823", `L2-L3 IT LINC00507` = "#07D8D8", 
                    `L3-L4 IT RORB` = "#09B2B2", `L4-L5 IT RORB` = "#69B199", `L6 CAR3` = "#898325", 
                    `L6 CT` = "#2D8CB8", `L6 IT` = "#19E3BE", L6b = "#7944AA", LAMP5 = "#D96C07", 
                    MICRO = "#513577", NDNF = "#f2798d", NP = "#93C43B", OLIGO = "#5E8A79", 
                    OPC = "#2E3E39", PAX6 = "#E96DF2", PVALB = "#FF2D4E", SST = "#F2B003", 
                    VIP = "#9C28CC", VLMC = "#697255")

file = qsFiles[[1]]
for (file in qsFiles) {
  cutData = qs::qread(file)

  chipID = file %>% basename() %>% str_remove("_RCTD_seurat.cut.qs*")
  cutData$first_type_prob <- purrr::map_dbl(1:nrow(cutData@meta.data), function(i){cutData@meta.data[i,cutData@meta.data$first_type[[i]] %>% as.character()]})
  cutData$max_type_prob <- purrr::map_dbl(1:nrow(cutData@meta.data), function(i){cutData@meta.data[i,cutData@meta.data$max_type[[i]] %>% as.character()]})
  
  cutData$first_major[cutData$first_type %in% non_subclass] = "non"
  
  cutData <- cutData[,cutData$spot_class != "reject"]
  cutData <- cutData[,cutData$nFeature_Spatial >= 100]
  
  nonData <- cutData[,cutData$first_major == "non"& cutData$first_type_prob > 0.01 & cutData$spot_class !="reject"]
  
  tmp = ggplot(nonData@meta.data %>% st_as_sf(wkt = "contour") , aes(x,y,fill = first_type, color = first_type) ) + geom_sf( linewidth = 0.01) + geom_sf(data = nonData@misc$cutRegionShape,fill = NA,inherit.aes = F)  + 
    facet_wrap(~first_type,ncol = 6) + scale_fill_manual(values = subclass_color) + scale_color_manual(values = subclass_color) +   
    cowplot::theme_map()  + 
    theme(panel.background = element_rect(fill= "black", color = "black"),plot.background = element_rect(fill= "black", color = "black"),legend.position = "none")
  ggsave(str_c(chipID, ".pdf"),tmp,width = 10,height = 2)
}
