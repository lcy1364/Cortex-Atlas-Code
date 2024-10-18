library(tidyverse)
library(qs)
library(parallel)
library(sf)
library(patchwork)
library(magrittr)
setwd("~/data/STEREO/AnalysisPlot/")
layerCut <- read.csv("./range230826.csv")


chipIDS <- read.delim("./cortex_selected") %$% chip
chipRegion <- read.delim("./cortex_selected") %>% {setNames(.$region,.$chip)}

chipID <- chipIDs[[1]]

mainLayerFiles <-
  list.files("./3_sf//", "qs", full.names = T)

subclasses <-
  c(
    "AST",
    "ENDO",
    "ET",
    "L2-L3 IT LINC00507",
    "L3-L4 IT RORB",
    "L4-L5 IT RORB",
    "L6 CAR3",
    "L6 CT",
    "L6 IT",
    "L6b",
    "LAMP5",
    "MICRO",
    "NDNF",
    "NP",
    "OLIGO",
    "OPC",
    "PAX6",
    "PVALB",
    "SST",
    "VIP",
    "VLMC"
  )

subclass_color <-
  c(
    AST = "#665C47",
    ENDO = "#604B47",
    ET = "#CEC823",
    `L2-L3 IT LINC00507` = "#07D8D8",
    `L3-L4 IT RORB` = "#09B2B2",
    `L4-L5 IT RORB` = "#69B199",
    `L6 CAR3` = "#898325",
    `L6 CT` = "#2D8CB8",
    `L6 IT` = "#19E3BE",
    L6b = "#7944AA",
    LAMP5 = "#D96C07",
    MICRO = "#513577",
    NDNF = "#f2798d",
    NP = "#93C43B",
    OLIGO = "#5E8A79",
    OPC = "#2E3E39",
    PAX6 = "#E96DF2",
    PVALB = "#FF2D4E",
    SST = "#F2B003",
    VIP = "#9C28CC",
    VLMC = "#697255"
  )
chipID = "B02221E6"

cellFiles <- list.files("./17_density/","qs", full.names = T)

cutData <-
  mclapply(chipIDS, mc.cores = length(chipIDs), function(chipID) {
    tmpRCTD <- cellFiles %>% stringr::str_subset(chipID) %>% qread()
    tmpRCTD %<>% filter(!out)
    mainLayer <- str_subset(mainLayerFiles, chipID) %>% qread
    mainLayer$geometry <- mainLayer$geometry / 2400 * 26460
    
    polygonMatrix <-
      layerCut %>% filter(sample == chipID) %>% select(x, y) %>% mutate(x = 1000 *
                                                                          x, y = 1000 * y)
    polygonMatrix[nrow(polygonMatrix) + 1, ] <- polygonMatrix[1, ]
    polygons <-
      polygonMatrix %>% as.matrix() %>% st_linestring() %>% {
        st_polygon(list(.))
      } %>% st_sfc()
    
    # annoFiles
    
    
    layerLineWithin <-
      st_intersection(mainLayer, polygons)  #%>% st_difference(polygons)
    
    
    glbCut <-
      ggplot() + geom_sf(
        data = tmpRCTD,
        mapping = aes(geometry = contour,fill = first_type, color = first_type) 
      ) + scale_color_manual(values = subclass_color) + scale_fill_manual(values = subclass_color)  + guides(color = "none", fill = "none") +  ggnewscale::new_scale_color() + 
      geom_sf(
        data = mainLayer,
        linewidth = 0.5,
        mapping = aes(geometry = geometry, color = color),
        fill = NA
      )  +  guides(color = "none") +
      geom_sf(
        data = polygons,
        fill = NA,
        color = "red",
        linewidth = 0.4
      ) + cowplot::theme_map() + guides(color = "none", fill = "none") 
    
    subCut <-  ggplot() + geom_sf(
      data = tmpRCTD[tmpRCTD$inCut,],
      mapping = aes(geometry = contour,fill = first_type, color = first_type)
    ) + scale_color_manual(values = subclass_color) + scale_fill_manual(values = subclass_color)  +  guides(color = guide_legend(override.aes = list(size = 1))) + ggnewscale::new_scale_color() + 
      geom_sf(
        data = layerLineWithin,
        linewidth = 0.5,
        mapping = aes(geometry = geometry, color = color),
        fill = NA
      )  +  guides(color = "none") +
      geom_sf(
        data = polygons,
        fill = NA,
        color = "red",
        linewidth = 0.4
      ) + cowplot::theme_map()+ Seurat::SpatialTheme()
    
    subclass <-  ggplot() + geom_point(
      data = tmpRCTD %>% filter(inCut,first_type %in% c("L6 IT","L6 CAR3","L6b","L6 CT","OLIGO")) %>% filter(first_type == "OLIGO" | mainLayer %in% c("L6","L5")),
      mapping = aes(x,y,fill = first_type, color = first_type), size = 0.5
    ) + scale_color_manual(values = subclass_color) + scale_fill_manual(values = subclass_color)  +  guides(color = guide_legend(override.aes = list(size = 5),nrow = 1)) + ggnewscale::new_scale_color() + 
      geom_sf(
        data = layerLineWithin,
        linewidth = 0.5,
        mapping = aes(geometry = geometry, color = color),
        fill = NA
      )  +  guides(color = "none") +
      geom_sf(
        data = polygons,
        fill = NA,
        color = "red",
        linewidth = 0.4
      ) + cowplot::theme_map()+ Seurat::SpatialTheme() + facet_wrap( ~ first_type, ncol = 5)
    
    ggsave(
      str_c(
        "./17_density//","fig8_",
        names(chipIDs)[str_detect(chipIDs, chipID)],
        "_",
        chipID,
        "_together.pdf"
      ),
      patchwork::wrap_plots(list(glbCut, subCut), nrow = 1) / subclass  + patchwork::plot_layout(heights = c(3,1.5)) +  plot_annotation(str_c(names(chipIDs)[str_detect(chipIDs, chipID)],  " ", chipID))
      ,
      width = 10,
      height = 7,
      dpi = 150
    )
    ggsave(
      str_c(
        "./17_density//","fig8_",
        names(chipIDs)[str_detect(chipIDs, chipID)],
        "_",
        chipID,
        "_together.png"
      ),
      patchwork::wrap_plots(list(glbCut, subCut), nrow = 1) / subclass  + patchwork::plot_layout(heights = c(3,1.5)) +  plot_annotation(str_c(names(chipIDs)[str_detect(chipIDs, chipID)],  " ", chipID))
      ,
      width = 10,
      height = 7,
      dpi = 150
    )
  })
