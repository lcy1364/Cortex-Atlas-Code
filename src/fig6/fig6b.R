library(qs)
library(magrittr)
library(MASS)
library(tidyverse)
library(sf)
library(Seurat)
library(patchwork)
setwd("/home/luomeng/data/STEREO/AnalysisPlot//")
fileFiles <- list.files("~/data/STEREO/cellbin","qs")

layerColors <-
  ggsci::pal_aaas()(10) %>% sort %>% {
    c("#CCCCCC", .)
  } %>% alpha(alpha = 1) %>% setNames(c("ARACHNOID", str_c("L", 1:6), "WM"))
# layerColors[[4]] <- "#FF2D4EFF"

mainLayerFiles <-
  list.files("/home/luomeng/data/STEREO/AnalysisPlot/1_mainLayer/",
             "_shape.qs",
             full.names = T)

{subclasses <-
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

layerColor <-
  c(
    "#e48171",
    "#bf3729",
    "#e69b00",
    "#f5bb50",
    "#ada43b",
    "#355828",
    "#17154f",
    "#2f357c",
    "#6c5d9e",
    "#9d9cd5",
    "#b0799a",
    "#f6b3b0"
  ) %>% setNames(c("L1", "L4", "L3", "L2", "L6", "L5", "WM"))

region <-
  c(
    'FPPFC',
    'DLPFC',
    'VLPFC',
    'M1',
    'S1',
    'S1E',
    'PoCG',
    'SPL',
    'SMG',
    'AG',
    'V1',
    'ITG',
    'STG',
    'ACC'
  )

region_cata <-
  list(
    c("FPPFC", "DLPFC", "VLPFC"),
    c("ACC"),
    c("S1E"),
    c("STG", "ITG"),
    c("V1"),
    c("SMG", "AG", "SPL", "PoCG"),
    c("S1", "M1")
  )

region_color <-
  c(
    '#3F4587',
    '#8562AA',
    '#EC8561',
    '#B97CB5',
    '#D43046',
    '#F0592B',
    '#ED4A96',
    '#593C97',
    '#A54486',
    '#FBDE13',
    '#299FAE',
    '#75CCE3',
    '#0C6939',
    '#0D9547'
  ) %>% setNames(region)

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
  )}

layerCut <- read.csv("../9_Layers_analysis/range230826.csv")

file = fileFiles[[1]]

countList <- parallel::mclapply(fileFiles, mc.cores = length(fileFiles), function(file){ 

  chipID = file %>% basename() %>% str_extract("[A-Z0-9]*")
  print(chipID)
  filterData <- qread(file)
  filterData %<>% as.data.frame()
  filterData %<>% filter(out == "FALSE")
  
  # update layer
  mainLayer <- qread(mainLayerFiles %>% str_subset(chipID))
  mainLayer$geometry <- mainLayer$geometry / 2400 * 26460
  
  polygonMatrix <-
    layerCut %>% filter(sample == chipID) %>% dplyr::select(x, y) %>% mutate(x = 1000 * x, y = 1000 * y) 
  
  polygonMatrix[nrow(polygonMatrix) + 1, ] <- polygonMatrix[1, ]
  polygons <-
    polygonMatrix %>% as.matrix() %>% st_linestring() %>% {
      st_polygon(list(.))
    } %>% st_sfc()
  

  # which cells inside
  filterData$inCut <-  filterData$contour %>%  st_within(polygons, sparse = F)
  cutRegionShape <- st_intersection(mainLayer, polygons)  
  
  ggplot(filterData %>% filter(inCut),aes(x,y)) + geom_point(size = 0.0) + facet_wrap(~first_type) + geom_sf(data = cutRegionShape,inherit.aes = F,fill = NA)
  
  cutData <- filterData %>% filter(inCut)
  region_count <- cutData %>% group_by(first_type) %>% filter(str_detect(Layer,"L[1-6]"))%>% summarise(count = n()) %>% mutate(per = count/sum(count))
  region_count$region <- cortexMeta[chipID]
  region_count$chip <- chipID  
  region_count
  
})
countList %<>% purrr::reduce(rbind)
count_region <- countList %>% purrr::reduce(rbind) %>% group_by(region,first_type) %>% summarise(per = mean(per))
count_region %<>% pivot_wider(values_from = "per",names_from = "first_type",values_fill = 0)
count_region %<>% column_to_rownames("region")

pheatmap::pheatmap(count_region,scale = "column",color = MetBrewer::MetPalettes$Benedictus[[1]] %>% rev,filename = "./fig7B.pdf",width = 15,height = 5)

count_region %>% write.csv("count_region.csv")

count_region <- read.csv("count_region.csv",row.names = 1)
