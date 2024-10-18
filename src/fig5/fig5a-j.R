library(igraph)
library(ggraph)
library(tidyverse)
library(magrittr)
library(tidygraph)
library(sf)
library(ggforce)
library(patchwork)
library(qs)
setwd("~/DATA/BRAIN/STEREO/frequentGraph/plot/")

# 15um ------------
project = "../soma15nn15.network"
wsize = 100
chipID_region <- read.delim("../../cortex") %>% {setNames(.$region,nm = .$chip)}

subclass_color <-c(AST = "#665C47", ENDO = "#604B47", ET = "#CEC823", `L2-L3 IT LINC00507` = "#07D8D8",`L3-L4 IT RORB` = "#09B2B2", `L4-L5 IT RORB` = "#69B199", `L6 CAR3` = "#898325",`L6 CT` = "#2D8CB8", `L6 IT` = "#19E3BE", L6b = "#7944AA", LAMP5 = "#D96C07",MICRO = "#513577", NDNF = "#f2798d", NP = "#93C43B", OLIGO = "#5E8A79",OPC = "#2E3E39", PAX6 = "#E96DF2", PVALB = "#FF2D4E", SST = "#F2B003",VIP = "#9C28CC", VLMC = "#697255")

subclass_class <- c(AST = "NON", ENDO = "NON", ET = "EXC", `L2-L3 IT LINC00507` = "EXC",
`L3-L4 IT RORB` = "EXC", `L4-L5 IT RORB` = "EXC", `L6 CAR3` = "EXC",
                    `L6 CT` = "EXC", `L6 IT` = "EXC", NP = "EXC", L6b = "EXC", LAMP5 = "INH",
                    MICRO = "NON", NDNF = "INH", OLIGO = "NON", OPC = "NON", PAX6 = "INH",
                    PVALB = "INH", SST = "INH", VIP = "INH", VLMC = "NON")

subclass_class.fct <- subclass_class %>% factor(levels = c("NON","INH","EXC"))


cellOrder <- c(  "ET",  "L6 CT","L6 CAR3","L6b", "L6 IT", "NP","L4-L5 IT RORB", "L3-L4 IT RORB", "L2-L3 IT LINC00507",
               "VLMC", "ENDO",  "OLIGO", "OPC","AST", "MICRO",
               "LAMP5", "NDNF","PAX6", "VIP","PVALB", "SST") # 16-21
cellOrder.num <-setNames(1:length(cellOrder),nm = cellOrder)

## some plot function------------

expandLmits <- function(x,z =1.2){
  r = (max(x)-min(x))/2
  m = (max(x)+min(x))/2
  c(m - z*r,m + z*r)
}


find_leaf_vertices <- function(g) {
  vertex_degrees <- degree(g)
  leaf_vertices <- which(vertex_degrees == min(vertex_degrees))
  leaf_vertex_names <- V(g)[leaf_vertices]
  sorted_leaf_vertex_names <- leaf_vertex_names[order(cellOrder.num[leaf_vertex_names %>% names])]
  return(sorted_leaf_vertex_names)
}


generate_graph_name <- function(g) {
  sorted_leaf_vertices <- find_leaf_vertices(g)
  start_vertex <- sorted_leaf_vertices[1]
  current_vertex <- start_vertex
  visited <- c()
  vertex_order <- c(current_vertex)

  while(TRUE) {
    visited <- c(visited, current_vertex)
    neighbors <- neighbors(g, current_vertex)
    next_vertex <- neighbors[!neighbors %in% visited]

    if (length(next_vertex) == 0) break
    next_vertex <- next_vertex[1]

    vertex_order <- c(vertex_order, next_vertex)
    current_vertex <- next_vertex
  }

  paste(vertex_order %>% names, collapse = "--")
}

find_max_distance_combination <- function(df, iterations = 100) {
  if (!all(c("x", "y", "category") %in% names(df))) {
    stop("DataFrame must contain x, y, and category columns")
  }
  euclidean_distance <- function(point1, point2) {
    sqrt(sum((point1 - point2)^2))
  }

  max_distance <- -Inf
  best_combination <- NULL

  for (i in 1:iterations) {
    selected_points <- df %>%
      group_by(category) %>%
      sample_n(1) %>%
      ungroup()

    total_distance <- sum(unlist(sapply(1:(nrow(selected_points) - 1), function(i) {
      sapply((i + 1):nrow(selected_points), function(j) {
        euclidean_distance(selected_points[i, c("x", "y")], selected_points[j, c("x", "y")])
      })
    })))

    if (total_distance > max_distance) {
      max_distance <- total_distance
      best_combination <- selected_points
    }
  }
  return(best_combination)
}
# end ---------

graphList <- qs::qread("graphList.qs")

#fig5a ----------
chipID_region <- read.delim("../../cortex") %>% {setNames(.$region,nm = .$chip)}

selected_network <- c("345","156","171","596")

tmpGraphList <- graphList
names(tmpGraphList) <- map(tmpGraphList,.f = ~.$org_id) %>% unlist
subnet = selected_network[[1]]

targeted_Units <- parallel::mclapply(selected_network,mc.cores = length(selected_network),function(subnet) {
  
  tmpGraph <- tmpGraphList[[subnet]]
  tmpName <- generate_graph_name(tmpGraph)
  majorType <- subclass_class[tmpGraph %>% V %>% names] %>% unique %>% sort %>% str_c(collapse = "-")
  belongings <- all.net[ as.numeric(tmpGraph$belongings) +1] 
  
  belongings <- purrr::map(1:length(belongings),function(x) {
    tmp = belongings[[x]]
    tmp$org_id = (as.numeric(tmpGraph$belongings) + 1)[[x]]
    tmp$patternID = subnet
    tmp$chip <- names(belongings)[[x]]
    tmp$region <- chipID_region[tmp$chip]
    tmp$VertexNum <- length(V(tmp))
    tmp$majorType <- majorType
    tmp$networkType <- tmpName
    layerDis <- all.chipData[[tmp$chip]][["mainLayer"]][V(tmp) %>% names %>% as.numeric()]  %>% table
    tmp$Layer <-names(layerDis)[which.max(layerDis)][[1]] # 增加一个layer
    tmp
  })
  
}) %>% purrr::reduce(c)

colorSet <- MetBrewer::met.brewer(name="OKeeffe1", n=length(selected_network), type="discrete")  %>% as.character()%>% setNames(targeted_Units %>% map(~.$networkType) %>% unlist %>% unique)

selectedChips <- c("B02222E1","B02204C5","B01012A3","B02223C3")
chip = selectedChips[[1]]

parallel::mclapply(selectedChips,mc.cores = length(selectedChips),function(chip) {
  
  cellData <- all.chipData[[chip]]
  
  
  locations <- targeted_Units[targeted_Units %>% purrr::map(~.$chip == chip) %>% unlist]
  
  cellData$shape_zoom <- (cellData$contour - st_centroid(cellData$contour))*1.5 + st_centroid(cellData$contour)
  
  selected_units <- 1:length(locations) %>% map(function(x){
    subLoc = locations[[x]]
    subLoc$targetedCell <- cellData[V(subLoc) %>% names %>% as.numeric(),]
    merged_object <- st_union(subLoc$targetedCell$contour)
    subLoc$unit_centroid <- st_centroid(merged_object)[[1]]
    subLoc$convex_hull <- st_convex_hull(merged_object)
    subLoc$id <- x
    subLoc
  })
  
  selected_units <- selected_units[selected_units %>% map_dfr(~data.frame(category = .$networkType, vertexNum = .$VertexNum, id = .$id, x = .$unit_centroid[[1]], y = .$unit_centroid[[2]]))  %>% group_by(category) %>%  slice_min(n = 1,order_by = vertexNum) %>%  find_max_distance_combination() %$% id]
  
  # library(arrow)
  p_all <- ggplot( ) +
    geom_sf(data = cellData, mapping = aes(geometry =  shape_zoom,fill = first_type),color = NA,inherit.aes = F) +
    scale_fill_manual(values = alpha(subclass_color)) + cowplot::theme_nothing() + theme(plot.background = element_rect(fill = "black"))
  p_all2 <- p_all
  
  p_unit <- list()
  for (subLoc in selected_units) {
    xmin <- subLoc$unit_centroid[[1]]  - wsize*2.5
    xmax <- subLoc$unit_centroid[[1]]  + wsize*2.5
    ymin <- subLoc$unit_centroid[[2]]  - wsize*2.5
    ymax <- subLoc$unit_centroid[[2]]  + wsize*2.5
    print(c(xmin,xmax,ymin,ymax))
    
    p_all2 <- p_all2 +
      annotate("rect",xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = NA, color = colorSet[subLoc$networkType], linewidth = 1) +
      annotate("text",
               x =  subLoc$unit_centroid[[1]],
               y = subLoc$unit_centroid[[2]],
               label = subLoc$networkType,size = 10, color = "white",hjust = 0.75, vjust = 2)
    
    p_all <- p_all +
      annotate("rect",xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = NA, color = colorSet[subLoc$networkType], linewidth = 1) #+
    
    p_unit[[subLoc$networkType]] <- ggplot() +
      geom_sf(data = cellData %>% filter(x < xmax, x > xmin, y < ymax, y > ymin), mapping = aes(geometry =  contour,fill = first_type),inherit.aes = F) +
      geom_sf(data = (subLoc$convex_hull[[1]]-st_centroid(subLoc$convex_hull[[1]]))*1.5 +st_centroid(subLoc$convex_hull[[1]]), fill = NA, color = "red",inherit.aes = F,linewidth = 2,linetype = "dashed") +
      geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = NA, color = colorSet[subLoc$networkType], linetype = "dashed",inherit.aes = F,linewidth = 5) +
      cowplot::theme_map() +
      scale_fill_manual(values = alpha(subclass_color,0.5)) + cowplot::theme_nothing()+ theme(plot.background = element_rect(fill = "white"))
    
    
  }
  
  ggsave(str_c(str_c(somaType,"/"),str_c(chipID_region[chip],chip,subLoc$networkType,subLoc$VertexNum,sep = "-"),".pdf",sep = ""),p_unit, height = 10,width = 10,dpi = 100)
  ggsave(str_c(str_c(somaType,"/"),str_c(chipID_region[chip],chip,sep = "-"),".png",sep = ""),p_all, height = 5,width = 5,dpi = 300)
  
  ggsave(str_c(str_c(somaType,"/"),str_c(chipID_region[chip],chip, sep = "-"),".pdf",sep = ""),p_all, height = 5,width = 5,dpi = 300)
  
  ggsave(str_c(str_c(somaType,"/"),str_c(chipID_region[chip],chip,sep = "-"),"_with_label.pdf",sep = ""),p_all2, height = 5,width = 5,dpi = 300)
  
})


# fig5b #==========================================

resTxt_all <- readLines(str_c("../",project,"0.fsm"))
# resTxt <- readLines("./all_15.network.fsm")

head <- ( 1:length(resTxt_all))[str_detect(resTxt_all,"^t")]
end <- ( 1:length(resTxt_all))[str_detect(resTxt_all,"^x")]
ctCodes <- c("AST", "ENDO", "ET", "L2-L3 IT LINC00507", "L3-L4 IT RORB", # 1-5
             "L4-L5 IT RORB", "L6 CAR3", "L6 CT", "L6 IT", "L6b",  # 6-10
             "LAMP5", "MICRO", "NDNF", "NP", "OLIGO", # 11-15
             "OPC", "PAX6", "PVALB", "SST", "VIP", "VLMC") # 16-21

# 1-1 统计
interaction_2 <- parallel::mclapply(1:length(head), mc.cores = 20,function(i){
  tmp <-resTxt_all[head[[i]]:end[[i]]]
  Vertice <- tmp %>% str_subset("^v") %>% map_dfr(.f = function(x){
    str_split_fixed(x," ",n = 3)[2:3] %>% as.numeric() %>% { data.frame(row.names = .[[1]] ,name = .[[2]])}})
  
  
  if (nrow(Vertice) == 2) {
    belongings <- tmp %>% str_subset("^x") %>% str_remove("x ") %>% str_split(" ",simplify = T)
    
    support <- length(belongings)
    
    Edg <- tmp %>% str_subset("^e") %>% map_dfr(.f = function(x){
      str_split_fixed(x," ",n = 4)[2:3] %>% map(~Vertice[.,"name"]) %>% unlist %>% sort %>%
        { data.frame(first = ctCodes[c(.[[1]],.[[2]])],
                     second =ctCodes[c(.[[2]],.[[1]])],
                     netID = i,
                     support = support#,
                     #belongings = tmp %>% str_subset("^x") %>% str_remove("x ")
        )}
    })
  } else { NULL}
  
})%>% purrr::reduce(rbind) %>% arrange(first,second)

interaction_2$first %<>% as.character() %>% factor(levels = cellOrder)
interaction_2$second %<>% as.character() %>% factor(levels = cellOrder)

{ggplot(interaction_2,aes(first,second,fill = support,size = support)) + geom_raster() + ggprism::theme_prism() + scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11,"Spectral") %>% rev,trans = "log") + scale_size_continuous( trans = "log") + theme(axis.text.x = element_text(angle = 90,hjust = 1))} %>% ggsave(filename = str_c("figA.interaction.pdf"),height = 10,width = 11)
# end ---------

# fig5c --------

tmpGraphList <- graphList %>% keep(~length(V(.))==2) %>% keep(~(V(.) %>% names %>% unique %>% length) >1)

names(tmpGraphList) <- map(tmpGraphList,.f = ~.$org_id) %>% unlist
subnet = tmpGraphList[[203]]

targeted_Units.1v1.layer <- parallel::mclapply(tmpGraphList,mc.cores = 20,function(subnet) {
  
  tmpGraph <- subnet
  tmpName <- generate_graph_name(tmpGraph)
  
  if(all(subclass_class[V(subnet) %>% names] == "EXC")) return(NULL)
  if(length(subclass_class[V(subnet) %>% names] %>% unique) < 2) return(NULL)
  
  # non - non
  # non - gaba
  # non - glut
  # gaba - gaba
  # gaba - glut
  
  belongings <- all.net[ as.numeric(tmpGraph$belongings) +1] 
  
  belongings <- purrr::map(1:length(belongings),function(x) {
    tmp = belongings[[x]]
    tmp$org_id = (as.numeric(tmpGraph$belongings) + 1)[[x]]
    tmp$chip <- names(belongings)[[x]]
    tmp$region <- chipID_region[tmp$chip]
    tmp$VertexNum <- length(V(tmp))
    tmp$networkType <- tmpName
    layerDis <- all.chipData[[tmp$chip]][["mainLayer"]][V(tmp) %>% names %>% as.numeric()]  %>% table
    tmp$Layer <-names(layerDis)[which.max(layerDis)][[1]] # 增加一个layer
    tmp
  })
  
}) %>% purrr::reduce(c) %>% keep(~!is.null(.))

# mark region with their four lobes

region_cata <-
  c(FPPFC = "Frontal", DLPFC = "Frontal", VLPFC = "Frontal",
    M1 = "Parietal",  S1 = "Parietal", S1E = "Parietal", PoCG = "Parietal", SPL = "Parietal", SMG = "Parietal", AG = "Parietal",
    V1 = "Occipital",
    ITG = "Temporal",
    STG = "Temporal",
    ACC = "ACC")

brainRegionOrder <- c("Frontal", "Parietal", "Occipital", "Temporal", "ACC")

regionUnitStat.1v1.layer <- targeted_Units.1v1.layer %>% purrr::map_dfr(function(x){
  data.frame(chip = x$chip,network = x$networkType, layer = x$Layer, region = region_cata[chipID_region[x$chip]])}) %>%
  filter(layer != "WM") %>%
  group_by(chip,network,region,layer) %>%
  summarise(freq = n())  %>%
  group_by(chip,network,region) %>%
  mutate(freq = freq/sum(freq))  %>%
  pivot_wider(names_from = "network",values_from = "freq", values_fill = 0) %>%
  mutate(region = factor(region,levels =  c("Frontal", "Parietal", "Occipital", "Temporal", "ACC")))  %>%
  ungroup %>%
  select(-chip) %>%
  group_by(region,layer)  %>%
  summarise_all(mean) %>%
  arrange(layer,region) %>%
  unite(col = "region_layer",layer,region,sep = "--") %>%
  column_to_rownames("region_layer")

colType <- regionUnitStat.1v1.layer %>% colnames %>% str_split_fixed("--",n = 2) %>% data.frame() %>% cbind( regionUnitStat.1v1.layer %>% colnames)
colnames(colType) <- c("first","second","colnames")

colType$first_type <-  factor(subclass_class,levels = c("NON","EXC","INH"))[colType$first]
colType$second_type <-  fct_collapse(factor(subclass_class,levels = c("INH","EXC","NON")) ,OTHER = c("INH", "NON"))[colType$second]
colType$first %<>% factor(levels = cellOrder)
colType$second %<>% factor(levels = cellOrder)

colOrder <- colType %>% arrange(first_type,second_type,first,second) %$% colnames

regionUnitStat.1v1.layer <-regionUnitStat.1v1.layer[regionUnitStat.1v1.layer %>% rownames() %>% str_detect("ARACHNOID|ACC",negate = T),regionUnitStat.1v1.layer %>% colnames() %>% str_detect("L6b|CAR3|ET|ENDO|VLMC",negate = T)]

{ggheatmap::ggheatmap(as.matrix(regionUnitStat.1v1.layer) %>% t,cluster_cols = F,cluster_rows = F,levels_rows = colOrder,show_cluster_rows = F,scale = "none",color =  MetBrewer::met.brewer("Demuth", 20) %>% rev)  +theme(axis.text.x = element_text(angle = 90,face = "bold",size = 10,hjust = 1))} %>%  ggsave(filename = str_c("./","layer_1v1_NN",".pdf",sep = ""), width = dim(regionUnitStat.1v1.layer)[[1]]*0.25, height=  dim(regionUnitStat.1v1.layer)[[2]]*0.15,dpi = 100)


# fig5d ---------
somaType <- "EXC_INH"

graphList.sub <- parallel::mclapply (
  unique(names(graphList)),
  mc.cores = 7,
  FUN = function(type) {
    print(type)
    tmp <- graphList[names(graphList) == type]
    if (length(tmp) > 200)
      tmp <- tmp[1:500]
    tmp <- tmp[(map(tmp,  ~ .$support) %>% unlist) > 30]
    
    org_id <- (1:length(graphList))[names(graphList) == type]  
    
    lapply(1:length(tmp), function(i) {
      layout_matrix <- layout_with_kk(tmp[[i]]) %>% as.data.frame()
      colnames(layout_matrix) <- c("x", "y")
    
    
    tmp2 <- list( lapply(1:length(tmp), function(i) {
      tmp2 <- tmp[[i]]
      tmp2$org_id <- org_id[[i]]
      tmp2
    }))
    names(tmp2) <- type
    tmp2
  }
)  %>% purrr::reduce(c)

selected_network <- c("44", "48",  "52","73",  "90", "130",  "141", "156", "171", "212", "269","294","345")
    
tmpGraphList <- graphList.sub[[somaType]]
names(tmpGraphList) <- map(tmpGraphList,.f = ~.$org_id) %>% unlist
subnet = selected_network[[1]]

dir.create(somaType)

targeted_Units <- parallel::mclapply(selected_network,mc.cores = length(selected_network),function(subnet) {
  tmpGraph <- tmpGraphList[[subnet]]
  
  # 生成图的名字
  tmpName <- generate_graph_name(tmpGraph)
  layout_matrix <- layout_with_kk(tmpGraph) %>% as.data.frame()
  colnames(layout_matrix) <- c("x", "y")
  p_graph <-
    ggraph(tmpGraph, layout_matrix) +
    geom_edge_arc(width = 5,
                  color = "grey80",
                  strength  = 0.05) + geom_node_point(size = 20, aes(color = name)) +
    geom_node_text(
      aes(label  = name),
      color = "darkblue",
      fontface = "bold",
      hjust = "inward",
      vjust = "inward",
      size = 10
    ) + labs(title = str_c(c(subnet, " ", tmpGraph$support), collapse = "")) + scale_color_manual(values = alpha(subclass_color, alpha = 0.9)) + cowplot::theme_map() + theme(legend.position = "none") + coord_cartesian(xlim =   layout_matrix$x %>% expandLmits,
                                                                                                                                                                                                                          ylim = layout_matrix$y %>% expandLmits)
  ggsave(str_c(str_c(somaType,"/"),somaType,subnet,tmpName,tmpGraph$support,".pdf",sep = "_"),p_graph + patchwork::plot_layout(widths = c(1,0.5)), height = 5,width = 5,dpi = 100)
  
  
  belongings <- all.net[ as.numeric(tmpGraph$belongings) +1] # 子网的编号要记得+1，这里面是从0开始计数的，上一级是all.net，all.net的编号芯片+对应meta里的细胞顺序
  
  belongings <- purrr::map(1:length(belongings),function(x) {
    tmp = belongings[[x]]
    tmp$org_id = (as.numeric(tmpGraph$belongings) + 1)[[x]]
    tmp$chip <- names(belongings)[[x]]
    tmp$region <- chipID_region[tmp$chip]
    tmp$VertexNum <- length(V(tmp))
    tmp$networkType <- tmpName
    tmp
  })
  
})%>% purrr::reduce(c)


targeted_Units %>% purrr::map(~.[[1]]$networkType)


qsave(targeted_Units,str_c(somaType,".qs"))

targeted_Units <- qread(str_c(somaType,".qs"))

chipUnitStat <- targeted_Units %>% purrr::map_dfr(function(x){
  data.frame(chip = x$chip,network = x$networkType)
}) %>% group_by(chip,network) %>% summarise(freq = n()) %>% pivot_wider(names_from = "network",values_from = "freq", values_fill = 0) %>% select(-chip)

regionUnitStat <- targeted_Units %>% purrr::map_dfr(function(x){
  data.frame(chip = x$chip,network = x$networkType, region = chipID_region[x$chip])
}) %>% group_by(chip,network,region) %>% summarise(freq = n()) %>% pivot_wider(names_from = "network",values_from = "freq", values_fill = 0) %>% select(-chip) %>% group_by(region) %>% summarise_all(mean) %>%  column_to_rownames("region")

# pheatmap::pheatmap(regionUnitStat)

ggheatmap::ggheatmap(as.matrix(regionUnitStat),cluster_rows = T,cluster_cols = T,scale = "row",color =  MetBrewer::met.brewer("Demuth", 20) %>% rev) %>% ggheatmap::ggheatmap_theme(1,theme = list(theme(axis.text.x = element_text(angle = 90,face = "bold",size = 10,hjust = 1))))  %>%  ggsave(filename = str_c(str_c(somaType,"/"),somaType,".pdf",sep = "_"), height = 8,width = 5,dpi = 100)

# fig5e --------
somaType <- "EXC_INH_NON"
selected_network <- c("408", "429",  "447","452",  "545",  "596", "606")

# selected_network <- c( "429",  "606")

tmpGraphList <- graphList.sub[[somaType]]
names(tmpGraphList) <- map(tmpGraphList,.f = ~.$org_id) %>% unlist
subnet = selected_network[[1]]
dir.create(somaType)

targeted_Units <- parallel::mclapply(selected_network,mc.cores = length(selected_network),function(subnet) {
  
  tmpGraph <- tmpGraphList[[subnet]]
  tmpName <- generate_graph_name(tmpGraph)
  
 belongings <- all.net[ as.numeric(tmpGraph$belongings) +1] 
  
  belongings <- purrr::map(1:length(belongings),function(x) {
    tmp = belongings[[x]]
    tmp$org_id = (as.numeric(tmpGraph$belongings) + 1)[[x]]
    tmp$chip <- names(belongings)[[x]]
    tmp$region <- chipID_region[tmp$chip]
    tmp$VertexNum <- length(V(tmp))
    tmp$networkType <- tmpName
    layerDis <- all.chipData[[tmp$chip]][["mainLayer"]][V(tmp) %>% names %>% as.numeric()]  %>% table
    tmp$Layer <-names(layerDis)[which.max(layerDis)][[1]] # 增加一个layer
    tmp
  })
  
}) %>% purrr::reduce(c)

regionUnitStat <- targeted_Units %>% purrr::map_dfr(function(x){
  data.frame(chip = x$chip,network = x$networkType, region = chipID_region[x$chip])
}) %>% group_by(chip,network,region) %>% summarise(freq = n()) %>% pivot_wider(names_from = "network",values_from = "freq", values_fill = 0) %>% select(-chip) %>% group_by(region) %>% summarise_all(mean) %>%  column_to_rownames("region")

ggheatmap::ggheatmap(as.matrix(regionUnitStat),cluster_rows = T,cluster_cols = T,scale = "row",color =  MetBrewer::met.brewer("Demuth", 20) %>% rev) %>% ggheatmap::ggheatmap_theme(1,theme = list(theme(axis.text.x = element_text(angle = 90,face = "bold",size = 10,hjust = 1))))  %>%  ggsave(filename = str_c(str_c(somaType,"/"),somaType,".pdf",sep = ""), height = 8,width = 5,dpi = 100)



# fig5f --------

graphlistIncludedOLG_NN <- graphList %>% keep(~length(V(.))==2)  %>% keep(~any(str_detect(names(V(.)),"OPC|OLIGO"))) %>%  keep(~any(str_detect(subclass_class[names(V(.))],"INH|EXC")))

RegionGraphlistIncludedOLG_NN <- purrr::map(graphlistIncludedOLG_NN,
                                            function(g){
                                              belongings <- all.net[ as.numeric(g$belongings) +1]  
                                              # tmpName <- generate_graph_name(g)
                                              tmpName <- str_c(ifelse(V(g) %>% names() %>% str_detect("OPC") %>% any,"OPC","OLIGO"),"-",V(g) %>% names() %>% keep(~ . != "OPC" & . != "OLIGO") )
                                              belongings <- purrr::map(1:length(belongings),function(x) {
                                                tmp = belongings[[x]]
                                                tmp$org_id = (as.numeric(g$belongings) + 1)[[x]]
                                                tmp$chip <- names(belongings)[[x]]
                                                tmp$region <- chipID_region[tmp$chip]
                                                tmp$VertexNum <- length(V(tmp))
                                                tmp$networkType <- tmpName
                                                tmp
                                                
                                              })
                                              
                                            })%>% purrr::reduce(c)

regionUnitStat.olg <- RegionGraphlistIncludedOLG_NN %>% purrr::map_dfr(function(x){
  data.frame(chip = x$chip,network = x$networkType, region = chipID_region[x$chip])
}) %>% group_by(chip,network,region) %>% summarise(freq = n()) %>% pivot_wider(names_from = "network",values_from = "freq", values_fill = 0) %>% select(-chip) %>% group_by(region) %>% summarise_all(mean) %>%  column_to_rownames("region")

# fig5f-1--------

ggheatmap::ggheatmap(as.matrix(regionUnitStat.olg) %>% t ,scale = "column",color =  MetBrewer::met.brewer("Demuth", 20) %>% rev) %>% ggheatmap::ggheatmap_theme(1,theme = list(theme(axis.text.x = element_text(angle = 90,face = "bold",size = 10,hjust = 1))))  %>%  ggsave(filename = str_c("./","olg_opc_NN",".pdf",sep = ""), height = 6,width = 8,dpi = 100)


order.olg <- regionUnitStat.olg %>% colMeans() %>% sort(decreasing = T) %>% names

regionUnitStat.olg.long <- regionUnitStat.olg %>% rownames_to_column("region") %>% pivot_longer(!region,names_to = "type",values_to = "freq") %>% mutate(type = factor(type,levels = order.olg))

# fig5f-2--------

{ggplot(regionUnitStat.olg.long,aes(x = type,y = freq)) +
    geom_line(mapping = aes( group = region,color = region)) +
    stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge") +
    ggprism::theme_prism() + theme(axis.text.x = element_text(angle = 90,hjust = 1)) + geom_smooth(mapping = aes(group=-1))
} %>% ggsave(filename = str_c("./","olg_opc_NN_line",".pdf",sep = ""), height = 8,width = 8,dpi = 100)





#  all network  ---------
all.net <- qs::qread(str_c("../",project,".qs"))

all.net.txt <- readLines(str_c("../",project))

head <- ( 1:length(all.net.txt))[str_detect(all.net.txt,"^t")]
end <- c(head[-c(1)]-2, length(all.net.txt))

qsFiles <- list.files("../../cellbin/filterCellbin/",".qs",full.names = T) %>% sort %>% dput

k =1
all.chipData <-  lapply(1:length(qsFiles),FUN = function(k){

  file_path = qsFiles[[k]]
  clmeta <- qs::qread(file_path)
  print(timestamp())
  tmp <- list(clmeta)
  names(tmp) <-  str_extract(basename(file_path), pattern = "[A-Z0-9]{1,}") %>% str_remove(pattern = "_")
  tmp
  }) %>% purrr::reduce(c)

# add the pathway -----------------

ccNet <- read.csv("../cell-cell_communications.all.csv",row.names = 1)

#fig5g---------

# ast
graphlistIncludedAST <- graphList %>% keep(~length(V(.))==2)  %>% keep(~any(str_detect(names(V(.)),"AST"))) %>%  keep(~any(str_detect(names(V(.)),"INH|EXC")))

ccc.ast <- parallel::mclapply(graphlistIncludedAST,mc.cores = min(length(graphlistIncludedAST),20),function(subnet) {
  
  tmpGraph <- subnet
  tmpDF <- tmpGraph %>% igraph::as_data_frame() 
  
  tmpDF.Left <- tmpDF %>% left_join(ccNet,by = c("from" = "source", "to" = "target"))
  tmpDF.Right <- tmpDF %>% left_join(ccNet,by = c("to" = "source", "from" = "target"))
  
  tmpCC <- rbind(tmpDF.Left,tmpDF.Right)
}) %>% purrr::reduce(rbind) %>% distinct()


ccc.ast$from_type <- subclass_class[ccc.ast$from]
ccc.ast$to_type <- subclass_class[ccc.ast$to]

cellpair.ast <- ccc.ast %>% select(from,to) %>% distinct()

# mgc
graphlistIncludedMICRO <- graphList %>% keep(~length(V(.))==2)  %>% keep(~any(str_detect(names(V(.)),"MICRO"))) %>%  keep(~any(str_detect(names(V(.)),"INH|EXC")))


ccc.micro <- parallel::mclapply(graphlistIncludedMICRO,mc.cores = min(length(graphlistIncludedMICRO),20),function(subnet) {

  tmpGraph <- subnet
  tmpDF <- tmpGraph %>% igraph::as_data_frame() 

  tmpDF.Left <- tmpDF %>% left_join(ccNet,by = c("from" = "source", "to" = "target"))
  tmpDF.Right <- tmpDF %>% left_join(ccNet,by = c("to" = "source", "from" = "target"))

  tmpCC <- rbind(tmpDF.Left,tmpDF.Right)
}) %>% purrr::reduce(rbind) %>% distinct()


ccc.micro$from_type <- subclass_class[ccc.micro$from]
ccc.micro$to_type <- subclass_class[ccc.micro$to]

cellpair.mgc <- ccc.micro %>% select(from,to) %>% distinct()

interaction.ast <-c("NRXN1 - LRRTM4", "NRG3 - ERBB4", "NRXN1 - NLGN1")

interaction.mgc <- c(
  "Glu-(SLC17A6+GLS) - GRIA2",
  "Glu-(SLC17A6+GLS) - GRIK2",
  "Glu-(SLC17A6+GLS2) - GRIA2",
  "Glu-(SLC17A6+GLS2) - GRIK2",

  "Glu-(SLC1A3+GLS) - GRM2",
  "Glu-(SLC1A3+GLS) - GRM3",
  "IL34 - CSF1R",
  "LRFN4 - PTPRD",

  "NRG4 - ERBB4",

  "PROS1 - MERTK",
  "SEMA4D - PLXNB2",

  "SLITRK1 - PTPRD",
  "SLITRK3 - PTPRD",

  "TGFB2 - (TGFBR1+TGFBR2)")

cellpair <- rbind(cellpair.ast,cellpair.mgc) %>% distinct()


tmpDF.Left <- cellpair %>% left_join(ccNet,by = c("from" = "source", "to" = "target"))
tmpDF.Right <- cellpair %>% left_join(ccNet,by = c("to" = "source", "from" = "target"))

tmpCC <- rbind(tmpDF.Left,tmpDF.Right)

tmpCC %<>% unite(col = "cellPair",from,to,remove = F,sep = "-")

tmpCC %<>% filter(interaction_name_2 %in% c(interaction.mgc,interaction.ast))

tmpCC$cellPair %<>% factor(levels = c("AST-L2-L3 IT LINC00507", "AST-L4-L5 IT RORB", "AST-SST", "AST-VIP","MICRO-ET", "MICRO-NP", "MICRO-L2-L3 IT LINC00507", "MICRO-L3-L4 IT RORB","MICRO-L4-L5 IT RORB", "MICRO-L6 CAR3", "MICRO-L6 IT", "MICRO-L6b","MICRO-LAMP5", "MICRO-NDNF", "MICRO-PAX6", "MICRO-PVALB","MICRO-VIP"))

{ggplot(tmpCC %>% filter(str_detect(cellPair,"NP|ET|CAR3",negate = T)),aes(x = cellPair, y = interaction_name_2,color = prob, size  = prob)) + geom_point() + ggprism::theme_prism() +MetBrewer::scale_color_met_c("Paquin",direction = -1) + theme(axis.text.x  = element_text(angle = 90, hjust = 1)) + coord_flip()
} %>% ggsave(filename = str_c("ast_mgc_LRpair.pdf"),height = 7, width =8)


# fig5h -------------

graphlistIncludedAST_NN <- graphList %>% keep(~length(V(.))==2)  %>% keep(~any(str_detect(names(V(.)),"AST"))) %>%  keep(~any(str_detect(subclass_class[names(V(.))],"INH|EXC")))

RegionGraphlistIncludedAST_NN <- purrr::map(graphlistIncludedAST_NN,
                                   function(g){
                                       belongings <- all.net[ as.numeric(g$belongings) +1] # 子网的编号要记得+1，这里面是从0开始计数的，上一级是all.net，all.net的编号芯片+对应meta里的细胞顺序
                                       tmpName <- generate_graph_name(g)
                                       belongings <- purrr::map(1:length(belongings),function(x) {
                                         tmp = belongings[[x]]
                                         tmp$org_id = (as.numeric(g$belongings) + 1)[[x]]
                                         tmp$chip <- names(belongings)[[x]]
                                         tmp$region <- chipID_region[tmp$chip]
                                         tmp$VertexNum <- length(V(tmp))
                                         tmp$networkType <- tmpName
                                         tmp

                                   })

                                       })%>% purrr::reduce(c)

regionUnitStat.ast <- RegionGraphlistIncludedAST_NN %>% purrr::map_dfr(function(x){
  data.frame(chip = x$chip,network = x$networkType, region = chipID_region[x$chip])
}) %>% group_by(chip,network,region) %>% summarise(freq = n()) %>% pivot_wider(names_from = "network",values_from = "freq", values_fill = 0) %>% select(-chip) %>% group_by(region) %>% summarise_all(mean) %>%  column_to_rownames("region")

# fig5h-1---------

ggheatmap::ggheatmap(as.matrix(regionUnitStat.ast),cluster_rows = T,show_cluster_rows = T,cluster_cols = F,scale = "row",color =  MetBrewer::met.brewer("Demuth", 20) %>% rev) %>% ggheatmap::ggheatmap_theme(1,theme = list(theme(axis.text.x = element_text(angle = 90,face = "bold",size = 10,hjust = 1))))  %>%  ggsave(filename = str_c("./","ast_NN",".pdf",sep = ""), height = 6,width = 8,dpi = 100)

order.ast <- regionUnitStat.ast %>% colMeans() %>% sort(decreasing = T) %>% names

regionUnitStat.ast.long <- regionUnitStat.ast %>% rownames_to_column("region") %>% pivot_longer(!region,names_to = "type",values_to = "freq") %>% mutate(type = factor(type,levels = order.ast))

# fig5h-2---------
{ggplot(regionUnitStat.ast.long,aes(x = type,y = freq)) +
    geom_line(mapping = aes( group = region,color = region)) +
    stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge") +
    ggprism::theme_prism() + theme(axis.text.x = element_text(angle = 90,hjust = 1)) + geom_smooth(mapping = aes(group=-1))} %>% ggsave(filename = str_c("./","ast_NN","_line.pdf",sep = ""), height = 6,width = 8,dpi = 100)


# fig5i --------

graphlistIncludedMGC_NN <- graphList %>% keep(~length(V(.))==2)  %>% keep(~any(str_detect(names(V(.)),"MICRO"))) %>%  keep(~any(str_detect(subclass_class[names(V(.))],"INH|EXC")))

RegionGraphlistIncludedMGC_NN <- purrr::map(graphlistIncludedMGC_NN,
                                            function(g){
                                              belongings <- all.net[ as.numeric(g$belongings) +1] 
                                              # tmpName <- generate_graph_name(g)
                                              tmpName <- str_c("MICRO","-",V(g) %>% names() %>% keep(~.!="MICRO"))
                                              belongings <- purrr::map(1:length(belongings),function(x) {
                                                tmp = belongings[[x]]
                                                tmp$org_id = (as.numeric(g$belongings) + 1)[[x]]
                                                tmp$chip <- names(belongings)[[x]]
                                                tmp$region <- chipID_region[tmp$chip]
                                                tmp$VertexNum <- length(V(tmp))
                                                tmp$networkType <- tmpName
                                                tmp

                                              })

                                            })%>% purrr::reduce(c)

regionUnitStat.mgc <- RegionGraphlistIncludedMGC_NN %>% purrr::map_dfr(function(x){
  data.frame(chip = x$chip,network = x$networkType, region = chipID_region[x$chip])
}) %>% group_by(chip,network,region) %>% summarise(freq = n()) %>% pivot_wider(names_from = "network",values_from = "freq", values_fill = 0) %>% select(-chip) %>% group_by(region) %>% summarise_all(mean) %>%  column_to_rownames("region")

# fig5i-1--------

ggheatmap::ggheatmap(as.matrix(regionUnitStat.mgc) %>% t ,scale = "row",color =  MetBrewer::met.brewer("Demuth", 20) %>% rev) %>% ggheatmap::ggheatmap_theme(1,theme = list(theme(axis.text.x = element_text(angle = 90,face = "bold",size = 10,hjust = 1))))  %>%  ggsave(filename = str_c("./","mgc_NN",".pdf",sep = ""), height = 6,width = 8,dpi = 100)

# fig5i-2--------

order.mgc <- regionUnitStat.mgc %>% colMeans() %>% sort(decreasing = T) %>% names

regionUnitStat.mgc.long <- regionUnitStat.mgc %>% rownames_to_column("region") %>% pivot_longer(!region,names_to = "type",values_to = "freq") %>% mutate(type = factor(type,levels = order.mgc))

{ggplot(regionUnitStat.mgc.long,aes(x = type,y = freq)) +
  geom_line(mapping = aes( group = region,color = region)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge") +
ggprism::theme_prism() + theme(axis.text.x = element_text(angle = 90,hjust = 1)) + geom_smooth(mapping = aes(group=-1))} %>% ggsave(filename = str_c("./","mgc_NN","_line.pdf",sep = ""), height = 6,width = 8,dpi = 100)


# fig5j -----------

graphlistIncludedMGCAST_NN <- graphList %>% keep(~length(V(.))==2)  %>% keep(~any(str_detect(names(V(.)),"MICRO|AST"))) %>%  keep(~any(str_detect(subclass_class[names(V(.))],"INH|EXC")))

RegionGraphlistIncludedMGCAST_NN <- purrr::map(graphlistIncludedMGCAST_NN,
                                            function(g){
                                              belongings <- all.net[ as.numeric(g$belongings) +1] 
                                              # tmpName <- generate_graph_name(g)
                                              tmpName <- str_c(ifelse(V(g) %>% names() %>% str_detect("MICRO") %>% any,"MICRO","AST"),"-",V(g) %>% names() %>% keep(~ . != "MICRO" & . != "AST") %>% {extract(subclass_class,.)})
                                              belongings <- purrr::map(1:length(belongings),function(x) {
                                                tmp = belongings[[x]]
                                                tmp$org_id = (as.numeric(g$belongings) + 1)[[x]]
                                                tmp$chip <- names(belongings)[[x]]
                                                tmp$region <- chipID_region[tmp$chip]
                                                tmp$VertexNum <- length(V(tmp))
                                                tmp$networkType <- tmpName
                                                tmp

                                              })

                                            })%>% purrr::reduce(c)


regionUnitStat.mgcast <- RegionGraphlistIncludedMGCAST_NN %>% purrr::map_dfr(function(x){
  data.frame(chip = x$chip,network = x$networkType, region = chipID_region[x$chip])
}) %>% group_by(chip,network,region) %>% summarise(freq = n()) %>% pivot_wider(names_from = "network",values_from = "freq", values_fill = 0) %>% select(-chip) %>% group_by(region) %>% summarise_all(mean) %>%  column_to_rownames("region")

ggheatmap::ggheatmap(as.matrix(regionUnitStat.mgcast) ,scale = "column",color =  MetBrewer::met.brewer("Demuth", 20) %>% rev) %>% ggheatmap::ggheatmap_theme(1,theme = list(theme(axis.text.x = element_text(angle = 90,face = "bold",size = 10,hjust = 1))))  %>%  ggsave(filename = str_c("./","mgc_ast_NN",".pdf",sep = ""), height = 4,width = 5,dpi = 100)




