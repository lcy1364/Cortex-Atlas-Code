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

# fig2b heatmap ==========================================
project = "../soma15nn15.network"

resTxt_all <- readLines(str_c("../",project,"0.fsm"))

head <- ( 1:length(resTxt_all))[str_detect(resTxt_all,"^t")]
end <- ( 1:length(resTxt_all))[str_detect(resTxt_all,"^x")]
ctCodes <- c("AST", "ENDO", "ET", "L2-L3 IT LINC00507", "L3-L4 IT RORB", # 1-5
             "L4-L5 IT RORB", "L6 CAR3", "L6 CT", "L6 IT", "L6b",  # 6-10
             "LAMP5", "MICRO", "NDNF", "NP", "OLIGO", # 11-15
             "OPC", "PAX6", "PVALB", "SST", "VIP", "VLMC") # 16-21

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
