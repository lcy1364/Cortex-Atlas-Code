 library(qs)
library(magrittr)
library(dplyr)
library(Seurat)
library(parallel)
library(stringr)
library(tidyverse)
library(data.table)
library(RANN)
library(igraph)
library(sf)
# library(spatstat)

setwd("~/DATA/BRAIN/STEREO/frequentGraph")

qsFiles <- list.files("../cellbin/filterCellbin/",".qs",full.names = T) %>% sort %>% dput

#应该把所有的放在一起算

# 获取命令行参数
args <- commandArgs(trailingOnly = TRUE)

# 打印参数
print(args)

# =============  GSPAN v2.40 - STATS =============
#   Number of graph in the input database: 39967
# Frequent subgraph count : 5026
# Total time ~ 19 s
# Minsup : 0 graphs
# Maximum memory usage : 4598.066482543945 mb
# ===================================================

#=============  GSPAN v2.40 - STATS =============
# Number of graph in the input database: 160648
# Frequent subgraph count : 4768
# Total time ~ 8 s
# Minsup : 0 graphs
# Maximum memory usage : 2515.8689041137695 mb
# ===================================================
dis_soma = 15
dis_nn = 15
support_rate = 0
filename = str_c("soma",dis_soma,"nn",dis_nn,"minup",support_rate,".network")

# 细胞类型要进行数字编码不能使用名称
ctCodes <- c("AST", "ENDO", "ET", "L2-L3 IT LINC00507", "L3-L4 IT RORB", # 1-5
             "L4-L5 IT RORB", "L6 CAR3", "L6 CT", "L6 IT", "L6b",  # 6-10
             "LAMP5", "MICRO", "NDNF", "NP", "OLIGO", # 11-15
             "OPC", "PAX6", "PVALB", "SST", "VIP", "VLMC") # 16-21

print(timestamp())
k = 1
# all.subgraph <- list()
all.subgraph <- mclapply(1:length(qsFiles),mc.cores = min( length(qsFiles),10),FUN = function(k){
# for (k in 28:length(qsFiles)) {

  file_path = qsFiles[[k]]
  clmeta <- qread(file_path)
  print(k)
  print(timestamp())
  # ----------------------------------------- cut cellular surrounding region ---------------------
  ChipId <-
    str_extract(basename(file_path), pattern = "[A-Z0-9]{1,}") %>% str_remove(pattern = "_")

  # clmeta <- seu@meta.data

  colnames(clmeta) <- make.unique(colnames(clmeta))


  clmeta$ct_code <- clmeta$first_type %>% as.character() %>% factor(levels = ctCodes) %>% as.numeric()



  clmeta$id <-1:nrow(clmeta)
  # clmeta$major <- subclass_class[clmeta$first_type]


  nn.res.non <- clmeta[c("x","y")]  %>% nn2(searchtype = "radius",radius = dis_soma *2,k = 20)#15um
  nn.idx.non <- nn.res.non$nn.idx
  nn.idx.non <- nn.idx.non[rowSums(nn.idx.non>0) > 1,]

  # 横幅变关联表
  network <- apply(nn.idx.non,MARGIN = 1, FUN = function(x){
    neighborS <- x[x>x[[1]]]
    if(length(neighborS) > 0) {
      data.frame(first = x[[1]], second =neighborS)
    }else{NULL}

  }) %>% purrr::reduce(rbind)


  # 删除重复的接触并且划分子网络，编号
  soma_net <- graph_from_data_frame(network,directed = F)

  # 去除重复的edge
  soma_net %<>% simplify()

  comps <- components(soma_net)

  subIDs <- comps$membership %>% unique %>% sort
  subgraphs <- list()
  print(length(subIDs))
  for (j in 1:length(subIDs)) {
    if(j%%1000 == 0) print(j)
    x = subIDs[j]
    tmp <- induced_subgraph(soma_net, which(comps$membership == x)) #%>% as_data_frame( what = "edges")
    # print(j)
    # print(clmeta$ct_code[V(tmp) %>% names() %>% as.numeric()])
    V(tmp)$type <- clmeta$ct_code[V(tmp) %>% names() %>% as.numeric()] %>% as.character()
    E(tmp)$type <- "1"
    subgraphs <- append(subgraphs,list(tmp))
    names(subgraphs)[[length(subgraphs)]] <- ChipId

  }
  subgraphs
  # all.subgraph[[k]] <- subgraphs
})
print(timestamp())
print(length(all.subgraph))

all.subgraph %<>% reduce(c)


qsave(all.subgraph,str_c(filename,".qs"))
# print(timestamp())
# all.subgraph <- qread(str_c(filename,".qs"))

mclapply(1:length(all.subgraph), mc.cores = min(ceiling(length(all.subgraph)/500),50) ,function(j){
  tmp = all.subgraph[[j]]
  netDF <- tmp %>% igraph::as_data_frame( what = "both")
  # cat(,file = filename, sep = ,append = T)
  paste0(paste0("t\ #\ ", format(j-1, scientific = FALSE) ,"\n"),
        apply(netDF$vertices, 1, function(v){
          paste0(paste0(c("v",paste0(unlist(v),collapse = " ")),collapse = " "),"\n")
        }) %>% paste0(collapse = ""),
        apply(netDF$edges, 1, function(v){
          paste0(paste0(c("e",paste0(unlist(v),collapse = " ")),collapse = " "),"\n")
        }) %>% paste0(collapse = ""),
        "\n")

}) %>% str_c(collapse = "") %>% cat(file = filename)

print(str_c("java -jar ~/DATA/BRAIN/STEREO/frequentGraph/spmf.jar run GSPAN ", filename, " ",str_c(filename,".fsm "),support_rate," 3 false false true"))

system(str_c("java -jar ~/DATA/BRAIN/STEREO/frequentGraph/spmf.jar run GSPAN ", filename, " ",str_c(filename,".fsm "),support_rate," 3 false false true"))
# system(str_c("java -jar ~/DATA/BRAIN/STEREO/frequentGraph/spmf.jar run GSPAN ", filename, " ",str_c(filename,".fsm "),0.001,"  3 false false true"))

# high frequent network ----------------------------
project = "../soma15nn15.network"

resTxt <- readLines(str_c("../",project,"0.fsm"))
# resTxt <- readLines("./all_15.network.fsm")

head <- ( 1:length(resTxt))[str_detect(resTxt,"^t")]

end <- ( 1:length(resTxt))[str_detect(resTxt,"^x")]
ctCodes <- c("AST", "ENDO", "ET", "L2-L3 IT LINC00507", "L3-L4 IT RORB", # 1-5
             "L4-L5 IT RORB", "L6 CAR3", "L6 CT", "L6 IT", "L6b",  # 6-10
             "LAMP5", "MICRO", "NDNF", "NP", "OLIGO", # 11-15
             "OPC", "PAX6", "PVALB", "SST", "VIP", "VLMC") # 16-21


graphList <-  parallel::mclapply(1:length(head), mc.cores = 20, function(i){
  tmp <-resTxt[head[[i]]:end[[i]]]
  
  Edg <- tmp %>% str_subset("^e") %>% map_dfr(.f = function(x){
    str_split_fixed(x," ",n = 4)[2:3] %>% as.numeric()   %>% { data.frame(first= .[[1]],second = .[[2]])}})
  
  Vertice <- tmp %>% str_subset("^v") %>% map_dfr(.f = function(x){
    str_split_fixed(x," ",n = 3)[2:3] %>% as.numeric() %>% { data.frame( id = .[[1]],name = .[[2]])}})
  
  g <- igraph::graph_from_data_frame(Edg,vertices = Vertice,directed = F)
  
  V(g)$name <- ctCodes[V(g) %>% names() %>% as.numeric()] %>% as.character()
  
  g$belongings <- tmp %>% str_subset("^x") %>% str_remove("x ") %>% str_split(" ",simplify = T)
  g$support <- length(g$belongings)
  
  g
})

graphList <- graphList[graphList %>% map(~.$support) %>% unlist %>% order %>% rev]

# Non inh exc

graphList <- lapply(1:length(graphList), function(i){
  tmp <- list(graphList[[i]])
  names(tmp) <-subclass_class[V(graphList[[i]])$name] %>% unique %>% sort %>% paste0(collapse = "_")
  tmp
}) %>% purrr::reduce(c)

qs::qsave(graphList,"graphList.qs")
