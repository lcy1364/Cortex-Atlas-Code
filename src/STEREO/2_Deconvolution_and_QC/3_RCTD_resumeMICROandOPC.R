setwd("~/cortex/STEREO/2_Deconvolution_and_QC/")

chipList <- read.delim("~/DATA/data/STEREO/AnalysisPlot/cortex") %>% {setNames(nm = .$chip,.$region)}
chipList2 <- read.delim("~/DATA/data/STEREO2/sampleMeta") %>% filter(cortex_category == "neocortex") %>% {setNames(nm = .$chipID,.$abbr)}
chipList2[str_detect(chipList2,"DLPFC")] = "DLPFC"
chipList2 <- chipList2[chipList2 %in% chipList]
chipList  <- c(chipList,chipList2)


region_to_depth = c(0,0.5,NA,1,2,3,4,5,6,7)
names(region_to_depth) <- c("ARACHNOID", "ARACHNOID L1", "EMPTY", "L1", "L2", "L3", "L4", "L5", "L6",  "WM")

STEREO = list.files("../1_SnRNA_preprocessing/","qs",full.names = T) %>% keep(~str_detect(.,names(chipList) %>% str_c(collapse = "|")))

x = STEREO1[[1]]
OPC_MICRO <- mclapply(STEREO,mc.cores = 10,FUN = function(x){
  
  seu <- readRDS(x)
  chipID =  basename(x) %>% str_extract("^[A-Z0-9]*")
  seu <- seu[,seu$cellType %in% c("MICRO","OPC")]
  seu$cluster <- ifelse(seu$cellType == "MICRO","Micro/PVM_1","OPC_1")
  tmp = seu@meta.data[c("x","y","mainLayer","cellType","cluster")]
  colnames(tmp) <- c("x","y","domain","subclass","cluster")
  tmp$domain %<>% str_to_upper()
  tmp$chip = chipID
  tmp$region = chipList[chipID]
  tmp$depth = region_to_depth[tmp$domain]
  tmp
}) %>% purrr::reduce(rbind)
# "x"             "y"             "chip"     "region"        "domain"        "depth"         "subclass"      "cluster"


OPC_MICRO[is.na(OPC_MICRO) %>% rowSums() > 0,] %$% domain %>% unique

OPC_MICRO %<>% drop_na()
OPC_MICRO %>% write_csv("OPC_MICRO.csv")
