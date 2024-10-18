setwd("~/data/STEREO/GEF/bins")
library(qs)

cortexMeta <- read.delim("~/data/STEREO/AnalysisPlot/cortex") 
chipRegion <- setNames(cortexMeta$region,cortexMeta$chip) 
bin100Files <- list.files(".","100.qs")
chip = cortexMeta$chip[[1]]
dataCollect <- list()
for ( chip in cortexMeta$chip){
  data <- bin100Files %>% str_subset(chip) %>% qread
  dataCollect[[chip]] <- data.frame(chip = chip,
             NumberOfbin100 = data %>% nrow,
             total_reads = data$nCount_Spatial %>% sum(),
             total_genes = ncol(data),
             median_reads =   data$nCount_Spatial %>% median(),
             mean_reads = data$nCount_Spatial %>% mean(),
             median_genes =   data$nFeature_Spatial %>% median(),
             mean_genes = data$nFeature_Spatial %>% mean()
)
}
qcdata <- dataCollect %>% purrr::reduce(rbind)
