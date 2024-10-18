#https://github.com/chanzuckerberg/single-cell-curation/blob/main/notebooks/curation_api/R/get_datasets_R.ipynb
setwd("~/cortex/SnRNA/2_codePreprocessingExternalData/")
library(readr)
library(httr)
library(stringr)
library(rjson)
library(tidyverse)
library(magrittr)


domain_name <- "cellxgene.cziscience.com"
site_url <- str_interp("https://${domain_name}")
api_url_base <- str_interp("https://api.${domain_name}")
collection_id <- "d17249d2-0e6e-4500-abb8-e6c93fa1ac6f"

collection_path <- str_interp("/curation/v1/collections/${collection_id}")
collection_url <- str_interp("${api_url_base}${collection_path}")
res <- GET(url=collection_url, `Content-Type`="application/json")
stop_for_status(res)
res_content <- content(res)
print(res_content)

alldatasets <- res_content$datasets
alldatasets %<>% keep(~str_detect(.$title,"Supercluster"))
 
parallel::mclapply(alldatasets,mc.cores = length(alldatasets),FUN = function(x){
  
  tmp <- x$assets %>% keep(~.$filetype == "H5AD")
  fileName = str_c(x$title,'.h5ad')
  system(str_c('axel -n 10 -o \'',fileName,"\' ",tmp[[1]]$url))
  })



