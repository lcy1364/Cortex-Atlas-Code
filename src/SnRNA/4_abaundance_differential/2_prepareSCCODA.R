library(tidyverse)
setwd("~/cortex/SnRNA/4_abaundance_differential/")
library(data.table)
library(magrittr)
# subclass neuron test
# cluster  test

wholeMeta <- read.csv("../3_mergingDatasets/SnRNA_Meta.csv")

region <- c("ACC", "AG", "DLPFC", "FPPFC", "ITG", "M1", "PoCG", "S1", "S1E", 
            "SMG", "SPL", "STG", "V1", "VLPFC")

sub_neu <- expand_grid( "./neuronsubclass.h5ad",
                        region, 
                        wholeMeta %>% filter(class != "non-neuronal") %$% subclass %>% unique %>% str_replace_all(" ","=")) %>% unite(col = "command",dplyr::everything(),sep = " ")

sub_non <- expand_grid( "./nonsubclass.h5ad",
                        region, 
                        wholeMeta %>% filter(class == "non-neuronal") %$% subclass %>% unique %>% str_replace_all(" ","=") ) %>% unite(col = "command",dplyr::everything(),sep = " ")

# cluster  test

cl_neu <- expand_grid("./neuroncluster.h5ad",
                      region,
                      wholeMeta %>% filter(class != "non-neuronal") %$% cluster %>% unique %>% str_replace_all(" ","=")) %>% unite(col = "command",dplyr::everything(),sep = " ")


cl_non <- expand_grid("./noncluster.h5ad",
                      region,
                      wholeMeta %>% filter(class == "non-neuronal") %$% cluster %>% unique %>% str_replace_all(" ","=") ) %>% unite(col = "command",dplyr::everything(),sep = " ")


list(sub_neu,sub_non,cl_neu,cl_non) %>% purrr::reduce(rbind) %>% unlist %>% cat(file = "cmd.txt",sep = "\n")
