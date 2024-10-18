library(purrr)
library(magrittr)
library(tidyverse)
library(Seurat)
library(harmony)
library(ape)
library(uwot)
library(ggtree)
library(treeio)
library(ggtree)
library(treeio)
# library(future)
setwd("~/cortex/SnRNA/3_mergingDatasets/")

qsFiles <- list.files(".", "merge.qs", full.names = T)
x = qsFiles[[1]]

datasets <- parallel::mclapply(
  qsFiles,
  mc.cores = length(qsFiles),
  FUN = function(x) {
    tmp = qs::qread(x, )
    tmp  %<>% DietSeurat()
  }
)

seu_sampled <- parallel::mclapply(
  datasets,
  mc.cores = length(datasets),
  FUN = function(x) {
    set.seed(44)
    tmp = x
    sampled_cells <- tmp@meta.data %>% rownames_to_column("cellID") %>%
      group_by(cluster, region) %>%  # Group by cluster and area
      sample_n(size = min(100, n()), replace = FALSE) %>%  # Sample up to 100 cells without replacement
      ungroup() %$% cellID
    
    tmp[, sampled_cells]
  }
)  %>% {
  merge(.[[1]], .[-c(1)])
}

data <- seu_sampled@assays$RNA@data

data <- data[rowSums(data) > 0, ]

vgs_all <- scrattch.hicat::find_vg(data, verbose = T)

vgs <- vgs_all %>% slice_max(n = 3000, order_by = loess.z)

cl <- setNames(seu_sampled@meta.data$cluster,
               rownames(seu_sampled@meta.data))

cl.dat <- scrattch.hicat::get_cl_medians(data[vgs$gene, ], cl = cl)

tree <- scrattch.hicat::build_dend(cl.dat, ncores = 5L, nboot = 500)

tree <- tree$pvclust.result$hclust %>% ape::as.phylo() %>% as_tibble()

tree$subclass <- tree$label %>% str_remove("_[0-9]*$")

tree %<>% as.treedata()

treeplot <- ggtree(tree) + geom_tiplab(
  mapping = aes(color = subclass),
  angle = 90,
  hjust = 1,
  vjust = 0.5
) + ggtree::layout_dendrogram() + scale_color_manual(
  values = c(
    AST = "#665C47",
    ENDO = "#604B47",
    ET = "#CEC823",
    CHANDELIER = "#E25691",
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
) + cowplot::theme_nothing() + theme(legend.position = "none",
                                     axis.text.x = element_text(size = 0.5)) + xlim(c(0.3, -0.5)) 

treeplot + geom_nodelab(mapping = aes(label = node))
treeplot <- ggtree::flip(treeplot,node1 = 186, node2 = 185)
treeplot <- ggtree::rotate(treeplot,node = 172)
treeplot <- ggtree::rotate(treeplot,node = 173)
treeplot <- ggtree::rotate(treeplot,node = 250)
treeplot <- ggtree::rotate(treeplot,node = 180)
treeplot <- ggtree::rotate(treeplot,node = 248)
treeplot <- ggtree::rotate(treeplot,node = 236)

treeplot + geom_nodelab(mapping = aes(label = node))

ggsave(
  str_c("5_dend.png"),
  treeplot,
  width = 20,
  height = 3
)


saveRDS(treeplot,"../dend.RDS")

qs::qsave(list(seu_sampled,vgs_all,treeplot),"5_sampled_tree.qs")


