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

qsFiles <- list.files(".","merge.qs",full.names = T)
x = qsFiles[[1]]

datasets <- parallel::mclapply(qsFiles,mc.cores = length(qsFiles),FUN = function(x){
  tmp = qs::qread(x,)
  tmp  %<>% DietSeurat()})

i = 1
trees <- parallel::mclapply(1:50,mc.cores = 25,FUN = function(i){
  
  seu_sampled <- parallel::mclapply(datasets,mc.cores = length(datasets),FUN = function(x){
    set.seed(i)
    tmp = x
    sampled_cells <- tmp@meta.data %>% rownames_to_column("cellID") %>% 
      group_by(cluster, region) %>%  # Group by cluster and area
      sample_n(size = min(100, n()), replace = FALSE) %>%  # Sample up to 100 cells without replacement
      ungroup() %$% cellID
    
    tmp[,sampled_cells]
  })  %>% {merge(.[[1]],.[-c(1)])}
  
  data <- seu_sampled@assays$RNA@data
  
  data <- data[rowSums(data)>0,]
  
  vgs_all <- scrattch.hicat::find_vg(data,verbose = T)

  # how many genes included
  
  geneTree = list()
  for(nGenes in c(2500,3000,3500,4000)){
    
    vgs <- vgs_all %>% slice_max(n = nGenes, order_by = loess.z)
    
    cl <- setNames(seu_sampled@meta.data$cluster,
                   rownames(seu_sampled@meta.data))
    
    cl.dat <- scrattch.hicat::get_cl_medians(data[vgs$gene, ], cl = cl)
    
    tree <- scrattch.hicat::build_dend(cl.dat, ncores = 5L, nboot = 500)
    
    tree_edlein <- tree$pvclust.result$hclust %>% ape::as.phylo() %>% as_tibble()
    
    tree_edlein$subclass <- tree_edlein$label %>% str_remove("_[0-9]*$")
    
    
    tree_plot_edlein <- ggtree(tree_edlein %>% as.treedata()) + geom_tiplab(
      mapping = aes(color = subclass),
      angle = 90,
      hjust = 1,
      vjust = 0
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
    ) + cowplot::theme_cowplot() + theme(legend.position = "none",
                                         axis.text.x = element_text(size = 0.5)) + xlim(c(0.5, -0.5))
    
    ggsave(
      str_c("5_dend_", i, "_", nGenes, ".png"),
      tree_plot_edlein,
      width = 20,
      height = 5
    )
    
    geneTree[[str_c("gene",nGenes)]] <- list(tree = tree_plot_edlein, vgs = vgs)
  }

    qs::qsave(list(cells = colnames(seu_sampled),vgs_all,geneTree),str_c("5_sampled_",i,".qs"))
  
})



# 
# 
# sampled_data <- mclapply(unique(subclasses),mc.cores = length(unique(subclasses)),function(subclass){
#   seu <- subclass_data[subclasses == subclass] %>% purrr::reduce(merge)
#   
#   # seu@meta.data %>% rownames_to_column("cellID") %>% group_by(region,cluster) %>% 
#   
#   sampled_cells <- seu@meta.data %>% rownames_to_column("cellID") %>% 
#     group_by(cluster, region) %>%  # Group by cluster and area
#     sample_n(size = min(500, n()), replace = FALSE) %>%  # Sample up to 100 cells without replacement
#     ungroup() %$% cellID
#   
#   seu[,sampled_cells]
#   
#   # qs::qsave(seu,str_c(str_replace(seu$subclass[[1]],"/","_"),"_subclass.qs"),nthreads = 10,preset = "fast")
#   
#   # sample 10% for annotating datasets globally
#   
#   # seu[,sample(colnames(seu),ceiling(length(colnames(seu))*0.2),replace = F)]
# })  
# 
# seu_edlein <- merge(sampled_data[[1]],sampled_data[-c(1)])
# 
# seu_edlein$batch = "edlein"
# 
# seu_edlein %>% qs::qsave("SnRNA_0.2DownSampled.qs",nthreads = 20,preset = "fast")
# 
# 
# # 获取细胞的簇标签
# Idents(seu_edlein) <- "cluster"
# clusters <- Idents(seu_edlein)
# 
# # 获取高度变异基因（HVGs），并选择前3000个
# hvg_genes <- FindVariableFeatures(seu_edlein,nfeatures = 3000) %>% VariableFeatures()
# 
# # 获取基因表达数据中高度变异基因的表达矩阵
# expression_data <- GetAssayData(seu_edlein, assay = "RNA", slot = "data")[hvg_genes, ]
# 
# # 定义函数来计算每个簇的基因表达中位值
# calculate_median_expression <- function(cluster) {
#   # 获取属于当前簇的细胞
#   cluster_cells <- WhichCells(seu_edlein, idents = cluster)
#   
#   # 提取当前簇的高度变异基因表达矩阵
#   cluster_expression <- expression_data[, cluster_cells]
#   
#   # 计算当前簇中每个基因的中位表达值
#   cluster_median_expression <- apply(cluster_expression, 1, median)
#   
#   return(cluster_median_expression)
# }
# 
# # 使用 mclapply 来并行计算每个簇的中位表达值
# median_expression_by_cluster <- mclapply(levels(clusters), calculate_median_expression, mc.cores = 22)
# 
# # 将结果转化为数据框，行是基因，列是簇
# median_expression_df <- do.call(cbind, median_expression_by_cluster)
# 
# # 为列名添加簇标签
# colnames(median_expression_df) <- levels(clusters)
# 
# # 打印结果
# head(median_expression_df)
# 
# pvclust_tree_edlein <- median_expression_df %>% pvclust::pvclust(parallel = 10L) 
# 
# tree_edlein <- pvclust_tree_edlein$hclust %>% ape::as.phylo() %>% as_tibble()
# 
# tree_edlein$subclass <- tree_edlein$label %>% str_remove("_[0-9]*$")
# 
# tree_plot_edlein <- ggtree(tree_edlein %>% as.treedata()) + geom_tiplab(mapping = aes(color = subclass),angle = 90, hjust = 1, vjust = 0) + ggtree::layout_dendrogram() + theme(legend.position = "none") + scale_color_manual(values = c(AST = "#665C47", ENDO = "#604B47", ET = "#CEC823", CHANDELIER = "#E25691", 
#                                                                                                                                                                                                                                           `L2-L3 IT LINC00507` = "#07D8D8", `L3-L4 IT RORB` = "#09B2B2", 
#                                                                                                                                                                                                                                           `L4-L5 IT RORB` = "#69B199", `L6 CAR3` = "#898325", `L6 CT` = "#2D8CB8", 
#                                                                                                                                                                                                                                           `L6 IT` = "#19E3BE", L6b = "#7944AA", LAMP5 = "#D96C07", MICRO = "#513577", 
#                                                                                                                                                                                                                                           NDNF = "#f2798d", NP = "#93C43B", OLIGO = "#5E8A79", OPC = "#2E3E39",                                                                                                                                                                                                                                           PAX6 = "#E96DF2", PVALB = "#FF2D4E", SST = "#F2B003", VIP = "#9C28CC", 
#                                                                                                                                                                                                                                           VLMC = "#697255"))
# 
# tree_plot_edlein
# 
# 
# 
# 
# aver_exp_edlein <- AverageExpression(seu_edlein,group.by = "cluster",assays = "RNA",features = FindVariableFeatures(seu_edlein,nfeatures = 2000) %>% VariableFeatures()) 
# 
# pvclust_tree_edlein <- aver_exp_edlein[[1]] %>% pvclust::pvclust(parallel = 10L) 
# 
# tree_edlein <- pvclust_tree_edlein$hclust %>% ape::as.phylo() %>% as_tibble()
# 
# tree_edlein$subclass <- tree_edlein$label %>% str_remove("_[0-9]*$")
# 
# tree_plot_edlein <- ggtree(tree_edlein %>% as.treedata()) + geom_tiplab(mapping = aes(color = subclass),angle = 90, hjust = 1, vjust = 0) + ggtree::layout_dendrogram() + theme(legend.position = "none") + scale_color_manual(values = c(AST = "#665C47", ENDO = "#604B47", ET = "#CEC823", CHANDELIER = "#E25691", 
#                                                                                                                                                                                                                                           `L2-L3 IT LINC00507` = "#07D8D8", `L3-L4 IT RORB` = "#09B2B2", 
#                                                                                                                                                                                                                                           `L4-L5 IT RORB` = "#69B199", `L6 CAR3` = "#898325", `L6 CT` = "#2D8CB8", 
#                                                                                                                                                                                                                                           `L6 IT` = "#19E3BE", L6b = "#7944AA", LAMP5 = "#D96C07", MICRO = "#513577", 
#                                                                                                                                                                                                                                           NDNF = "#f2798d", NP = "#93C43B", OLIGO = "#5E8A79", OPC = "#2E3E39",                                                                                                                                                                                                                                           PAX6 = "#E96DF2", PVALB = "#FF2D4E", SST = "#F2B003", VIP = "#9C28CC", 
#                                                                                                                                                                                                                                           VLMC = "#697255"))
# 
# tree_plot_edlein
# 
# 
# aver_exp_edlein <- AverageExpression(seu_edlein,group.by = "subclass",assays = "RNA",features = FindVariableFeatures(seu_edlein,nfeatures = 3000) %>% VariableFeatures()) 
# 
# pvclust_tree_edlein <- aver_exp_edlein[[1]] %>% pvclust::pvclust(parallel = 10L) 
# 
# tree_edlein <- pvclust_tree_edlein$hclust %>% ape::as.phylo() %>% as_tibble()
# 
# tree_edlein$subclass <- tree_edlein$label %>% str_remove("_[0-9]*$")
# 
# tree_plot_edlein <- ggtree(tree_edlein %>% as.treedata()) + geom_tiplab(mapping = aes(color = subclass),angle = 90, hjust = 1, vjust = 0) + ggtree::layout_dendrogram() + theme(legend.position = "none") + scale_color_manual(values = c(AST = "#665C47", ENDO = "#604B47", ET = "#CEC823", CHANDELIER = "#E25691", 
#                                                                                                                                                                                                                                           `L2-L3 IT LINC00507` = "#07D8D8", `L3-L4 IT RORB` = "#09B2B2", 
#                                                                                                                                                                                                                                           `L4-L5 IT RORB` = "#69B199", `L6 CAR3` = "#898325", `L6 CT` = "#2D8CB8", 
#                                                                                                                                                                                                                                           `L6 IT` = "#19E3BE", L6b = "#7944AA", LAMP5 = "#D96C07", MICRO = "#513577", 
#                                                                                                                                                                                                                                           NDNF = "#f2798d", NP = "#93C43B", OLIGO = "#5E8A79", OPC = "#2E3E39",                                                                                                                                                                                                                                           PAX6 = "#E96DF2", PVALB = "#FF2D4E", SST = "#F2B003", VIP = "#9C28CC", 
#                                                                                                                                                                                                                                           VLMC = "#697255"))
# 
# tree_plot_edlein
# 
