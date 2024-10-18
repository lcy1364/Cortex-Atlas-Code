options(bitmapType = "cairo")
library(Seurat)
library(scrattch.hicat)
library(dendextend)
library(tidyverse)
library(matrixStats)
library(Matrix)
library(magrittr)
library(RColorBrewer)
library(ranger)
library(ggheatmap)
library(patchwork)
library(ggtree)
library(ggcorrplot)
library(bioDist)
library(harmony)
library(rtracklayer)

geneid_name <- read.csv("../SnRNA/1_SnRNA_preprocessing/gene_kept.csv") %>% {setNames(.$gene_name,.$gene_id)}
setwd("~/cortex/figS1-6/")


#function change
{
  scale_size_2 <-
    function (name = waiver(),
              breaks = waiver(),
              labels = waiver(),
              limits = NULL,
              range = c(1, 6),
              trans = "identity",
              guide = "legend")
    {
      continuous_scale(
        "size",
        "area",
        area_pal(range),
        name = name,
        n.breaks = 4,
        breaks = breaks,
        labels = labels,
        limits = limits,
        trans = trans,
        guide = guide
      )
    }
  environment(scale_size_2) <- environment(ggplot2::scale_size)
  
}


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
  )
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
regionID <-
  c('A1',
    'A2',
    'A3',
    'A4',
    'A5',
    'A6',
    'A7',
    'A8',
    'A9',
    'A10',
    'A11',
    'A12',
    'A13',
    'A14')
names(region_color) <- region
names(region) <- regionID

qsFiles <- list.files("../SnRNA/3_mergingDatasets/", "merge.qs", full.names = T) %>% keep(~str_detect(.,"IT"))
x = qsFiles[[1]]
datasets <- parallel::mclapply(
  qsFiles,
  mc.cores = length(qsFiles),
  FUN = function(x) {
    tmp = qs::qread(x)
    set.seed(3)
    sampled_cells <- tmp@meta.data %>% rownames_to_column("cellID") %>%
      group_by(cluster, region) %>%  # Group by cluster and area
      sample_n(size = min(200, n()), replace = FALSE) %>%  # Sample up to 100 cells without replacement
      ungroup() %$% cellID
    
    tmp[, sampled_cells]
  }
)  

x = datasets[[1]]

results <- parallel::mclapply(datasets, mc.cores = length(datasets), function(x) {
  print(date())
  print(x)
  norm.dat <- x
  print(dim(norm.dat))
  cl.clean <- norm.dat$region
  
  tmp.med <- get_cl_medians(norm.dat@assays$RNA@data,
                            cl.clean)
  
  cl.cor <- cor(tmp.med)
  cl.cor.long <-
    cl.cor %>%
    as.data.frame() %>%
    rownames_to_column("sample") %>%
    pivot_longer(cols = !sample,
                 names_to = "to_sample",
                 values_to = "correlation")
  
  print(date())
  print("starting randomForest...")
  set.seed(123)
  norm.dat %<>% FindVariableFeatures()
  norm.dat %<>% RunPCA(npcs = 30, verbose = FALSE)
  norm.dat %<>% RunHarmony(group.by.vars = "region")
  norm.dat %<>% RunUMAP(dim = 1:30, reduction = "harmony")
  
  train_data <-
    norm.dat@reductions$pca@cell.embeddings %>% as.data.frame()
  
  train_data["label"] <- norm.dat$region %>% factor()
  
  train.forest <- ranger(label ~ .,
                               data = train_data,num.threads = 20)
  
  print(date())
  print("starting de")
  Idents(norm.dat) <- norm.dat$region
  
  DE <- parallel::mclapply(norm.dat$region %>% unique,mc.cores = 14,FUN = function(x){
    tmp <- FindMarkers(norm.dat,ident.1 = x,ident.2 = NULL,group.by = "region")
    tmp$cluster <- x
    tmp$gene <- rownames(tmp)
    tmp
  }) %>% purrr::reduce(rbind)
  
  # DE <- FindAllMarkers(norm.dat)
  
  
  subclass_result <-
    list(
      cl.cor.long = cl.cor.long,
      train.forest = train.forest,
      DE = DE,
      seu = norm.dat,
      subclass = x$subclass[[1]]
    )
  
  qs::qsave(subclass_result, str_c("figS2e", x$subclass[[1]], "_res.qs"))

  subclass_result
  
})

subclass_result = results[[1]]
plotresult <- parallel::mclapply(results, mc.cores = length(results), function(subclass_result) {
  
  train.forest <- subclass_result$train.forest
  tmp_plot <-
    train.forest$confusion %>% as.data.frame() %>%  group_by(true) %>% mutate(per = Freq/sum(Freq)) %>% 
    pivot_wider(id_cols = !Freq, names_from = predicted,values_from = per) %>% column_to_rownames("true") %>%
    as.matrix
  
  p_dend <- try({
    pvclust::pvclust(tmp_plot %>% t())
  })
  
  p_dend <- p_dend$hclust %>% ape::as.phylo() %>% ggtree()
  
  subclass_result$p_dend <- p_dend
  
  region_order <- get_taxa_name(p_dend)
  
  # dend
  plotresult <- list()
  
  plotresult[[subclass_result$subclass]]$dend <- p_dend + layout_dendrogram()
  
  # 0_umap
  plotresult[[subclass_result$subclass]]$umap <- DimPlot(
    subclass_result$seu,
    reduction = "umap",
    group.by = "region",
    shuffle = T,
    cols = region_color
  ) +
    theme(title = element_blank(), legend.position = "none") + labs(title = subclass_result$subclass)
  
  plotresult[[subclass_result$subclass]]$umap_splited <- DimPlot(
    subclass_result$seu,
    reduction = "umap",
    group.by = "region",
    split.by = "region",
    ncol = 3,
    shuffle = T,
    cols = region_color
  ) +
    theme(title = element_blank())
  
  # 1 expression similarity
  cl.cor.long <- subclass_result$cl.cor.long
  cl.cor.long$sample %<>% factor(levels = region_order)
  cl.cor.long$to_sample %<>% factor(levels = region_order)
  
  
  plotresult[[subclass_result$subclass]]$cl.cor <-
    ggplot(cl.cor.long,
           aes(x = sample,
               y = to_sample,
               fill = correlation)) +
    geom_tile(color = "black") +
    scale_fill_gradientn(colors = rev(brewer.pal(11, "Spectral"))) +
    labs(x = "", y = "", title = "Expression Similarity") +
    cowplot::theme_cowplot() %+replace% theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      legend.text = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      legend.title = element_text(vjust = 1),
      legend.title.align = 0
    )
  
  # 2 separatibility heatmap
  tmp_plot_long <-
    tmp_plot %>% t() %>% as.data.frame %>% rownames_to_column("region") %>% pivot_longer(cols = !region,
                                                                                         names_to = "region2",
                                                                                         values_to = "percentage")
  tmp_plot_long$region %<>% factor(levels = region_order)
  tmp_plot_long$region2 %<>% factor(levels = region_order)
  
  plotresult[[subclass_result$subclass]]$confusion_matrix <-
    ggplot(tmp_plot_long,
           aes(x = region,
               y = region2,
               fill = percentage)) +
    geom_tile(color = "black") +
    # scale_size( range = c(1, 6)*2)+
    scale_fill_gradientn(colors = rev(brewer.pal(11, "Spectral"))) +
    labs(x = "Predicted Region", y = "Original Region", title = "Consufion Matrix for Separatibility", ) +
    cowplot::theme_cowplot() %+replace% theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      legend.text = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      legend.title = element_text(vjust = 1),
      legend.title.align = 0
    )
  
  
  DE <- subclass_result$DE
  DE %>%
    filter(!is.na(geneid_name[gene]), avg_log2FC > 0.2) %>%
    group_by(cluster) %>%
    top_n(n = 3, wt = avg_log2FC) %>%
    ungroup() -> top10
  
  top10$cluster %<>% as.character() %>% factor(levels = region_order)
  top10 %>% arrange(cluster, desc(avg_log2FC)) %>% distinct(gene) %>% `$`("gene") -> top_genes
  
  plotresult[[subclass_result$subclass]]$deg_filter <- DotPlot(subclass_result$seu,
                                                      features = top_genes ,
                                                      group.by = "region") + theme(
                                                        axis.title.x = element_blank(),
                                                        axis.title.y = element_blank(),
                                                        axis.text.x = element_text(
                                                          angle = 90,
                                                          hjust = 1,
                                                          vjust = 0.5
                                                        )
                                                      ) +
    scale_color_gradientn(colors = rev(brewer.pal(11, "Spectral")), n.breaks = 4) +
    scale_size_2() +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "vertical",
      legend.box.just = "right"
    ) +
    guides(
      size = guide_legend("Percent Expressed", order = 1),
      color = guide_colorbar("Average Expression", order = 2)
    ) + coord_flip() + scale_x_discrete(labels = geneid_name[top_genes])
  plotresult[[subclass_result$subclass]]$deg_filter$data$id %<>% factor(levels = region_order)
  plotresult
}) 

plotresult <- plotresult[c("L2-L3 IT LINC00507",
                           "L3-L4 IT RORB",
                           "L4-L5 IT RORB",
                           "L6 IT")]
{
  S3index <-
    names(plotresult)[names(plotresult) %in% c("L2-L3 IT LINC00507",
                                               "L3-L4 IT RORB",
                                               "L4-L5 IT RORB",
                                               "L6 IT")]
  # names(plotresult[[1]])
  plot_id <- c(2, 1, 5, 4, 6)
  tmp <-
    as_tibble(plotresult)[plot_id, S3index] %>% as.list() %>% unlist(recursive = F)
  
  ggsave(
    "./figS2e_h.pdf",
    wrap_plots(tmp, nrow = 5, byrow = F) + plot_layout(heights = c(1, 0.3, 1, 1, 3)),
    dpi = 600,
    width = 18,
    height = 23,
    limitsize = F
  )
  
}
