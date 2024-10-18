library(tidyverse)
set("~/cortex/figS1-6/")
devtools::load_all("~/scrattch.bigcat/")
library(scrattch.hicat)

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
names(region_color) <- region

merge_seu <- readRDS("../SnRNA/SnRNA_seurat.RDS")
seu <-  merge_seu[,merge_seu$class == "excitatory"]
cl.region <- seu@meta.data$region %>% setNames(seu %>% colnames())

source("./knn.R")

knn_graph.rg <-
  get_knn_graph(seu@reductions$pca@cell.embeddings, cl = cl.region[colnames(seu)])

knn.cl.df.rg <- knn_graph.rg$knn.cl.df %>% arrange(desc(frac))

threshold <- knn.cl.df.rg %>% filter(cl.from != cl.to) %$% frac %>% quantile(probs = 0.92)

region.df <-
  data.frame(
    cluster_label = region,
    cluster_color = region_color[region],
    cluster_id = region,
    cl = region,
    row.names = region
  )

regionLoc <-
  data.frame(
    x = c(87, 142, 207, 225, 340, 353, 384, 350, 410, 425, 460, 500, 480, 626),
    y = c(50, 301, 159, 251, 50, 187, 100, 233, 300, 50, 112, 75, 172, 270),
    row.names = c(
      "ACC",
      "FPPFC",
      "DLPFC",
      "VLPFC",
      "M1",
      "S1E",
      "S1",
      "STG",
      "ITG",
      "PoCG",
      "SMG",
      "SPL",
      "AG",
      "V1"
    )
  )

purrr::walk(1:nrow(region.df), function(i) {
  region.df$x[i] <<- regionLoc[region.df$cluster_label[i], "x"] / 10
  region.df$y[i] <<- regionLoc[region.df$cluster_label[i], "y"] / 10
  region.df$cluster_size[i] <<-
    sum(cl.region == region.df$cluster_label[i])
})

backgroud <- png::readPNG("./brain.png")
backgroud[,,4] <- 0.8 # x = 100, 650, y = 20, 图像是反过来的
xlim = c(0, dim(backgroud)[2])/10
ylim = c(0, dim(backgroud)[1])/10

plot_t <-
  plot_constellation_2(knn.cl.df = knn.cl.df.rg %>% filter(frac > threshold) ,
                       cl.center.df = region.df,
                       exxageration = 1,
                       curved = T,
                       label.size = 8,
                       max_size = 40,
                       node.dodge = T) +
  scale_y_reverse() +
  coord_fixed(expand = F, xlim = xlim, ylim = rev(ylim)) +
  annotation_custom(grob = grid::rasterGrob(backgroud)) + cowplot::theme_map()

plot_t$layers <- c(plot_t$layers[[4]],plot_t$layers[1:3])
plot_t

ggsave("figS2d.pdf",plot_t,height = 15,width = 15,dpi = 50)

