library(Seurat)
library(tidyverse)
library(ranger)
set("~/cortex/figS1-6/")

seu <- readRDS('../SnRNA/SnRNA_seurat.RDS')
seu <- seu[,seu$class == "excitatory"]


train_data <-
  seu@reductions$pca@cell.embeddings %>% as.data.frame()

train_data["label"] <- seu$region %>% factor

train.forest <- ranger(label ~ .,
                             data = train_data)


tmp_plot <-
  train.forest$confusion %>% as.data.frame() %>%  group_by(true) %>% 
  mutate(predicted_rate = Freq/sum(Freq))

tmp_plot %<>% pivot_wider(id_cols = !Freq,names_from = predicted,values_from = predicted_rate) 
tmp_plot %<>% column_to_rownames("true")
tmp_plot <- scale(t(tmp_plot)) %>% t


p_dend <- try({pvclust::pvclust(tmp_plot %>% t())} )

p_dend <- p_dend$hclust %>% ape::as.phylo() %>% ggtree()
region_order <- get_taxa_name(p_dend)

pdf("figS2c_1.pdf",width = 10,height = 3)
p_dend + layout_dendrogram() + geom_tiplab() + ylim(c(20,-5)) + xlim(c(1,-1))
dev.off()

tmp_plot_long <-
  tmp_plot %>% as.data.frame %>% rownames_to_column("original") %>% pivot_longer(cols = !original,
                                                                                names_to = "predicted",
                                                                                values_to = "percentage")
tmp_plot_long$original %<>% factor(levels = region_order)
tmp_plot_long$predicted %<>% factor(levels = region_order)


pdf("figS2c_2.pdf",width = 5,height = 4)
ComplexHeatmap::Heatmap(tmp_plot,row_order = region_order,column_order = region_order,col = c("black","purple","white","red"))
dev.off()

