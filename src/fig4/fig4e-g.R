library(qs)
library(tidyverse)
library(Seurat)
library(magrittr)
setwd("~/cortex/fig4/")

sst <- qread("sstRNA.qs")
Idents(sst) <- "depth_cluster"
geneID_name <- read_csv("../SnRNA/1_SnRNA_preprocessing/gene_kept.csv") %>% {setNames(object = .$gene_uni,nm = .$gene_id)}
clusterMarkers <- FindAllMarkers(sst,only.pos = T)
clusterMarkers %>% write_csv("sst_degs.csv")
clusterMarkers <- read_csv("sst_degs.csv")

clusterMarkers$sym <- geneID_name[clusterMarkers$gene]

sst$cl = str_c(sst$depth_cluster,sst$cluster)

clusterMarkers.selected <- clusterMarkers  %>% group_by(cluster) %>% slice_max(n = 20,order_by = avg_log2FC)  %>% group_by(gene) %>% slice_max(order_by = avg_log2FC) 
clusterMarkers.selected <- rbind(clusterMarkers.selected,clusterMarkers %>% filter(sym %in% c("PDZD2", "GNAL", "GRIA4", "CALB1", "DCC", "TRHDE", "GRIN3A"))) %>% distinct()
sst %<>% ScaleData(features = clusterMarkers.selected$gene)
exp <- AverageExpression(sst,slot = "data",group.by = "cl",features = clusterMarkers.selected$gene) %$% RNA

exp <- as.matrix(exp)

exp <- exp[,colnames(exp) %>% sort]
rownames(exp) <- geneID_name[rownames(exp)]
colSplit <- str_extract(colnames(exp),"shallow|deep")
exp <- t(scale(t(exp)))
exp %<>% scales::oob_squish(range = c(-3,3))
colnames(exp) %<>% str_remove("deep|shallow")

rowLabels <- rownames(exp)
# rowLabels[!rownames(exp) %in% c("GNAL","GRIA4","DCC","TRHDE","GRIN3A","CALB1","PDZD2")] = ""


pdf("fig4e.pdf",height = 13,width = 10)
ComplexHeatmap::Heatmap(exp, cluster_rows = T,row_labels = rowLabels,row_gap = unit(5, "mm"),column_gap = unit(5, "mm"),
                        column_split = colSplit,border = T,row_split = clusterMarkers.selected$cluster,
                        cluster_columns = T,col = c("black","purple","yellow"))
dev.off()

# fig4f
{FeaturePlot(sst,features = "ENSG00000104327") + scale_color_gradientn(colours = MetBrewer::MetPalettes$Greek[[1]] %>% rev) + labs(title = "CALB1") + coord_fixed() + cowplot::theme_map() } %>% # CALB1
  {ggsave(filename = "fig4f.pdf",width = 5,height = 5)}

# fig4g

library(enrichR)

deep <- enrichR::enrichr(genes = clusterMarkers %>% filter(cluster == "deep", p_val_adj < 0.05) %$% sym, databases = "GO_Biological_Process_2021" )

shallow <- enrichR::enrichr(genes = clusterMarkers %>% filter(cluster == "shallow", p_val_adj < 0.05) %$% sym, databases = "GO_Biological_Process_2021" )

deep$GO_Biological_Process_2021$type = "deep"
shallow$GO_Biological_Process_2021$type = "shallow"

enrichtable = rbind(deep$GO_Biological_Process_2021,shallow$GO_Biological_Process_2021)
enrichtable %<>% group_by(Term) %>% slice_min(Adjusted.P.value)
enrichtable %<>% group_by(type) %>% slice_min(Adjusted.P.value,n = 6)
enrichtable$Term %<>% str_remove("\\(.*\\)")
enrichtable$Combined.Score = enrichtable$Combined.Score * ifelse(enrichtable$type == "deep", -1,1)
enrichtable %<>% arrange(type,Combined.Score)
enrichtable$Term  %<>%  str_wrap(width = 40) 

enrichtable$Term %<>% factor(levels = enrichtable$Term)

{
  ggplot(enrichtable, aes(x = Term, y = Combined.Score, fill = type)) + 
    geom_bar(stat = "identity") +  cowplot::theme_cowplot() +  MetBrewer::scale_fill_met_d("VanGogh2") + ggh4x::coord_axes_inside()  + theme(axis.text.x = element_text(
    angle = 90,
    hjust = 1,
    vjust = 0.5,
    lineheight = 0.6
  ))
} %>%
  {ggsave(filename = "fig4g.pdf" ,height = 8,width = 6)}

