library(qs)
library(parallel)
library(magrittr)
library(tidyverse)
library(org.Hs.eg.db)
setwd("~/cortex/fig6/")
devtools::load_all("~/seurat/")
devtools::load_all("~/ClusterGVis-main/")
library(scRNAtoolVis)

seu <- qread("../STEREO/st_domain_seu_44slides.qs")

genes <- read.csv("../../ensemble93gtf_rmXY.csv")
 
region_color <- c(FPPFC = "#3F4587", DLPFC = "#8562AA", VLPFC = "#EC8561", M1 = "#B97CB5", 
                  S1 = "#D43046", S1E = "#F0592B", PoCG = "#ED4A96", SPL = "#593C97", 
                  SMG = "#A54486", AG = "#FBDE13", V1 = "#299FAE", ITG = "#75CCE3", 
                  STG = "#0C6939", ACC = "#0D9547")
 
seu <- seu[,seu$majorDomain == "L4"]
 
Idents(seu) <- seu$region %>% factor(levels = c("FPPFC", "DLPFC", "VLPFC", "M1", "S1", "S1E", "PoCG", "SPL", "SMG", "AG", "V1", "ITG", "STG", "ACC"))

Markers_different_region_L4 <- mclapply(Idents(seu)  %>% as.character %>% unique,mc.cores = min(length(levels(Idents(seu))),20), FUN = function(domain){
  tmp2 <- FindMarkers(seu,assay = "Spatial",ident.1 = domain,ident.2 = NULL,min.pct = 0.1,logfc.threshold = 0.1, only.pos = T, group.by = "region")
  tmp2$domain <- domain
  tmp2$gene <- rownames(tmp2)
  tmp2
})
# 
# Markers_different_region_L4_conserved <- mclapply(Idents(seu)  %>% as.character %>% unique,mc.cores = min(length(levels(Idents(seu))),20), FUN = function(domain){
#   tmp2 <- FindConservedMarkers(seu,assay = "Spatial",ident.1 = domain,ident.2 = NULL, grouping.var = "donor")
#   tmp2$domain <- domain
#   tmp2$gene <- rownames(tmp2)
#   tmp2
# })


Markers_different_region_L4 %>% purrr::reduce(rbind)  %>% write_csv("fig6f_Layer4_regional_marker.csv")

Markers_different_region_L4 <- read_csv("fig6f_Layer4_regional_marker.csv") %>% filter(abs(avg_log2FC) > 1, p_val_adj < 0.05) %>%  dplyr::group_by(gene) %>%
  dplyr::top_n(n = 1,wt = avg_log2FC)


Markers_different_region_L4$cluster <- Markers_different_region_L4$domain
Markers_different_region_L4$cluster %<>% as.character %>% factor(levels =  c("FPPFC", "DLPFC", "VLPFC", "M1", "S1", "S1E", "PoCG", "SPL", "SMG", "AG", "V1", "ITG", "STG", "ACC"))
Markers_different_region_L4 %<>% arrange(cluster)
# Markers_different_region_L4$avg_log2FC %<>% scales::oob_squish(range = c(-5,5))
Markers_different_region_L4 %<>% filter(str_detect(gene,"\\.[0-9]$",negate = T))

{jjVolcano(diffData = Markers_different_region_L4,
           log2FC.cutoff = 1, size  = 3.5,
           pSize = 0.2,base_size = 5,legend.position = c(0.1,0.9), tile.col = region_color,
           topGeneN = 2,polar = T,cluster.order = c( "DLPFC","VLPFC", "M1", "S1", "S1E", "PoCG", "SPL", "SMG", "AG", "V1", "ITG", "STG", "ACC","FPPFC") %>% rev)} %>% ggsave(filename = "fig6f_L4_gene_downup_circle.pdf",height = 10,width = 10)

# {jjVolcano(diffData = Markers_different_region_L4,
           # log2FC.cutoff = 0.5,
           # pSize = 0.2,base_size = 2,legend.position = c(0.1,0.9), tile.col = region_color,
           # topGeneN = 2,polar = T)} %>% ggsave(filename = "fig6f_L4_gene_downup_circle.pdf",height = 10,width = 10)

# fig 6g


seu_L4 <- seu[,seu$region == "V1"]

Idents(seu_L4) <- seu_L4$domain %>% factor()

Markers_different_region_L4abc <- mclapply(Idents(seu_L4)  %>% as.character %>% unique,mc.cores = min(length(levels(Idents(seu_L4))),20), FUN = function(domain){
  tmp2 <- FindMarkers(seu_L4,assay = "Spatial",ident.1 = domain,ident.2 = NULL,min.pct = 0.1)
  tmp2$domain <- domain
  tmp2$gene <- rownames(tmp2)
  tmp2
  
}) %>% purrr::reduce(rbind)

write.csv(Markers_different_region_L4abc,"./fig6g_L4abc_DEGs.csv")

Markers_different_region_L4abc <- read_csv("fig6g_L4abc_DEGs.csv") %>% filter(p_val_adj<0.05) %>% filter(str_detect(gene,"\\.[0-9]$",negate = T))

library(scRNAtoolVis)

Markers_different_region_L4abc$cluster <- Markers_different_region_L4abc$domain
Markers_different_region_L4abc %<>% arrange(domain)
Markers_different_region_L4abc$cluster %<>% factor()
Markers_different_region_L4abc$avg_log2FC %<>% scales::oob_squish(range = c(-5,5))

{jjVolcano(diffData = Markers_different_region_L4abc,
           log2FC.cutoff = 0.2,
           pSize = 0.3,base_size = 10,legend.position = c(0.1,0.9), tile.col = c("L4a" = "#EC8561", "L4b" = "#B97CB5", "L4c" = "#3F4587"),
           topGeneN = 2,polar = T)} %>% ggsave(filename = "fig6f_L4abc_gene_downup_circle.pdf",height = 7,width = 7)

# jjVolcano(diffData = domainMarkers.all, #%>% filter(str_detect(gene,"[0-9]{1,}$",negate = T)),
#           log2FC.cutoff = 0.1,pSize = 0.5,base_size = 10,legend.position = c(0.1,0.9), tile.col = colorPallete,
#           topGeneN = 3,polar = T) %>% ggsave(filename = "L4abc_gene_downup_circle.pdf",f5_3,height = 7,width = 7)

# cluster_DEGs <- Markers_different_region_L4abc %>% group_by(cluster) %>% summarise(up = sum(avg_log2FC > 0), down = -sum(avg_log2FC < 0))

# f6h
domainMarkers.up <- Markers_different_region_L4abc %>% filter(avg_log2FC> 0)



library(ClusterGVis)
st.data <- prepareDataFromscRNA(object = seu_L4,
                                assays = "Spatial",
                                diffData = domainMarkers.up,
                                showAverage = TRUE)



library(org.Hs.eg.db)
enrich <- enrichCluster(object = st.data,
                        OrgDb = org.Hs.eg.db,
                        type = "BP",
                        organism = "hsa",
                        pvalueCutoff = 0.5,
                        topn = 10,
                        seed = 5201314)


enrich %>% write_csv("fig6h_L4abc_GO.csv")

markGenes <- domainMarkers.up %>% 
  filter(str_detect(gene,negate = T,"^A[A-Z].*\\.[0-9]") ) %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 3, wt = avg_log2FC) %$% gene 

markGenes <- c(markGenes,c("EGR1", "ENC1", "HOPX", "NEFH", "SYT2", "SMYD2", "CD74", "RAB3B", "CNTN5")) %>% unique 


pdf('fig6gL4abc_enrichGO.pdf',height = 9,width = 11,onefile = F)
visCluster(object = st.data,
           plot.type = "both",
           column_names_rot = 45,
           show_row_dend = F,
           markGenes = markGenes,
           markGenes.side = "left",
           annoTerm.data = enrich,
           line.side = "left",
           add.box = T,
           ctAnno.col = c("L4a" = "#EC8561", "L4b" = "#B97CB5", "L4c" = "#3F4587"),
           sample.col = c("L4a" = "#EC8561", "L4b" = "#B97CB5", "L4c" = "#3F4587"),
           cluster.order = c(1:14),
           go.col = c("L4a" = "#EC8561", "L4b" = "#B97CB5", "L4c" = "#3F4587")[enrich$group %>% str_extract("[0-9]{1,}") %>% as.numeric()],
           add.bar = T)
dev.off()

FetchData()

tmp <- domainMarkers.up %>% filter(cluster == "L4b") 

seu_L4 %<>% ScaleData(features = unique(domainMarkers.up$gene),split.by = "chip")

cl = "L4a"

mg <- c("EGR1","ENC1","HOPX")

tmp <- SpatialFeaturePlot(seu_L4,images = c("V1_B01012C5"),stroke = 0,slot = "scale.data",pt.size.factor = 3,features = mg,alpha = c(0,1),combine = F)  %>% purrr::map(~. + coord_flip() + scale_y_reverse() )
tmp <- patchwork::wrap_plots(tmp, nrow = 1,byrow = T) 
ggsave(str_c("fig6h_L4abc_",cl,"_V1_B01012C5_markers.png"),tmp, height = 21/8,width =  92/8/2,dpi = 150,limitsize = FALSE )

tmp <- SpatialFeaturePlot(seu_L4,images = c("V1_B02222E1"),stroke = 0,slot = "scale.data",pt.size.factor = 3,features = mg,alpha = c(0,1),combine = F)  %>% purrr::map(~. + coord_flip() + scale_y_reverse() )
tmp <- patchwork::wrap_plots(tmp, nrow = 1,byrow = T) 
ggsave(str_c("fig6h_L4abc_",cl,"_V1_B02222E1_markers.png"),tmp, height = 21/8,width =  92/8/2,dpi = 150,limitsize = FALSE )

cl = "L4b"
mg <- c("NEFH","SYT2","NEFM")

tmp <-{ SpatialFeaturePlot(seu_L4,images = c("V1_B01012C5"),stroke = 0,slot = "scale.data",pt.size.factor = 3,features = mg,alpha = c(0,1),combine = F) } %>% purrr::map(~. + scale_fill_gradientn(limits = c(-3,3),oob = scales::oob_squish,colours = RColorBrewer::brewer.pal(name = "Spectral",n = 11) %>% rev) + coord_flip() + scale_y_reverse() )
tmp <- patchwork::wrap_plots(tmp, nrow = 1,byrow = T) 

ggsave(str_c("fig6h_L4abc_",cl,"_V1_B01012C5_markers.png"),tmp, height = 21/8,width =  92/8/2,dpi = 150,limitsize = FALSE )

tmp <- SpatialFeaturePlot(seu_L4,images = c("V1_B02222E1"),stroke = 0,slot = "scale.data",pt.size.factor = 3,features = mg,alpha = c(0,1),combine = F)  %>% purrr::map(~.  + scale_fill_gradientn(limits = c(-3,3),oob = scales::oob_squish,colours = RColorBrewer::brewer.pal(name = "Spectral",n = 11) %>% rev) + coord_flip() + scale_y_reverse() )
tmp <- patchwork::wrap_plots(tmp, nrow = 1,byrow = T) 
ggsave(str_c("fig6h_L4abc_",cl,"_V1_B02222E1_markers.png"),tmp, height = 21/8,width =  92/8/2,dpi = 150,limitsize = FALSE )

cl = "L4c"

mg <- c("EGR1","ENC1","HOPX","NEFH","SYT2","CD74","RAB3B","CNTN5")


tmp <- SpatialFeaturePlot(seu_L4,images = c("V1_B01012C5"),stroke = 0,slot = "scale.data",pt.size.factor = 3,features = mg,alpha = c(0,1),combine = F)  %>% purrr::map(~. + coord_flip() + scale_y_reverse() )
tmp <- patchwork::wrap_plots(tmp, nrow = 1,byrow = T) 
ggsave(str_c("fig6h_L4abc_",cl,"_V1_B01012C5_markers.png"),tmp, height = 21/8,width =  92/8/2,dpi = 150,limitsize = FALSE )

tmp <- SpatialFeaturePlot(seu_L4,images = c("V1_B02222E1"),stroke = 0,slot = "scale.data",pt.size.factor = 3,features = mg,alpha = c(0,1),combine = F)  %>% purrr::map(~. + coord_flip() + scale_y_reverse() )
tmp <- patchwork::wrap_plots(tmp, nrow = 1,byrow = T) 
ggsave(str_c("fig6h_L4abc_",cl,"_V1_B02222E1_markers.png"),tmp, height = 21/8,width =  92/8/2,dpi = 150,limitsize = FALSE )





{ ggplot(plotData,aes(x,y,color = scale.exp, alpha = scale.exp)) + geom_point(size = 0.1) + facet_wrap(~chip+gene,ncol = 3) + scale_color_gradientn(colors = RColorBrewer::brewer.pal(n = 11, "Spectral") %>% rev,limits = c(-3,3)) + cowplot::theme_map()} %>% ggsave(filename = "fig6i.pdf",height = 3.2*27,width = 8,limitsize = F)
