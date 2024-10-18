library(qs)
library(parallel)
library(magrittr)
library(tidyverse)
library(org.Hs.eg.db)

setwd("~/cortex/figS1-6/")

devtools::load_all("~/ClusterGVis-main/")
devtools::load_all("~/seurat/")

seu <- qread("../STEREO/st_domain_seu_44slides.qs")

genes <- read.csv("../STEREO/ensemble93gtf_rmXY.csv")

colorPallete <- c(ggsci::pal_aaas()(10),
                  ggsci::pal_jama()(10),
                  ggsci::pal_npg()(10),
                  ggsci::pal_lancet()(10),
                  ggsci::pal_frontiers()(10),
                  ggsci::pal_nejm()(10)) %>% unique %>% str_subset(negate = T, "1B1919")

seu_merge <- seu[,seu$domain == "L6b"]

Markers_different_region_L6b <- parallel::mclapply(seu_merge$region %>% unique,mc.cores = 14,FUN = function(x){
  tmp <- FindMarkers(seu_merge,ident.1 = x,ident.2 = NULL,group.by = "region")
  tmp$cluster <- x
  tmp$gene <- rownames(tmp)
  tmp
}) %>% purrr::reduce(rbind)


Markers_different_region_L6b %>% write_csv("./figS6_Layer6b_regional_marker.csv")

Markers_different_region_L6b <- read_csv("./figS6_Layer6b_regional_marker.csv") %>% filter(abs(avg_log2FC) > 1, p_val_adj < 0.05) %>%  dplyr::group_by(gene) %>%
  dplyr::top_n(n = 1,wt = avg_log2FC)


Markers_different_region_L6b$cluster %<>% as.character %>% factor(levels =  c("FPPFC", "DLPFC", "VLPFC", "M1", "S1", "S1E", "PoCG", "SPL", "SMG", "AG", "V1", "ITG", "STG", "ACC"))
Markers_different_region_L6b %<>% arrange(cluster)
Markers_different_region_L6b$avg_log2FC %<>% scales::oob_squish(range = c(-3,3))
Markers_different_region_L6b %<>% filter(str_detect(gene,"\\.[0-9]$",negate = T))

{jjVolcano(diffData = Markers_different_region_L6b,
           log2FC.cutoff = 0.5,
           pSize = 0.2,base_size = 2,legend.position = c(0.1,0.9), tile.col = region_color,
           topGeneN = 2)} %>% ggsave(filename = "figS6a_L6_gene_region.png",height = 3,width = 10)

{jjVolcano(diffData = Markers_different_region_L6b,
           log2FC.cutoff = 0.5,
           pSize = 0.2,base_size = 2,legend.position = c(0.1,0.9), tile.col = region_color,
           topGeneN = 2)} %>% ggsave(filename = "figS6a_L6_gene_region.pdf",height = 3,width = 10)

Idents(seu_merge) <- seu_merge$region %>% factor(levels = c("FPPFC", "DLPFC", "VLPFC", "M1", "S1", "S1E", "PoCG", "SPL", "SMG", "AG", "V1", "ITG", "STG", "ACC"))

st.data <- prepareDataFromscRNA(object = seu_merge,
                                assays = "Spatial",
                                diffData = Markers_different_region_L6b,
                                showAverage = TRUE)

enrich <- enrichCluster(object = st.data,
                        OrgDb = org.Hs.eg.db,
                        type = "BP",
                        organism = "hsa",
                        pvalueCutoff = 0.5,
                        topn = 10,
                        seed = 5201314)


write_csv(enrich,"figS6L6b_GO.csv")

enrichPlot <- read.csv("figS6L6b_GO.csv") %>% filter(str_detect(Description,"animal|toxin|viral|polysaccharide|tumor|eye|development|leukocyte|cardiac|muscle|leukocyte|fever|respiratory|heart|respiration|lymphocyte|vitamin",negate = T)) %>% group_by(group) %>% top_n(n = 3, pvalue) %>% dplyr::select(group,Description)

markGenes <- Markers_different_region_L6b %>% 
  filter(str_detect(gene,negate = T,"^A[A-Z].*\\.[0-9]") ) %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 3, wt = avg_log2FC) %$% gene 


region_color <- c("#2F4587", "#8562AA", "#EC8561", "#B97CB5", "#D43046", "#F0592B", 
                  "#ED4A96", "#593C97", "#A54486", "#FBDE13", "#299FAE", "#75CCE3", 
                  "#0C6939", "#0D9547") %>% setNames(c("FPPFC", "DLPFC", "VLPFC", "M1", "S1", "S1E", "PoCG", "SPL",
                                                       "SMG", "AG", "V1", "ITG", "STG", "ACC"))
pdf('figS6b.pdf',height = 10.5,width = 12,onefile = F)

visCluster(object = st.data,
           plot.type = "both",
           column_names_rot = 45,
           show_row_dend = F,
           markGenes = markGenes,
           sample.col = region_color,
           ctAnno.col = region_color[c("FPPFC", "DLPFC", "VLPFC", "M1", "S1", "S1E", "PoCG", "SPL",
                                       "SMG", "AG", "V1", "ITG", "STG", "ACC")],
           markGenes.side = "left",
           annoTerm.data = enrichPlot,
           line.side = "left",
           cluster.order = st.data$wide.res$cluster %>% unique %>% as.numeric() %>% sort,
           go.col = region_color[c("FPPFC", "DLPFC", "VLPFC", "M1", "S1", "S1E", "PoCG", "SPL",
                                   "SMG", "AG", "V1", "ITG", "STG", "ACC")[enrichPlot$group %>% str_extract("[0-9]{1,}") %>% as.numeric()]])

dev.off()

