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
 
seu_merge <- seu[,seu$majorDomain == "L4"]


Markers_different_region_L4 <- read_csv("../fig6/fig6f_Layer4_regional_marker.csv") %>% filter(abs(avg_log2FC) > 1, p_val_adj < 0.05) %>%  dplyr::group_by(gene) %>%
  dplyr::top_n(n = 1,wt = avg_log2FC)


Markers_different_region_L4$cluster <- Markers_different_region_L4$domain
Markers_different_region_L4$cluster %<>% as.character %>% factor(levels =  c("FPPFC", "DLPFC", "VLPFC", "M1", "S1", "S1E", "PoCG", "SPL", "SMG", "AG", "V1", "ITG", "STG", "ACC"))
Markers_different_region_L4 %<>% arrange(cluster)
Markers_different_region_L4$avg_log2FC %<>% scales::oob_squish(range = c(-5,5))
Markers_different_region_L4 %<>% filter(str_detect(gene,"\\.[0-9]$",negate = T))

Idents(seu_merge) <- seu_merge$region %>% factor(levels = c("FPPFC", "DLPFC", "VLPFC", "M1", "S1", "S1E", "PoCG", "SPL", "SMG", "AG", "V1", "ITG", "STG", "ACC"))

st.data <- prepareDataFromscRNA(object = seu_merge,
                                assays = "Spatial",
                                diffData = Markers_different_region_L4,
                                showAverage = TRUE)

enrich <- enrichCluster(object = st.data,
                        OrgDb = org.Hs.eg.db,
                        type = "BP",
                        organism = "hsa",
                        pvalueCutoff = 0.5,
                        topn = 10,
                        seed = 5201314)


# enrich %<>% purrr::reduce(rbind)
write_csv(enrich,"figS5L4_GO.csv")

enrichPlot <- read.csv("figS5L4_GO.csv") %>% filter(str_detect(Description,"viral|polysaccharide|tumor|eye|development|leukocyte|cardiac|leukocyte|fever|respiratory|heart|respiration|lymphocyte|vitamin",negate = T)) %>% group_by(group) %>% top_n(n = 5, pvalue) %>% dplyr::select(group,Description)

markGenes <- Markers_different_region_L4 %>% 
  filter(str_detect(gene,negate = T,"^A[A-Z].*\\.[0-9]") ) %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 3, wt = avg_log2FC) %$% gene 


region_color <- c("#2F4587", "#8562AA", "#EC8561", "#B97CB5", "#D43046", "#F0592B", 
                  "#ED4A96", "#593C97", "#A54486", "#FBDE13", "#299FAE", "#75CCE3", 
                  "#0C6939", "#0D9547") %>% setNames(c("FPPFC", "DLPFC", "VLPFC", "M1", "S1", "S1E", "PoCG", "SPL",
                                                       "SMG", "AG", "V1", "ITG", "STG", "ACC"))
pdf('figS5.pdf',height = 15,width = 12,onefile = F)

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

