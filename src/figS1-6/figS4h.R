library(org.Hs.eg.db)
library(Seurat)
library(magrittr)
library(tidyverse)
library(clusterProfiler)
library(enrichR)
setwd("~/cortex/figS1-6/")
devtools::load_all("~/ClusterGVis-main/")

# add cell type
geneid_name <- read.csv("../SnRNA/1_SnRNA_preprocessing/gene_kept.csv") %>% {setNames(.$gene_name,.$gene_id)}

region_color <- c(FPPFC = "#3F4587", DLPFC = "#8562AA", VLPFC = "#EC8561", M1 = "#B97CB5", 
                  S1 = "#D43046", S1E = "#F0592B", PoCG = "#ED4A96", SPL = "#593C97", 
                  SMG = "#A54486", AG = "#FBDE13", V1 = "#299FAE", ITG = "#75CCE3", 
                  STG = "#0C6939", ACC = "#0D9547")


seu <- qs::qread("../SnRNA/3_mergingDatasets/2_OPC_merge.qs")

Idents(seu) <- "region"
Idents(seu) %<>% fct_relevel(c(
  'FPPFC',  'DLPFC',  'VLPFC',  'M1',  'S1',  'S1E',  'PoCG',  'SPL',  'SMG',  'AG',  'V1',  'ITG',  'STG',  'ACC'))


markers <- parallel::mclapply(seu$region %>% unique,mc.cores = 14,FUN = function(x){
  tmp <- FindMarkers(seu,ident.1 = x,ident.2 = NULL,group.by = "region",only.pos = T)
  tmp$cluster <- x
  tmp$gene <- rownames(tmp)
  tmp
}) %>% purrr::reduce(rbind)

table(markers$gene %in% names(geneid_name))
# TRUE 
# 55352 

markers$sym <- geneid_name[markers$gene]

markers %>% write_csv("figS4h_OPC_region_DEGs.csv")
markers <- read_csv("figS4h_OPC_region_DEGs.csv")

markers$cluster %<>% fct_relevel(c(
  'FPPFC',  'DLPFC',  'VLPFC',  'M1',  'S1',  'S1E',  'PoCG',  'SPL',  'SMG',  'AG',  'V1',  'ITG',  'STG',  'ACC'))

markers %<>% filter(p_val_adj < 0.05)%>%
  dplyr::group_by(gene) %>%
  dplyr::top_n(n = 1,wt = avg_log2FC)

markers %<>% arrange(cluster)

st.dataPlot <- prepareDataFromscRNA(object = seu,
                                    diffData = markers,
                                    showAverage = T)

st.dataPlot$wide.res$gene <- geneid_name[st.dataPlot$wide.res$gene]

enrich <- ClusterGVis::enrichCluster(st.dataPlot,type = "BP",OrgDb = org.Hs.eg.db)

table(duplicated(st.dataPlot$wide.res$gene))

write.csv(enrich,"figS4h_OPC_GO_Term.csv")

enrich <- read.csv("figS4h_OPC_GO_Term.csv") %>%  filter(str_detect(Description,"skin|bone",negate = T)) %>% group_by(group) %>% arrange(pvalue) %>% slice_head(n = 4) %>% select(group,Description)

markGenes <- markers %>% filter(avg_log2FC >0) %>% mutate(geneSym = geneid_name[gene]) %>% filter(str_detect(geneSym,"\\.[0-9]$",negate = T)) %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 2, wt = avg_log2FC) %>% {geneid_name[.$gene]}

table(duplicated(markGenes))


pdf('figS4h.pdf',height = 10.5,width = 11,onefile = F)

visCluster(object = st.dataPlot,
           plot.type = "both",
           add.line = F,
           column_names_rot = 45,
           markGenes = markGenes,
           ctAnno.col = region_color[c('FPPFC',  'DLPFC',  'VLPFC',  'M1',  'S1',  'S1E',  'PoCG',  'SPL',  'SMG',  'AG',  'V1',  'ITG',  'STG',  'ACC')],
           sample.col = region_color[c('FPPFC',  'DLPFC',  'VLPFC',  'M1',  'S1',  'S1E',  'PoCG',  'SPL',  'SMG',  'AG',  'V1',  'ITG',  'STG',  'ACC')],
           markGenes.side = "left",
           annoTerm.data = enrich,
           line.side = "left",
           cluster.order = 1:14,
           go.col = region_color[c('FPPFC',  'DLPFC',  'VLPFC',  'M1',  'S1',  'S1E',  'PoCG',  'SPL',  'SMG',  'AG',  'V1',  'ITG',  'STG',  'ACC')][enrich$group %>% str_extract("[0-9]{1,}") %>% as.numeric()],
 )

dev.off()
