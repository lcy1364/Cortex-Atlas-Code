library(Seurat)
library(tidyverse)
library(parallel)
devtools::load_all("~/ClusterGVis-main/")
setwd("~/cortex/fig3")
merge_seu <- readRDS("../SnRNA/SnRNA_seurat.RDS")

IT_seu <- merge_seu[,merge_seu$subclass %in% c("L2-L3 IT LINC00507", "L3-L4 IT RORB", "L4-L5 IT RORB", "L6 IT") ]

IT_markers <- mclapply(IT_seu$subclass %>% unique,mc.cores = 4,FUN = function(x){
  tmp <- FindMarkers(IT_seu,ident.1 = x,ident.2 = NULL,group.by = "subclass")
  tmp$subclass <- x
  tmp$gene <- rownames(tmp)
  tmp
})

IT_markers %<>% purrr::reduce(rbind)

IT_markers$subclass %<>% factor(levels = c("L2-L3 IT LINC00507", "L3-L4 IT RORB", "L4-L5 IT RORB", "L6 IT"))

Idents(IT_seu) <- IT_seu$subclass %>% factor(levels = c("L2-L3 IT LINC00507", "L3-L4 IT RORB", "L4-L5 IT RORB", "L6 IT"))

# get top 10 genes
markers <- IT_markers %>% filter(avg_log2FC >0, p_val_adj < 0.05)%>%
  dplyr::group_by(gene) %>%
  dplyr::top_n(n = 1,wt = avg_log2FC)
# dplyr::top_n(n = 15, wt = avg_log2FC)

markers$cluster <- markers$subclass

markers %<>% arrange(cluster)

st.dataPlot <- prepareDataFromscRNA(object = IT_seu,
                                    diffData = markers,
                                    showAverage = T)


geneKept <- read_csv("../SnRNA/1_SnRNA_preprocessing/gene_kept.csv") %>% {setNames(object = .$gene_uni,nm = .$gene_id)}

st.dataPlot$wide.res$gene <- geneKept[st.dataPlot$wide.res$gene]

# table(duplicated(st.dataPlot$wide.res$gene))
# FALSE
# 240

# markers$geneSym <- geneKept[markers$gene]
# table(duplicated(markers$geneSym))
# FALSE  TRUE
# 240    45
markGenes <- markers %>% filter(avg_log2FC >0) %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 10, wt = avg_log2FC) %>% {geneKept[.$gene]}

table(duplicated(markGenes))

 
pdf("fig3e.pdf",height = 10,width = 7,onefile = F)
visCluster(object = st.dataPlot,
           plot.type = "both",
           column_names_rot = 45,
           show_row_dend = F,
           markGenes = markGenes,
           sample.col = c("#07D8D8","#09B2B2","#69B199","#19E3BE"),
           ctAnno.col =  c("#07D8D8","#09B2B2","#69B199","#19E3BE"),
           markGenes.side = "left",
            line.side = "left",
           cluster.order = c(1:4)
 )
dev.off()

# fig3f

scData <- FetchData(IT_seu,vars = c(str_c("PC_",1:50),"subclass","region"))

mod11 <-  aov(PC_1 + + PC_2 + PC_3 + PC_4 + PC_5 + PC_6 + PC_7 + PC_8 + PC_9 + PC_10 + PC_11 + PC_12 + PC_13 + PC_14 + PC_15 + PC_16 + PC_17 + PC_18 + PC_19 + PC_20 + PC_21 + PC_22 + PC_23 + PC_24 + PC_25 + PC_26 + PC_27 + PC_28 + PC_29 + PC_30 + PC_31 + PC_32 + PC_33 + PC_34 + PC_35 + PC_36 + PC_37 + PC_38 + PC_39 + PC_40 + PC_41 + PC_42 + PC_43 + PC_44 + PC_45 + PC_46 + PC_47 + PC_48 + PC_49 + PC_50~ region * subclass ,
              data = scData) %>% summary
tmp <- mod11[[1]] %>% as.data.frame() %>% rownames_to_column("class") %>% mutate(class = str_trim(class))
ggsave(ggplot(tmp %>% filter(class != "Residuals"),mapping = aes(x = fct_reorder(class,`Sum Sq`, .desc = T), y = `Mean Sq`, fill = class)) + geom_col() + ggprism::theme_prism() + ggsci::scale_fill_aaas() + labs(x = str_c("P value < 2e-16")),filename = "fig3f.pdf")

                                                                                                                                                                                                            