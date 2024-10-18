library(Seurat)
library(magrittr)
library(ggplot2)
library(harmony)
setwd("~/cortex/figS1-6/")
seu <- qs::qread("../STEREO/Seu_merge_31_slides.qs")

seu %<>% NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData() %>% RunPCA()

seu %<>% RunUMAP(dims = 1:40)
seu %<>% RunHarmony(ncores = 20,group.by.vars = "chipID", dims.use = 1:40)
seu %<>% RunUMAP(dims = 1:40,reduction = "harmony")
seu %>% qs::qsave("./Seu_merge_31_slides.qs")

seu_grey_matter <- seu[,seu$mainLayer %in% c("ARACHNOID","WM")]
seu_grey_matter %<>% ScaleData(features = rownames(seu_grey_matter))

ACC.markers <- FindMarkers(seu_grey_matter,group.by = "dissection",ident.1 = "ACC", ident.2 = NULL)
geneName_ID <- read.csv("../../ensemble93gtf_rmXY.csv") %>% {setNames(.$gene_id ,nm = .$gene_name)} 
geneID_name <- read.csv("../../ensemble93gtf_rmXY.csv") %>% {setNames(nm = .$gene_id ,.$gene_name)} 

ACC.markers$gene <- rownames(ACC.markers) 
ACC.markers$geneID <- geneName_ID[ACC.markers$gene]
ACC.markers %<>% drop_na()
ACC.markers %<>% filter(avg_log2FC >0 )
qs::qsave(ACC.markers,"ACC.markers.qs")

ast_markers <- qs::qread("~/DATA/data/SnRNA/ASTRO_marker_edlein.qs")

sharedGenes <- intersect(ast_markers$gene,ACC.markers$geneID)
ast_markers %<>% filter(gene %in% sharedGenes) 
ast_markers$geneName <- geneID_name[ast_markers$gene]

ast_markers %<>% filter(cluster %in% c("Astro_4","Astro_5"),avg_log2FC > 0.25,p_val_adj <0.05) %>% group_by(gene) %>% slice_max(order_by = avg_log2FC,n = 1) %>% arrange(cluster,desc(avg_log2FC)) %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC,n = 30)  

exp <- AverageExpression(seu_grey_matter,slot = "scale.data",group.by = "dissection",features = ast_markers$geneName) %$%Spatial
exp <- as.matrix(exp)
exp %<>% scales::oob_squish(range = c(-3,3))

colSplit = ifelse(colnames(exp) == "ACC",yes = "ACC", no = "others")
pdf("FigS4m.pdf",height = 15,width = 9)
ComplexHeatmap::Heatmap(exp, cluster_rows = F,row_gap = unit(5, "mm"),column_gap = unit(5, "mm"),column_split = colSplit,border = T,row_split = ast_markers$cluster,cluster_columns = F,col = c("black","purple","yellow"))
dev.off()
