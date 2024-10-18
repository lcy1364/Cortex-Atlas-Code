library(qs)
library(parallel)
library(magrittr)
library(tidyverse)
library(enrichR)
library(org.Hs.eg.db)
library(rrvgo) 
library(ggplot2)
library(ggrastr)
devtools::load_all("~/seurat/") # too many bugs in seurat v5
devtools::load_all("~/ClusterGVis-main/")

setwd("~/cortex/fig2/")

seu_merge <- qread("../STEREO/st_domain_seu_44slides.qs")

domainColor <- c(
  ARACHNOID = "#8a3b35",  
  L1 = "#8ba28e",         
  L2 = "#9ec87e",         
  L3 = "#669c68",         
  L4 = "#67b8bb",         
  L5 = "#5687ac",         
  L6 = "#5d5d8d",         
  WM = "#b5b3bb")

# colorPallete <- c(
#   ggsci::pal_aaas()(10),
#   ggsci::pal_jama()(10),
#   ggsci::pal_npg()(10),
#   ggsci::pal_lancet()(10),
#   ggsci::pal_frontiers()(10),
#   ggsci::pal_nejm()(10)
# ) %>% unique %>% str_subset(negate = T, "1B1919")
# 

Idents(seu_merge) <-
  seu_merge$region %>% factor(
    levels = c(
      "FPPFC",
      "DLPFC",
      "VLPFC",
      "M1",
      "S1",
      "S1E",
      "PoCG",
      "SPL",
      "SMG",
      "AG",
      "V1",
      "ITG",
      "STG",
      "ACC"
    )
  )


# conserved markers
Idents(seu_merge) <- "majorDomain"
domain <- levels(seu_merge)[[2]]

majorMarkers <- mclapply(levels(seu_merge),mc.cores = min(length(levels(Idents(seu_merge))),20), FUN = function(domain){
  tmp <- FindMarkers(seu_merge, ident.1 = domain,assay = "Spatial", ident.2 = NULL,min.cells.group = 5,only.pos = T)
  tmp$domain <- domain
  tmp$gene <- rownames(tmp)
  tmp
}) %>% purrr::reduce(rbind)

majorMarkers %>% qsave("fig2c_mainLayerConservedMarkers.qs") 

majorMarkers <- qread("fig2c_mainLayerConservedMarkers.qs")

majorMarkers %<>% filter(avg_log2FC > 0, p_val_adj < 0.05) %>%  dplyr::group_by(gene) %>%
  dplyr::top_n(n = 1, wt = avg_log2FC)


majorMarkers$cluster <- majorMarkers$domain
majorMarkers$cluster %<>% factor(
  levels =  c("ARACHNOID", "L1", "L2", "L3", "L4", "L5", "L6", "WM")
)
majorMarkers %<>% arrange(cluster)

Idents(seu_merge)  %<>%  as.character %>% factor(levels =  c("ARACHNOID", "L1", "L2", "L3", "L4", "L5", "L6", "WM"))

st.data <- prepareDataFromscRNA(
  object = seu_merge,
  assays = "Spatial",
  diffData = majorMarkers,
  showAverage = TRUE
)

rgWithDegs <- majorMarkers$cluster %>% unique

i = 1

for (i in 1:length(rgWithDegs)) {
  if (i == 1) enrich <- list()
  rg = rgWithDegs[i]
  print(rg)
  timestamp()
  mgs <-
    majorMarkers %>% filter(cluster == rg, avg_log2FC > 0, p_val_adj < 0.05) %$% gene
  print(length(mgs))
  
  if (length(mgs) == 0) {
    next
  }
  
  tmp <- enrichR::enrichr(mgs,"GO_Biological_Process_2023")
  tmp <- tmp[[1]]
  tmp$Description <- tmp$Term
  tmp$group <- str_c("C", i)
  tmp$region <- rg
  enrich[[rg]] <- tmp
}  

enrich %<>% purrr::reduce(rbind)

enrich %>% write.csv("fig2c_layer_GOTerm.csv")
# enrich <- read.csv("fig2c_layer_GOTerm.csv")

# go term to display
selected_go_terms <- list(
  ARACHNOID = c("GO:0030155", "GO:0006930", "GO:2000145", "GO:0048514", "GO:0072359"),
  
  L1 = c("GO:0099177", "GO:0048172", "GO:0051057", "GO:2000311", "GO:0035249"),
  
  L2 = c("GO:0048168", "GO:0035249", "GO:2000310", "GO:0050804", "GO:1900273"),
  L3 = c("GO:0007416", "GO:0007268", "GO:0046928", "GO:2000145", "GO:0050769"),
  
  L4 = c("GO:0007268", "GO:0010976", "GO:0007409", "GO:0010975", "GO:0098916"),
  
  L5 = c("GO:0021895", "GO:0010975", "GO:2000179", "GO:0051924", "GO:0098916"),
  
  L6 = c("GO:0050804", "GO:0051966", "GO:0050806", "GO:1903169", "GO:0043406"),
  
  WM = c("GO:0042552", "GO:0007409", "GO:0048713", "GO:0010001", "GO:0010975")
)

 print(selected_go_terms)


 print(selected_go_terms)

rg = names(selected_go_terms)[[1]]

for (rg in names(selected_go_terms)) {
  enrich %<>% group_by(region) %>% filter(region != rg | str_detect(Term,str_c(selected_go_terms[[rg]],collapse = "|")))
}
enrichPlot <- enrich[c("group","Description")]

markGenes <- majorMarkers %>%
  filter(str_detect(gene, negate = T, "^A[A-Z].*\\.[0-9]")) %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 5, wt = avg_log2FC) %$% gene
 
markGenes <- c(markGenes, c("CXCL14", "HPCAL1", "NEFM", "NEFH", "SNCG", "NPTX1", "DIRAS2", "ERMN", "BCAS1")) %>% unique

pdf(
  'fig2c.pdf',
  height = 15,
  width = 12,
  onefile = F
)

visCluster(
  object = st.data,
  plot.type = "both",
  column_names_rot = 45,
  show_row_dend = F,
  markGenes = markGenes,
  sample.col = domainColor,
  ctAnno.col = domainColor[rgWithDegs],
  markGenes.side = "left",
  annoTerm.data = enrichPlot,
  line.side = "left",
  cluster.order = st.data$wide.res$cluster %>% unique %>% as.numeric() %>% sort,
  go.col = domainColor[enrichPlot$group %>% str_extract("[0-9]{1,}") %>% as.numeric()],
)

dev.off()


# fig 2d

seu_merge

geneToPlot <- majorMarkers %>%
  filter(str_detect(gene, negate = T, "^A[A-Z].*\\.[0-9]")) %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 20, wt = avg_log2FC)

seu_merge <- ScaleData(seu_merge,features = majorMarkers$gene)

seu_merge@assays$Spatial$scale.data <- LayerData(seu_merge,layer = "scale.data") %>% as.matrix %>% scales::squish(range = c(-3,3))


foo = "L1"
# mclapply(names(domainColor), mc.cores = length(domainColor), function(foo) {
  # print(foo)
  # markerG = geneToPlot %>% filter(cluster == foo) %$% gene
  
  # foo2 = "A00797C3"
  markerG <- c("CXCL14", "HPCAL1", "NEFM", "SCN1B", "RORB-AS1", "TMEM233", "DIRAS2", "ERMN", "BCAS1")
  dataPlot <- FetchData(seu_merge,vars = c(,"x","y","chip"),layer = "data" )
  
  mclapply(seu_merge$chip %>% unique,mc.cores = 10,FUN = function(foo2){
    
    dataPlotTMP <- dataPlot %>% rownames_to_column("cellID") %>% filter(chip == foo2) %>% pivot_longer(names_to = "gene",values_to = "exp",cols = any_of(markerG)) %>% group_by(gene) %>% mutate(scale.exp = scales::oob_squish(scale(exp),c(-2,2))) 
    dataPlotTMP$gene %<>% factor(levels = markerG)
    tmp = ggplot(dataPlotTMP,aes(x,y,color = scale.exp,alpha = scale.exp)) + geom_point(size = 0.5) + facet_wrap(~gene,nrow = 1) + scale_alpha_continuous(range = c(0,1)) + ggplot2::scale_color_gradientn(colors = RColorBrewer::brewer.pal(name = "Spectral",n = 11) %>% rev,limits = c(-2,2)) + coord_fixed() + cowplot::theme_nothing()
    
    ggsave(
      file.path(
        str_c(foo2, "_mainLayer_markers.pdf")
      ),plot =   tmp,
      height = 4,
      width = 3.5 * length(markerG),
      dpi = 300,
      limitsize = FALSE
    )
    
    ggsave(
      file.path(
        str_c(foo2, "_mainLayer_markers.png")
      ),plot =   tmp,
      height = 4,
      width = 3.5 * length(markerG),
      dpi = 200,
      limitsize = FALSE
    )
  })
# })
