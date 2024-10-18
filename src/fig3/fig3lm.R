library(qs)
library(tidyverse)
library(vegan)
library(ggrepel)
library(cowplot)
library(magrittr)
library(scrattch.hicat)
library(ggtree)
library(Seurat)
setwd("~/cortex/fig3/")

sym_id <- readRDS("../SnRNA/1_SnRNA_preprocessing/geneSym_to_geneID.RDS")
id_sym <- read_csv("../SnRNA/1_SnRNA_preprocessing/gene_kept.csv") %>% {setNames(object = .$gene_uni,nm = .$gene_id)}
sub.seu.rna.ETall <- readRDS("../SnRNA/SnRNA_seurat.RDS") %>% subset(subclass == "ET")

sub.seu.rna.ETall$ACC <- sub.seu.rna.ETall$region == "ACC"
sub.seu.rna.ET <- sub.seu.rna.ETall %>% SplitObject("ACC")

sub.seu.rna.ET$`FALSE`$region_cluster <-  str_c("others",sub.seu.rna.ET$`FALSE`$cluster,sep = "_")
sub.seu.rna.ET$`TRUE`$region_cluster <-  str_c("ACC",sub.seu.rna.ET$`TRUE`$cluster,sep = "_")

region_color <- c(FPPFC = "#3F4587", DLPFC = "#8562AA", VLPFC = "#EC8561", M1 = "#B97CB5", 
                  S1 = "#D43046", S1E = "#F0592B", PoCG = "#ED4A96", SPL = "#593C97", 
                  SMG = "#A54486", AG = "#FBDE13", V1 = "#299FAE", ITG = "#75CCE3", 
                  STG = "#0C6939", ACC = "#0D9547")


VEN_markers <-
  c("ADRA1A",
    "TCIM",
    "GABRQ",
    "VMAT2",
    "HTR2B",
    "DRD3",
    "DISC1",
    "FEZF2",
    "CTIP2",
    "SLC40A1",
    "ITGA4",
    "ACTG1P18",
    "SFRP2",
    "DSG2",
    "BMP3",
    "OCAR",
    "RREB1",
    "TEC",
    "TIPARP",
    "GABRQ",
    "VIM",
    "LINC01320",
    "NEUROD1",
    "DNAI1",
    "ANKRD34B",
    "ANGPT1",
    "POU3F1",
    "GYG2",
    "HAPLN4",
    "ADAMTS2",
    "EPHX2",
    "LINC00861",
    "TXK",
    "CTXN3",
    "LYPD1",
    "DHRS9",
    "FAM84B",
    "GNA15"
  ) %>% unique

venGenes <- intersect(sym_id[VEN_markers],rownames(sub.seu.rna.ETall)) %>% keep(~!is.na(.))
ETplot <- list()

for (cl in names(sub.seu.rna.ET)) {

  ET.cl_means <- get_cl_means(sub.seu.rna.ET[[cl]]@assays$RNA@counts[sub.seu.rna.ET[[cl]]@assays$RNA@var.features,] %>% logCPM,
                              setNames(sub.seu.rna.ET[[cl]]$region_cluster, colnames(sub.seu.rna.ET[[cl]])))
  if(ncol(ET.cl_means) > 1) {ET_tree <- build_dend(ET.cl_means) %$%  ggtree(pvclust.result) + scale_x_reverse() + geom_tiplab()
  labelord <- get_taxa_name(ET_tree)}else{
    labelord <- colnames(ET.cl_means)
  }
  Idents(sub.seu.rna.ET[[cl]]) = "region_cluster"
  
  
  tmp <- VlnPlot(sub.seu.rna.ET[[cl]], venGenes, stack = T) + theme(legend.position = "none", axis.text.x = element_text( vjust = 1, hjust = 1))
  tmp$data$feature <- id_sym[tmp$data$feature %>% as.character()]
  tmp$data$feature %<>% factor(levels = VEN_markers)
  tmp$data$ident %<>% as.character() %>% factor(levels = rev(labelord))
  if(ncol(ET.cl_means) > 1) {
  ETplot[[cl]] <- tmp + ET_tree + patchwork::plot_layout(widths = c(0.5,0.1))
  } else {
    ETplot[[cl]] <- tmp
  }
}

ggsave("fig3l.pdf",ETplot[[1]] + ETplot[[2]],width = 15,height = 5)

# fig3m
Idents(sub.seu.rna.ETall) <- "region"
ET.rg.markers <- parallel::mclapply(sub.seu.rna.ETall$region %>% unique,mc.cores = 14,FUN = function(x){
  tmp <- FindMarkers(sub.seu.rna.ETall,ident.1 = x,ident.2 = NULL,group.by = "region")
  tmp$region <- x
  tmp$gene <- rownames(tmp)
  tmp
})

ET.rg.markers %<>% keep(~is.data.frame(.))  %>% purrr::reduce(rbind)
ET.rg.markers %<>% filter(p_val_adj <= 0.05)

ET.rg.markers %>% write_csv("./fig3m_et_deg.csv")

library(org.Hs.eg.db)
library(clusterProfiler)
library(DOSE)
data(geneList)

regions <- ET.rg.markers$region %>% unique
enrichment_rg_res <- list()

for(rg in regions){
  print(rg)
  timestamp()
  rg.ET.rg.markers <- ET.rg.markers %>% filter(region == rg,avg_log2FC >0)
  Go.rg <- enrichR::enrichr(id_sym[rg.ET.rg.markers$gene],databases = "GO_Biological_Process_2021")[[1]]
  Go.rg$region = rg
  Go.rg %<>% filter(Adjusted.P.value <0.05)
  if(nrow(Go.rg) >0 ) enrichment_rg_res[[rg]] <- Go.rg
}

enrichment_rg_res %<>% purrr::reduce(rbind)

enrichment_rg_res %>% write.csv("./fig3m_et_deg_go.csv")

enrichment_rg_res.selected <- enrichment_rg_res %>% group_by(region) %>%  slice_min(order_by = Adjusted.P.value,n = 5)

enrichment_rg_res.selected$y = str_c(enrichment_rg_res.selected$Term,enrichment_rg_res.selected$region, sep  = " ")

enrichment_rg_res.selected$x <- purrr::map(1:nrow(enrichment_rg_res.selected),~eval(parse(text = enrichment_rg_res.selected[.,"Overlap"]))) %>% unlist

enrichment_rg_res.selected$logPadj <- -log(enrichment_rg_res.selected$Adjusted.P.value,base = 10)

ven_go_term = enrichment_rg_res.selected%>% arrange(region,x, desc(logPadj))
ven_go <- ggplot(ven_go_term , aes(x,y, size = logPadj, color = logPadj )) +
  annotation_raster(rev(alpha(region_color[ven_go_term$region],0.7)),xmin = min(ven_go_term$x)*0.7,xmax = Inf, ymin = Inf,ymax = -Inf) +
  geom_point() + scale_y_discrete(limits = ven_go_term$y, labels = str_wrap(ven_go_term$Term, width = 30)) +
  geom_text(aes(label = region,x = 1.1*max(x)),color =  region_color[ven_go_term$region], size = 5, angle = 90) +
  scale_color_viridis_c(space = "Lab") +
  coord_cartesian(clip = "off",expand = F, xlim = c(min(ven_go_term$x)*0.7,max(ven_go_term$x)*1.1)) +
  guides(size = guide_legend(title = "logPadj"), color = guide_legend(title = "logPadj")) + cowplot::theme_cowplot() + labs(title = "GO Term") +
  theme(axis.text.y = element_text(colour = region_color[ven_go_term$region],lineheight = 0.7))

ggsave("fig3m_exc_ven_go.pdf",
       ven_go,
       height = 10,
       width = 8)

