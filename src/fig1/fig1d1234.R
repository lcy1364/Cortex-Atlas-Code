library(tidyverse)
library(Seurat)
library(pvclust)
library(ggtree)
library(magrittr)
library(dendextend)
library(patchwork)
setwd("~/cortex/fig1/")

subclass_color <-c(AST = "#665C47", ENDO = "#604B47", ET = "#CEC823", CHANDELIER = "#E25691", 
                   `L2-L3 IT LINC00507` = "#07D8D8", `L3-L4 IT RORB` = "#09B2B2", 
                   `L4-L5 IT RORB` = "#69B199", `L6 CAR3` = "#898325", `L6 CT` = "#2D8CB8", 
                   `L6 IT` = "#19E3BE", L6B = "#7944AA", L6b = "#7944AA",LAMP5 = "#D96C07", MICRO = "#513577", 
                   NDNF = "#f2798d", NP = "#93C43B", OLIGO = "#5E8A79", OPC = "#2E3E39", 
                   PAX6 = "#E96DF2", PVALB = "#FF2D4E", SST = "#F2B003", VIP = "#9C28CC", 
                   VLMC = "#697255")

region_color <- c(FPPFC = "#3F4587", DLPFC = "#8562AA", VLPFC = "#EC8561", M1 = "#B97CB5", 
                  S1 = "#D43046", S1E = "#F0592B", PoCG = "#ED4A96", SPL = "#593C97", 
                  SMG = "#A54486", AG = "#FBDE13", V1 = "#299FAE", ITG = "#75CCE3", 
                  STG = "#0C6939", ACC = "#0D9547")


dend <- readRDS("../SnRNA/dend.RDS")
p1 <- dend
ggsave("./fig1d_dend_1.pdf",dend,width = 30 ,height = 5)


meta <- read.csv("../SnRNA/3_mergingDatasets/SnRNA_Meta.csv")

# fig1d2 sample ratio

meta$donor %>% unique
meta$cluster %<>% as.character() %>% factor(levels = rev(get_taxa_name(dend) %>% str_to_upper()))
pdf("./fig1d_dend_2.pdf",width = 25,height = 3)
p2 <- ggplot(meta,aes(x = cluster,fill = donor)) + geom_bar(position = "fill") + 
  scale_fill_manual(values = MetBrewer::met.brewer(name = "Redon",n = 10)) +  cowplot::theme_nothing() + 
  theme(axis.text.x = element_blank(),axis.title = element_blank())
dev.off()


# fig1d3 regional distribution

pdf("./fig1d_dend_3.pdf",width = 15,height = 3)
p3 <- ggplot(meta,aes(x = cluster,fill = region)) + geom_bar(position = "fill") + scale_fill_manual(values = region_color) + cowplot::theme_nothing() + 
  theme( axis.text.y = element_text(hjust = 1,vjust = 0.5), axis.text.x = element_blank(), axis.title = element_blank())
dev.off()


# fig1d4

spatialCellMeta <- read_csv("../STEREO/spatialCellMeta.csv")
spatialCellMeta$cluster %<>% as.character() %>% str_to_upper %>% factor(levels = get_taxa_name(dend) %>% str_to_upper() %>% rev)

pdf("fig1d_dend_4.pdf",height = 5,width = 25)
p4 <- ggplot(spatialCellMeta, aes(y = cluster,  x = depth, fill = subclass )) +  ggridges::geom_density_ridges(scale = 3,bandwidth = 0.5) +
  ggridges::theme_ridges() + cowplot::theme_minimal_hgrid()+ 
  theme(legend.position = "none", axis.text.y = element_text(hjust = 1,vjust = 0.5), axis.text.x = element_blank(),axis.title = element_blank()) + 
  scale_fill_manual(values = subclass_color) + 
  coord_flip(xlim = c(-0.5,8),expand = F)  + 
  scale_x_continuous(breaks = 0:7,labels = c("ARACHNOID", "L1", "L2", "L3", "L4", "L5", "L6", "WM")) 

dev.off()

ggsave("fig1d.pdf",p1/p3/p4 + patchwork::plot_layout(heights = c(2,1,0.7)),height = 8,width = 20)
