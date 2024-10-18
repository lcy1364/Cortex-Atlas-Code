library(qs)
library(tidyverse)
library(vegan)
library(ggrepel)
library(cowplot)
library(ggh4x)
library(magrittr)
setwd("~/data/STEREO//AnalysisPlot/")

output_n = "."


region_color <- c(FPPFC = "#3F4587", DLPFC = "#8562AA", VLPFC = "#EC8561", M1 = "#B97CB5", 
                  S1 = "#D43046", S1E = "#F0592B", PoCG = "#ED4A96", SPL = "#593C97", 
                  SMG = "#A54486", AG = "#FBDE13", V1 = "#299FAE", ITG = "#75CCE3", 
                  STG = "#0C6939", ACC = "#0D9547")

subclass_color <-
  c(
    AST = "#665C47",
    ENDO = "#604B47",
    ET = "#CEC823",
    `L2-L3 IT LINC00507` = "#07D8D8",
    `L3-L4 IT RORB` = "#09B2B2",
    `L4-L5 IT RORB` = "#69B199",
    `L6 CAR3` = "#898325",
    `L6 CT` = "#2D8CB8",
    `L6 IT` = "#19E3BE",
    L6b = "#7944AA",
    LAMP5 = "#D96C07",
    MICRO = "#513577",
    NDNF = "#f2798d",
    NP = "#93C43B",
    OLIGO = "#5E8A79",
    OPC = "#2E3E39",
    PAX6 = "#E96DF2",
    PVALB = "#FF2D4E",
    SST = "#F2B003",
    VIP = "#9C28CC",
    VLMC = "#697255"
  )

scRNA <- qread("../../SnRNA/merge_seu.qs")

# IT relationship with region ----------------------------------------------
cell_por_per_region <- scRNA@meta.data %>% filter(str_detect(subclass,"IT")) %>% group_by(sample) %>%
  slice_sample(n = 20000, replace = T) %>%  group_by(region,subclass) %>% summarise(count = n()) %>%
  pivot_wider(names_from = subclass, values_from = count, values_fill = 0) %>% column_to_rownames("region") %>% as.data.frame()

cell_por_per_region %<>% scale

subclass.rg <- fortify(vegan::rda(cell_por_per_region),scaling = 1,axes = 1:2)

t3 <- ggplot(data = subclass.rg, mapping = aes(PC1,PC2,color = label,shape = score, label = str_wrap(label,10))) + geom_point() +
  geom_label_repel(data = subclass.rg %>% filter(score == "sites")) +
  geom_text_repel(data = subclass.rg %>% filter(score == "species"),lineheight = .7,size = 3) +
  # geom_point(data = fc.rg %>% filter(score = "sites"),mapping = aes( color = score)) +
  scale_colour_manual(values = c(region_color,subclass_color)) +
  coord_axes_inside() + theme_cowplot() + theme(legend.position = "none")

ggsave("10_fmri/3_region_IT_PCA.png",t3,height = 5,width = 5)
ggsave("10_fmri/3_region__IT_PCA.pdf",t3,height = 5,width = 5)


data <- qread("./9_Layers_analysis/cutData.qs")
cellMeta <- purrr::map_dfr(data, ~ .$cells)

 cell_por_per_region <- cellMeta %>% filter(major == "EXC",withinRegion == TRUE)  %>% group_by(chip) %>%   
  slice_sample(n = 20000, replace = T) %>%  group_by(region.x) %>% slice_sample(n = 20000, replace = T) %>%   group_by(region.x,subclass) %>% summarise(count = n()) %>%
  pivot_wider(names_from = subclass, values_from = count, values_fill = 0) %>% column_to_rownames("region.x") %>% as.data.frame()
 
cell_por_per_region %<>% scale
PCA_res <- vegan::rda(cell_por_per_region)
subclass.rg <- fortify(PCA_res,scaling = 1,axes = 1:2)

EigenPer <- PCA_res$CA$eig/sum(PCA_res$CA$eig)
 
 
t4 <- ggplot(data = subclass.rg, mapping = aes(PC1,PC2,color = label,shape = score, label = str_wrap(label,10))) + geom_point() +
  geom_label_repel(data = subclass.rg %>% filter(score == "sites")) +
  geom_text_repel(data = subclass.rg %>% filter(score == "species"),lineheight = .7,size = 3) +
  # geom_point(data = fc.rg %>% filter(score = "sites"),mapping = aes( color = score)) +
  scale_colour_manual(values = c(region_color,subclass_color)) +
  coord_axes_inside() + theme_cowplot() + theme(legend.position = "none") + xlab(str_c("PC1(",round(EigenPer[[1]]*100,digits = 0),")%")) + 
  ylab(str_c("PC2(",round(EigenPer[[2]]*100,digits = 0),")%"))

ggsave("fig6c.png",height = 5,width = 5)


