library(Seurat)
library(tidyverse)
setwd("~/cortex/figS1-6/")

Meta <- read_csv("../SnRNA/3_mergingDatasets/SnRNA_Meta.csv")
submeta <- Meta %>% filter(str_detect(class,"exc"))

# figS1b-1
{ggplot(submeta, aes(x = region, fill = subclass)) + geom_bar(position = "fill") + scale_fill_manual(values = subclass_color) + cowplot::theme_cowplot() + theme(axis.text.x = element_text(angle = 45,hjust = 1, vjust = 1))} %>%
ggsave(filename = "figS2b1.pdf",height = 5,width = 8)

# figS1b-2
{ggplot(submeta %>% filter(str_detect(subclass,"IT")), aes(x = region, fill = subclass)) + geom_bar(position = "fill") + scale_fill_manual(values = subclass_color) + cowplot::theme_cowplot() + theme(axis.text.x = element_text(angle = 45,hjust = 1, vjust = 1)) } %>% 
ggsave(filename = "figS2b2.pdf",height = 5,width = 8)

# figS1b-3
{ggplot(submeta %>% filter(!str_detect(subclass,"IT")), aes(x = region, fill = subclass)) + geom_bar(position = "fill") + scale_fill_manual(values = subclass_color) + cowplot::theme_cowplot() + coord_cartesian(expand = F)+ theme(axis.text.x = element_text(angle = 45,hjust = 1, vjust = 1))} %>% 
  {ggsave(filename = "figS2b3.pdf",height = 5,width = 8)}
