library(qs)
library(tidyverse)
library(vegan)
library(ggrepel)
library(cowplot)
library(ggh4x)
library(magrittr)
library(vegan)
library(ggrepel)
setwd("~/cortex/fig3/")

taskScore <- data.frame(
  stringsAsFactors = FALSE,
  row.names = c("working memory",
                "visuospatial","visual semantics","visual perception",
                "visual attention","verbal semantics","social cognition",
                "reward-based decision making","reading","pain",
                "numerical cognition","multisensory processing","motor","language",
                "inhibition","face/affective processing","eye movements",
                "emotion","declarative memory","cued attention",
                "cognitive control","autobiographical memory",
                "auditory processing","action"),
  SMG = c(0.518313938,0.510914088,
          0.497169209,0.499690174,0.506875651,0.501363047,
          0.485649048,0.466527368,0.492202943,0.524847636,0.541496893,
          0.491313586,0.507270723,0.492845959,0.495510931,
          0.466863174,0.518483707,0.475281493,0.485610868,0.520856037,
          1,0.515259099,0.514299939,0.567425482),
  PoCG = c(0.497949242,0.532038227,
           0.480526295,0.508920598,0.479202828,0.417571097,
           0.481624112,0.445064213,0.473181551,0.655965967,0.488123918,
           0.508229237,0.545311045,0.475208897,0.501377892,
           0.451495143,0.533199371,0.432125953,0.395147589,0.470048689,
           1,0.422147502,0.452451623,0.514699401),
  VLPFC = c(0.533362404,0.474569025,
            0.463825896,0.458791849,0.507912825,0.596526898,
            0.450152788,0.44802403,0.494840937,0.582321484,0.511315185,
            0.499593028,0.554827778,0.607584642,0.47279679,
            0.465815312,0.473676902,0.471063774,0.477135959,0.497979148,1,
            0.414852378,0.595775623,0.541691893),
  S1 = c(0.483500331,0.526135797,
         0.480974424,0.496902216,0.484905492,0.451256846,
         0.478074467,0.467696436,0.486552996,0.606382918,0.455908366,
         0.520986638,0.56782429,0.495615798,0.486922626,
         0.452157062,0.518374409,0.45529447,0.448931976,0.518433029,
         0.46875,0.466393492,0.479039323,0.57581881),
  V1 = c(0.522821458,0.580000528,
         0.528787433,0.561330772,0.48928388,0.497757901,0.473012636,
         0.494712634,0.534649611,0.38415173,0.511669859,
         0.553659951,0.500973627,0.521263073,0.477253974,
         0.529143601,0.562406883,0.506879235,0.504876481,0.528327412,
         0.273822563,0.524034953,0.503767272,0.526653828),
  FPPFC = c(0.512934069,0.455359044,
            0.459534961,0.424764974,0.493669865,0.460619896,
            0.527814706,0.528685071,0.496322947,0.505213519,0.46290408,
            0.462737098,0.481837186,0.469187832,0.508242241,
            0.502925239,0.455689296,0.522338526,0.519006209,0.517724272,
            0.534627983,0.563934176,0.463105947,0.483666536),
  SPL = c(0.541317688,0.54839282,
          0.495512054,0.519316612,0.507901299,0.48370745,0.495784628,
          0.471199509,0.487081227,0.447724128,0.572504259,
          0.507105492,0.487127414,0.505521268,0.500900966,
          0.468898121,0.545266989,0.455200016,0.508955978,0.523155302,
          0.617940837,0.511718632,0.438807255,0.537849554),
  S1E = c(0.486491674,0.456294013,
          0.435580609,0.420710627,0.446319132,0.517211704,
          0.393683461,0.404010108,0.470286079,0.709349822,0.540403,
          0.54003684,0.568304421,0.538708294,0.394667227,0.400259581,
          0.507163557,0.419137978,0.44724335,0.428803107,1,
          0.273116489,0.549959111,0.461563714),
  M1 = c(0.512461566,0.500796339,
         0.475289088,0.493193754,0.482601497,0.505950496,
         0.476183773,0.475736419,0.486286831,0.52482022,0.508195844,
         0.496747398,0.534266723,0.510218848,0.503816833,
         0.465954521,0.546629654,0.468563435,0.477166035,0.513966285,
         0.428009969,0.495970409,0.502920738,0.554195549),
  STG = c(0.483380878,0.467266514,
          0.482031202,0.489623736,0.488628497,0.611200253,
          0.488346124,0.42738236,0.493503171,0.48905351,0.44607305,
          0.520180683,0.484397738,0.572313781,0.466673162,
          0.501392136,0.514525467,0.486324657,0.473554279,0.466839644,1,
          0.503846036,0.681902701,0.519516384),
  DLPFC = c(0.572046938,0.503781361,
            0.513700839,0.460227897,0.511025045,0.517471132,
            0.474069878,0.477761888,0.511540073,0.460330416,0.559381176,
            0.494171282,0.482098923,0.562140494,0.505588453,
            0.505118568,0.449207689,0.515213272,0.529620695,0.524232792,
            0.34936212,0.512648864,0.484225036,0.511312334),
  AG = c(0.503906824,0.515662149,
         0.531070918,0.530041846,0.486855411,0.571882957,
         0.530104171,0.471418453,0.459086664,0.43718692,0.53106337,
         0.500513385,0.460172783,0.527959973,0.470011306,
         0.506083313,0.535090033,0.46977407,0.532473163,0.483616803,1,
         0.567996029,0.481908955,0.536479353),
  ACC = c(0.490878691,0.438283994,
          0.46922636,0.417305497,0.462780466,0.416991096,0.488056229,
          0.5067212,0.489735828,0.586335954,0.430990977,
          0.483902241,0.506215864,0.481713226,0.49161774,0.462221564,
          0.511542583,0.4847342,0.478152261,0.507752421,1,
          0.539212532,0.457047935,0.519477054),
  ITG = c(0.487130452,0.528877414,
          0.585129523,0.455886355,0.479057031,0.604121293,
          0.486366929,0.451849336,0.445995738,0.37653077,0.450740875,
          0.460620745,0.403343937,0.538995798,0.380312895,
          0.522103682,0.491534767,0.47031267,0.560941396,0.424116644,1,
          0.576928733,0.410022811,0.486537806)
) %>% t %>% as.data.frame()

taskScore <- taskScore[c("working memory",
                         "visual semantics","visual perception",
                         "verbal semantics","social cognition",
                         "reward-based decision making","reading","pain",
                         "numerical cognition","multisensory processing","motor","language",
                         "autobiographical memory","auditory processing","action")]

taskColor <- c(ggsci::pal_aaas()(10), ggsci::pal_nejm()(8),ggsci::pal_npg()(10)) %>% unique %>% sort %>%  setNames( str_replace_all(colnames(taskScore),"\\."," "))

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

scRNA <- read.csv("../SnRNA/3_mergingDatasets/SnRNA_Meta.csv")

cell_por_per_region <- Meta %>% filter(str_detect(class, "exc")) %>% group_by(library_prep) %>% sample_n(10000, replace = T) %>%
  group_by(region,subclass,library_prep) %>% summarise(count = n()) %>% 
  pivot_wider(names_from = subclass, values_from = count, values_fill = 0) %>% data.frame()


errorNameCor <- c( "ET", "L2-L3 IT LINC00507", "L3-L4 IT RORB", 
                   "L4-L5 IT RORB", "L6 CAR3", "L6 CT", "L6 IT", "L6b", "NP") %>% setNames(c( "ET", "L2.L3.IT.LINC00507", "L3.L4.IT.RORB", 
                                                                                              "L4.L5.IT.RORB", "L6.CAR3", "L6.CT", "L6.IT", "L6b", "NP"))

rownames(cell_por_per_region) <- cell_por_per_region$region %>% make.unique()

decorana(cell_por_per_region[-c(1:2)])


cell_por_per_region.hel <- decostand(cell_por_per_region[-c(1:2)],"hellinger")

rda.res <- vegan::rda(Y = cell_por_per_region.hel, X = taskScore[cell_por_per_region$region,],scale = T)

rda.pos <- ggvegan:::fortify.cca(rda.res,scaling = 2 ,axes = 1:2) %>%  filter(!score %in% c("constraints","sites")) 

rda.pos$label_right <- errorNameCor[rda.pos$label] 
rda.pos$label_right[rda.pos$label_right %>% is.na()] <- rda.pos$label[rda.pos$label_right %>% is.na()]

rg.subclass.function.plt <- ggplot() +
  geom_segment(data = rda.pos %>% filter(score == "biplot"),mapping = aes(x = 0, xend = RDA1, y =  0, yend = RDA2, color = score),color = "grey30", alpha = 0.5) +
  geom_point(data = rda.pos %>% filter(score == "biplot"),mapping = aes(x = RDA1, y = RDA2, color = label_right)) +
  geom_text_repel(data = rda.pos %>% filter(score %in% c("species","biplot")),mapping = aes(x = RDA1,y = RDA2,  label  = label_right %>% str_wrap(10), color = label_right),lineheight = .6,max.overlaps = 70,size = 3) +
  geom_point(data = rda.pos %>% filter(score == "species"), mapping = aes(x = RDA1,  y = RDA2,color = label_right ), shape = 21) +
  cowplot::theme_cowplot() + theme(legend.position = "none") + 
  scale_color_manual(values = c(subclass_color,taskColor),na.value = "grey20") +  
  xlab( str_c("RDA 1 (",round(100*rda.res$CCA$eig[1]/sum(rda.res$CCA$eig),digits = 2),"%)")) + 
  ylab( str_c("RDA 2 (",round(100*rda.res$CCA$eig[2]/sum(rda.res$CCA$eig),digits = 2),"%)")) + coord_axes_inside() 

rg.subclass.function.plt

ggsave("fig3j.pdf",rg.subclass.function.plt,height = 5,width =5)
