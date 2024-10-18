library(CellChat)
library(tidyverse)
library(Seurat)
library(qs)
library(magrittr)
options(stringsAsFactors = FALSE)
setwd("~/DATA/BRAIN/STEREO/frequentGraph")

geneIDtoSym <- read.csv("../../SnRNA/gene_kept.csv") %>% {setNames(.$gene_name,nm = .$gene_id)}


scRNA <-  qread("../../SnRNA/3_5_seu.rna.rmXY.qs")
Idents(scRNA) <-"subclass"

cell.list <- WhichCells(scRNA,downsample = 5000,seed = 123)
scRNA <- scRNA[, cell.list]

gc()

scRNA$samples = scRNA$sample
meta <- scRNA@meta.data
exp_normalizated <- Seurat::GetAssayData(object = scRNA,slot = "data")
exp_normalizated@Dimnames[[1]] %>% duplicated() %>% sum()

rownames(exp_normalizated) <- geneIDtoSym[rownames(exp_normalizated)]

exp_normalizated <- exp_normalizated[!duplicated(rownames(exp_normalizated)),]

meta <- scRNA@meta.data

rm(scRNA)
gc()
meta %>%  qsave("cellchatMeta.qs")
exp_normalizated %>% qsave("exp_normalizated.qs")

meta <- qread("cellchatMeta.qs")
exp_normalizated <- qread("exp_normalizated.qs")

cellchat <- CellChat::createCellChat(exp_normalizated, meta  = meta,group.by = "subclass")
rm(exp_normalizated)
gc()

cellchat %>% qsave("cellchat.qs")
 
cellchat <- qread("cellchat.qs")

showDatabaseCategory(CellChatDB.human)

cellchat@DB <- CellChatDB.human

cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

options(future.globals.maxSize = 20000 * 1024 ^ 2)  # Set maximum allowed size to 3000 MiB (3 GiB)

future::plan("multisession", workers = 4) # do parallel
cellchat <- computeCommunProb(cellchat, type = "triMean",nboot = 20)
cellchat %>% qsave("cellchat.qs")

cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat
