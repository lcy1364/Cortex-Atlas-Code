library(Matrix)
library(Seurat)
library(stringr)
library(magrittr)
# library(SeuratDisk)
library(tidyverse)
library(Matrix)
library(qs)
library(parallel)
devtools::load_all("~/spacexr-master/")
setwd("~/cortex/STEREO/1_Deconvolution_and_QC/")


qsFiles <- list.files("../cellbin/", ".qs$", full.names = T)

snRNA = readRDS("~/cortex/SnRNA/3_mergingDatasets/SnRNA_seurat.RDS")
reference = Reference(GetAssayData(snRNA,slot = "counts"),cell_types = snRNA$subclass,nUMI = snRNA$nCount_RNA)

gtf <- read.csv("../ensemble93gtf_rmXY.csv")
name_id <- gtf$gene_id %>% setNames(gtf$gene_name)

mclapply (qsFiles, mc.cores = length(qsFiles), function(filePath) {
  print(filePath)

  chipName <-
    str_extract(filePath %>% basename(), "[A-Z0-9]*")

  ST <- qs::qread(filePath)

  countSP <- GetAssayData(ST,slot = "count")
  countSP <- countSP[rownames(countSP) %in% gtf$gene_name,]

  rownames(countSP) <- name_id[rownames(countSP)]

  umiSP <-  countSP %>% Matrix::colSums()

  coords <-
    ST@meta.data[c("x", "y")] #%>% dplyr::distinct() %>% column_to_rownames("label")

  puck <- SpatialRNA(coords, countSP, nUMI = umiSP)

  print("3")

  # reference updated
  reference_tmp <- reference
  reference_tmp@counts <-
    reference_tmp@counts[rownames(reference_tmp@counts) %in% rownames(countSP), ]
  reference_tmp@nUMI <- reference_tmp@counts %>% colSums()


  # RCTD
  myRCTD <- create.RCTD(puck, reference_tmp, max_cores = 20)

  qsave(myRCTD, str_c(chipName, "_RCTD.qs"), nthreads = 10)

  NULL
  # myRCTD
})

rctdFiles <- list.files("./","qs$")

mclapply (rctdFiles, mc.cores = 8, function(filePath) {
  #filePath = rctdFiles[[1]]
  myRCTD <- qread(filePath)
  print(filePath)
    print(stringr::str_c(filePath," not done"))
    myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')
    qsave(myRCTD, filePath)
})

