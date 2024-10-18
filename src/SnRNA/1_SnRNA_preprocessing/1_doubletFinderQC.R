library(Seurat)
library(tidyverse)
library(ggplot2)
library(parallel)
library(magrittr)
library(DoubletFinder)
library(rtracklayer)
library(scater)
setwd("~/cortex/SnRNA/1_SnRNA_preprocessing")
matrixFiles <- list.files("0_matrix/","filtered_feature_bc_matrix$",recursive = T,include.dirs = T,full.names = T)
out_dir <- "./QC"
set.seed(1234)
dir.create(out_dir)

gtf <- import("../1_SnRNA_preprocessing/gencode.v32.primary_assembly.annotation.gtf") %>% as.data.frame() %>% 
  select(seqnames,gene_name,gene_id) %>% mutate(gene_id = str_remove(gene_id,"\\..*")) %>% distinct() %>% filter(str_detect(seqnames,"chr[0-9]*$")) %>% {setNames(.$gene_id,nm = .$gene_name)} 

gtf %>% saveRDS("../1_SnRNA_preprocessing/geneSym_to_geneID.RDS")

matrixFile = matrixFiles[[1]] 

while (length(list.files("./QC","qs")) != 42) {
mclapply(matrixFiles,mc.cores = floor(length(matrixFiles)/2),FUN = function(matrixFile){
  sample <- matrixFile %>% str_extract("S[0-9]{4}-A[0-9]{1,2}")
  outputFile <- file.path(out_dir, str_c(sample, '.qs'))
  if (!file.exists(outputFile)) {

  counts <- Read10X(matrixFile)
  colnames(counts) <- str_c(sample, colnames(counts), sep = "_")
  seu <- CreateSeuratObject(counts)
  seu$donor <- str_extract(sample,"^S[0-9]{4}")

  seu <-
    PercentageFeatureSet(seu, "^MT-", col.name = "percent.mt")

  seu <-
      PercentageFeatureSet(seu, "^RP[SL]", col.name = "percent.ribo")

  seu <-
      PercentageFeatureSet(seu, "^HB[^(P)]", col.name = "percent.hb")
  
  seu <-
      PercentageFeatureSet(seu, "PECAM1|PF4", col.name = "percent.plat")

  seu <-
      subset(seu, subset = nFeature_RNA > 500 &
               percent.mt < 1 & percent.ribo < 1)

    seu <-
      NormalizeData(seu,
                    normalization.method = "LogNormalize",
                    scale.factor = 10000)
    seu <-
      FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
    seu <- ScaleData(seu, features = rownames(seu))
    seu <- RunPCA(seu, features = VariableFeatures(object = seu))
    seu <- FindNeighbors(seu, dims = 1:15)
    seu <- RunUMAP(seu, dims = 1:15)
    seu <- FindClusters(seu, dims = 1:15)
    
    ################### Doublet ###########################
    
    # system.time({
    sweep.res.list <- paramSweep(seu, PCs = 1:15,num.cores = 5)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    mpK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
    # setting the appropriate nEXP(by 10x reagents)
    annotations <- seu$seurat_clusters
    homotypic.prop <- modelHomotypic(annotations)
    nExp_poi <- round(1/100 *(ncol(seu)/1000) * ncol(seu@assays$RNA@data))
    # https://kb.10xgenomics.com/hc/en-us/articles/360054599512-What-is-the-cell-multiplet-rate-when-using-the-3-CellPlex-Kit-for-Cell-Multiplexing-
    # 0.8% per 1000 cell
    nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
    
    seu <- doubletFinder(
      seu,
      PCs = 1:15,
      pN = 0.25,
      pK = mpK,
      nExp = nExp_poi,
      reuse.pANN = FALSE)
    
    
    # 3 MAD
    
    seu.nCount <-
      isOutlier(seu$nCount_RNA,
                nmads = 3,
                type = "lower",
                log = TRUE)
    
    seu.nFeatures <-  isOutlier(seu$nFeature_RNA,
                                nmads = 3,
                                type = "lower",
                                log = TRUE)
    
    seu.mito <- isOutlier(seu$percent.mt,
                          nmads = 3,
                          type = "higher",
                          log = TRUE)
    
    seu.doublet <-
      seu@meta.data[, str_detect(colnames(seu@meta.data), "DF.classifications")] == "Doublet"
    
    cells.out <- seu.mito | seu.nFeatures | seu.nCount | seu.doublet
    
    print(table(cells.out))
    
    
    counts <- Read10X()
    
    seu <- seu[, !cells.out]
    
    exp <- GetAssayData(object = seu, assay = "RNA",layer = "counts")
    
    rownames(exp) <- gtf[rownames(exp)]
    
    qs::qsave(x = CreateSeuratObject(exp[!is.na(rownames(exp)),],meta.data = seu@meta.data),outputFile,preset = "fast")
  }
  }) }


resList <- mclapply(list.files("./QC","qs",full.names = T),mc.cores = length(list.files("./QC","qs")),FUN = qs::qread)
seu <- merge(resList[[1]],resList[-c(1)]) 
seu %>% qs::qsave("./SnRNA_seurat.qs",preset = "fast")

