library(Seurat)
library(parallel)
library(stringr)
library(tidyverse)
library(magrittr)
library(hdf5r)
library(dplyr)
library(rjson)
library(Seurat)
library(ggplot2)
library(argparser)
library(SeuratDisk)


setwd("~/cortex/STEREO/1_cellbin/")
h5adFiles <- list.files("./batch_gem/", ".h5ad$", full.names = T)

mclapply(h5adFiles, mc.cores = min(length(h5adFiles), 50), function (file) {
  infile <- file
  chipID <- str_extract(basename(file), "[A-Z0-9]*")
  outfile <- str_c(chipID, ".qs")
  
  # convert h5ad as h5seurat, which means a seurat-object format stored in h5
  Convert(infile,
          dest = "h5seurat",
          assay = "Spatial",
          overwrite = TRUE)
  
  h5file <- paste(paste(unlist(strsplit(
    infile, "h5ad", fixed = TRUE
  )), collapse = 'h5ad'), "h5seurat", sep = "")
  print(paste(
    c("Finished! Converting h5ad to h5seurat file at:", h5file),
    sep = " ",
    collapse = NULL
  ))
  
  object <- LoadH5Seurat(h5file, assays = "Spatial")
  print(paste(
    c("Successfully load h5seurat:", h5file),
    sep = " ",
    collapse = NULL
  ))
  
  # spatial already transform to `Spatial`` in assays
  if (!is.null(object@reductions$spatial)) {
    object@reductions$spatial <- NULL
  }
  
  assay.used <- 'Spatial'
  
  
  # TODO follow with old code, don't touch
  print("Start add image...This may take some minutes...(~.~)")
  # add image
  cell_coords <- unique(object@meta.data[, c('x', 'y')])
  cell_coords['cells'] <- row.names(cell_coords)
  
  # generate object @images$slice1
  generate_BGI_spatial <- function(cell_coords) {
    return(new(Class = 'SlideSeq', coordinates = cell_coords))
  }
  
  BGI_spatial <- generate_BGI_spatial(cell_coords = cell_coords)
  
  # can be thought of as a background of spatial
  # import image into seurat object
  object@images[['image']] <- BGI_spatial
  object@images$image@key <- "image_"
  object@images$image@assay <- "Spatial"
  
  # conversion done, save
  print("Finished add image...Start to saveRDS...")
  
  qsave(object, outfile)
})
