library(Seurat)
library(parallel)
library(stringr)
library(tidyverse)
library(magrittr)
library(hdf5r)
library(qs)

setwd("~/cortex/STEREO/GEM")
h5adFiles <- list.files("~/cortex/STEREO/GEM/bin200", "\\.h5ad", full.names = T)

file = h5adFiles[[1]]
mclapply(h5adFiles, mc.cores = 5, function (file) {
  print(file)
  data <- H5File$new(file)
  # get cell index
  meta.data = data.frame(
    cellId = data[["obs"]][["_index"]]$read(),
    x = data[["obs"]][["x"]]$read(),
    y = data[["obs"]][["y"]]$read()
  )
  
  print(nrow(meta.data))
  rownames(meta.data) <- meta.data$cellId
  
  for (cl in names(data[["obs"]])[str_detect(names(data[["obs"]]), "domain")]) {
    meta.data[[cl]] = data[["obs"]][[cl]][["codes"]]$read() %>% factor(levels = data[["obs"]][[cl]][["categories"]]$read())
  }
  
  # create Seurat object
  seu <- CreateSeuratObject(
    assay = "Spatial",
    counts =  Matrix::sparseMatrix(
      i = data[["X"]][["indices"]]$read() + 1,
      x = data[["X"]][["data"]]$read(),
      p = data[["X"]][["indptr"]]$read(),
      dimnames = list(data[["var"]][["_index"]]$read(), meta.data$cellId)
    ),
    meta.data = meta.data
  )
  
  cell_coords <- unique(seu@meta.data[, c('x', 'y')])
  cell_coords['cells'] <- row.names(cell_coords)
  
  seu@images[['image']] <- new(Class = 'SlideSeq', coordinates = cell_coords)
  seu@images$image@key <- "image_"
  seu@images$image@assay <- "RNA"
  
  data$close_all()

  qsave(seu, file.path("bin200",str_replace(basename(file), "h5ad", "qs")))
  
})
