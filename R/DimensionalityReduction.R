
#' GetPCA
#'
#' @param object Seurat object.
#' @param dims Dimensions to obtain.
#' @param verbose Print verbose.
#' @param ... Arguments passed to other methods.
#'
#' @return PCA representation.
#' @export
GetPCA <- function(object = NULL, dims = 10, verbose = FALSE, ...){
  object <- Seurat::ScaleData(object, ...)
  object <- Seurat::RunPCA(object, npcs = dims, ...)

  return(object)
}

#' GetUMAP
#'
#' @param object Seurat object.
#' @param dims Dimensions to use.
#' @param reduction Reduction to use.
#' @param PCA Obtain PCA.
#' @param scale Scale the data.
#' @param verbose Print verbose.
#' @param seed Set a random seed. By default, sets the seed to 42.
#' @param ... Arguments passed to other methods.
#'
#' @return UMAP representation.
#' @export
GetUMAP <- function(object = NULL, dims = 10, reduction = "pca", PCA = TRUE, scale = TRUE, seed = 42, verbose = FALSE, ...){
  if(scale)
    object <- Seurat::ScaleData(object, verbose = verbose, ...)
  if(PCA)
    object <- Seurat::RunPCA(object, npcs = dims, verbose = verbose, ...)

  object <- Seurat::RunUMAP(object, reduction = reduction, dims = 1:dims, seed.use = seed, verbose = verbose, ...)

  return(object)
}
