

#' GetUMAP
#'
#' @param object Seurat object.
#' @param dims Dimensions to use.
#' @param reduction Reduction to use.
#' @param PCA Obtain PCA.
#' @param scale Scale the data.
#' @param verbose Print verbose.
#' @param ... Arguments passed to other methods.
#'
#' @return
#' @export
#'
#' @examples
GetUMAP <- function(object = NULL, dims = 10, reduction = "pca", PCA = TRUE, scale = TRUE, verbose = FALSE, ...){
  if(scale)
    object <- ScaleData(object, verbose = verbose, ...)
  if(PCA)
    object <- RunPCA(object, npcs = dims, verbose = verbose, ...)

  object <- RunUMAP(object, reduction = reduction, dims = 1:dims, verbose = verbose, ...)

  return(object)
}
