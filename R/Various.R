#' SeuratPreprocessing
#'
#' @param object Seurat object.
#' @param verbose Print verbose.
#' @param ... Arguments passed to other methods.
#'
#' @return Normalized Seurat object and its variable features.
#' @export
SeuratPreprocessing <- function(object = NULL, verbose = FALSE, ...){
  object <- Seurat::NormalizeData(object, verbose = FALSE, ...)
  object <- Seurat::FindVariableFeatures(object, verbose = FALSE, ...)

  return(object)
}
