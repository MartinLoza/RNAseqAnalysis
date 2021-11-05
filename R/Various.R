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

#' Title
#'
#' @param object Object to sample
#' @param frac Fraction of samples. Must be a value between 0 and 1.
#' @param seed Seed to use on sampling
#' @param ... Arguments passed to other methods.
#'
#' @return Sampled data
#' @export
SampleData <- function(object = NULL, frac = NULL, seed = 777, ...){
  set.seed(seed)
  nCells <- ncol(object)
  samples <- floor(frac*nCells)
  idx <- sample(x = nCells, size = samples, ... )

  return(object[,idx])
}

SortDF <- function(df = NULL, sort.by = NULL, decreasing = TRUE){
  idx <- sort(df[[sort.by]], decreasing = decreasing, index.return = TRUE)$ix
  tmp <- df[idx,]
  return(tmp)
}
