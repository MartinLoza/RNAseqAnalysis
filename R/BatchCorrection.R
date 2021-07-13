#' RunComBatseq
#'
#' @param object A seurat object to correct batch effects.
#' @param batch Batch labels.
#' @param runningTime Return the running time.
#' @param verbose Print verbose.
#' @param ... Arguments passed to other methods.
#'
#' @return Corrected and normalized Seurat object.
#' @export
RunComBatseq <- function(object = NULL, batch = "batch", runningTime = FALSE, verbose = FALSE, ...){

  features <- Seurat::VariableFeatures(object)
  if(length(features) == 0){
    warning("Variable features not defined. Running 'FindVariableFeatures' function.", call. = TRUE)
    features <- Seurat::VariableFeatures(Seurat::FindVariableFeatures(object, verbose = verbose))
  }

  counts <- as.matrix(Seurat::GetAssayData(object, assay = "RNA", slot = "counts")[features,])
  md <- object[[]]

  if(!(batch %in% colnames(md)))
    stop(paste0(batch, "not found in object's metadata. Check the batch label."))

  time <- system.time({
    corrCounts <- sva::ComBat_seq(counts = counts, batch = md[[batch]], full_mod = FALSE)
  })

  object[["integrated"]] <- Seurat::CreateAssayObject(counts = corrCounts)
  Seurat::DefaultAssay(object) <- "integrated"

  object <- Seurat::NormalizeData(object = object, assay = "integrated", verbose = verbose, ...)
  Seurat::VariableFeatures(object) <- features

  if(runningTime == FALSE)
    return(object)
  else
    return(list(object = object, time = time))
}

#' RunScMerge
#'
#' @param object A seurat object to correct batch effects.
#' @param batch Batch labels.
#' @param ks A vector indicates the kmeans's K for each batch, which length needs to be the same as the number of batches.
#' @param runningTime Return the running time.
#' @param verbose Print verbose.
#' @param ... Arguments passed to other methods.
#'
#' @return Corrected Seurat object.
#' @export
RunScMerge <- function(object = NULL, batch = "batch", ks = NULL, runningTime = FALSE, verbose = FALSE, ...){

  features <- Seurat::VariableFeatures(object)
  if(length(features) == 0){
    warning("Variable features not defined. Running 'FindVariableFeatures' function.", call. = TRUE)
    features <- Seurat::VariableFeatures(Seurat::FindVariableFeatures(object, verbose = verbose))
  }

  data <- as.matrix(Seurat::GetAssayData(object, assay = "RNA", slot = "data")[features,])
  md <- object[[]]

  if(!(batch %in% colnames(md)))
    stop(paste0(batch, "not found in object's metadata. Check the batch label."))

  if(is.null(ks)){
    nBatches <- length(unique(md[,batch]))
    ks <- rep(5, nBatches)
  }

  tmp <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = data, logcounts = data), colData = md)

  time <- system.time({
    tmp <- scMerge::scMerge(sce_combine = tmp, ctl = features, assay_name = "scMerge",
                            kmeansK = ks, batch_name = batch, plot_igraph = FALSE, verbose = FALSE, ...)
  })

  # Seurat assay
  corrData <- as.matrix(SummarizedExperiment::assay(tmp, "scMerge"))
  object[["integrated"]] <- Seurat::CreateAssayObject(counts = corrData)
  Seurat::DefaultAssay(object) <- "integrated"
  Seurat::VariableFeatures(object) <- features

  if(runningTime == FALSE)
    return(object)
  else
    return(list(object = object, time = time))
}

#' RunMNN
#'
#' @param object A seurat object to correct batch effects.
#' @param batch Batch labels.
#' @param runningTime Return the running time.
#' @param verbose Print verbose.
#' @param ... Arguments passed to other methods.
#'
#' @return Seurat object with the corrected data in the 'integrated' assay.
#' @export
RunMNN <- function(object = NULL, batch = "batch", runningTime = FALSE, verbose = FALSE, ...){

  features <- Seurat::VariableFeatures(object)
  if(length(features) == 0){
    warning("Variable features not defined. Running 'FindVariableFeatures' function.", call. = TRUE)
    features <- Seurat::VariableFeatures(Seurat::FindVariableFeatures(object, verbose = verbose))
  }

  data <- as.matrix(Seurat::GetAssayData(object, assay = "RNA", slot = "data")[features,])
  md <- object[[]]

  if(!(batch %in% colnames(md)))
    stop(paste0(batch, "not found in object's metadata. Check the batch labels."))

  time <- system.time({
    corrData <- batchelor::mnnCorrect(data, batch = md[[batch]], ...)
  })

  corrData <- SummarizedExperiment::assay(corrData, "corrected")
  object[["integrated"]] <- Seurat::CreateAssayObject(counts = corrData)
  Seurat::DefaultAssay(object) <- "integrated"
  Seurat::VariableFeatures(object) <- features

  if(runningTime == FALSE)
    return(object)
  else
    return(list(object = object, time = time))
}
