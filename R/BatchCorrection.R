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
