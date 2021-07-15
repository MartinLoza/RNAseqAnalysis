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

#' RunScanorama
#'
#' @param object A seurat object to correct batch effects.
#' @param batch Batch labels.
#' @param runningTime Return the running time.
#' @param verbose Print verbose.
#' @param ... Arguments passed to other methods.
#'
#' @return Seurat object with the corrected data in the 'integrated' assay.
#' @export
RunScanorama <- function(object = NULL, batch = "batch", runningTime = FALSE, verbose = FALSE, ...){

  Scanorama <- reticulate::import("scanorama")
  datal <- list()
  genel <- list()
  features <- Seurat::VariableFeatures(object)

  if(length(features) == 0){
    warning("Variable features not defined. Running 'FindVariableFeatures' function.", call. = TRUE)
    features <- Seurat::VariableFeatures(Seurat::FindVariableFeatures(object, verbose = verbose))
  }

  if(!(batch %in% colnames(object[[]])))
    stop(paste0(batch, "not found in object's metadata. Check the batch labels."))

  objectl <- Seurat::SplitObject(object, split.by = batch)

  for(i in seq_len(length(objectl))){
    datal[[i]] <- Seurat::GetAssayData(objectl[[i]], assay = "RNA", slot = "data")[features,] # Normalized counts
    datal[[i]] <- as.matrix(datal[[i]])
    datal[[i]] <- t(datal[[i]]) # Cell x genes

    genel[[i]] <- features
  }

  time <- system.time({
    corrDatal <- Scanorama$correct(datasets_full = datal, genes_list = genel, return_dense = TRUE)
  })

  corrData <- Reduce(rbind, corrDatal[[1]])
  corrData <- t(corrData)
  rownames(corrData) <- corrDatal[[2]]
  colnames(corrData) <- unlist(sapply(objectl,colnames))

  # Same cell names as the original object
  corrData <- corrData[,colnames(object)]

  ## Create Seurat assay
  object[["integrated"]] <- Seurat::CreateAssayObject(counts = corrData)
  Seurat::DefaultAssay(object) <- "integrated"
  Seurat::VariableFeatures(object) <- features

  if(runningTime == FALSE)
    return(object)
  else
    return(list(object = object, time = time))
}

#' RunLiger
#'
#' @param object A seurat object to correct batch effects.
#' @param batch Batch labels.
#' @param k Inner dimension of factorization (number of factors)
#' @param runningTime Return the running time.
#' @param verbose Print verbose.
#' @param ... Arguments passed to other methods.
#'
#' @return Seurat object with the corrected data in the 'Liger' reduction.
#' @export
RunLiger <- function(object = NULL, batch = "batch", k = 30, runningTime = FALSE, verbose = FALSE, ...){

  features <- Seurat::VariableFeatures(object)

  if(length(features) == 0){
    warning("Variable features not defined. Running 'FindVariableFeatures' function.", call. = TRUE)
    features <- Seurat::VariableFeatures(Seurat::FindVariableFeatures(object, verbose = verbose))
  }

  if(!(batch %in% colnames(object[[]])))
    stop(paste0(batch, "not found in object's metadata. Check the batch labels."))

  tmp <- object[features,]

  time <- system.time({
    tmp <- Seurat::ScaleData(tmp, split.by = "batch", do.center = FALSE, verbose = verbose, ...)
    tmp <- SeuratWrappers::RunOptimizeALS(tmp, k = k, split.by = "batch", ...)
    tmp <- SeuratWrappers::RunQuantileNorm(tmp, split.by = "batch", ...)
  })

  object[["Liger"]] <- tmp[["iNMF"]]

  if(runningTime == FALSE)
    return(object)
  else
    return(list(object = object, time = time))
}

#' RunComBat
#'
#' @param object A seurat object to correct batch effects.
#' @param batch Batch labels.
#' @param runningTime Return the running time.
#' @param verbose Print verbose.
#' @param ... Arguments passed to other methods.
#'
#' @return Seurat object with the corrected data in the 'integrated' assay.
#' @export
RunComBat <- function(object = NULL, batch = "batch", runningTime = FALSE, verbose = FALSE, ...){

  features <- Seurat::VariableFeatures(object)
  if(length(features) == 0){
    warning("Variable features not defined. Running 'FindVariableFeatures' function.", call. = TRUE)
    features <- Seurat::VariableFeatures(Seurat::FindVariableFeatures(object, verbose = verbose))
  }

  data <- as.matrix(Seurat::GetAssayData(object, assay = "RNA", slot = "data")[features,])
  md <- object[[]]

  if(!(batch %in% colnames(md)))
    stop(paste0(batch, "not found in object's metadata. Check the batch label."))

  time <- system.time({
    corrData <- sva::ComBat(dat = data, batch = md[[batch]], ...)
  })

  object[["integrated"]] <- Seurat::CreateAssayObject(counts = corrData)
  Seurat::DefaultAssay(object) <- "integrated"
  Seurat::VariableFeatures(object) <- features

  if(runningTime == FALSE)
    return(object)
  else
    return(list(object = object, time = time))
}

#' RunHarmony
#'
#' @param object A seurat object to correct batch effects.
#' @param batch Batch labels.
#' @param dims Dimensions to use in the correction.
#' @param runningTime Return the running time.
#' @param verbose Print verbose.
#' @param ... Arguments passed to other methods.
#'
#' @return Seurat object with the corrected data in the 'harmony' reduction.
#' @export
RunHarmony <- function(object = NULL, batch = "batch", dims = 10, runningTime = FALSE, verbose = FALSE, ...){

  if(!("pca" %in% Seurat::Reductions(object))){
    if(verbose)
      print("Running PCA.")
    time <- system.time({
      object <- GetPCA(object = object, dims = dims, verbose = verbose, ...)
      object <- harmony::RunHarmony(object = object, group.by.vars = batch, dims.use = 1:dims, verbose = verbose, ...)
    })
  }else{
    time <- system.time({
      object <- harmony::RunHarmony(object = object, group.by.vars = batch, dims.use = 1:dims, verbose = verbose, ...)
    })
  }

  if(runningTime == FALSE)
    return(object)
  else
    return(list(object = object, time = time))
}

#' Title
#'
#' @param object A seurat object to correct batch effects.
#' @param batch Batch labels.
#' @param reduction Reduction to use.
#' @param dims Number of dimensions to use.
#' @param per Percentages of the mean batch size.
#' @param acceptance Return the acceptance rate.
#' @param verbose Print verbose.
#'
#' @return kBET mean score.
#' @export
RunKBET <- function(object = NULL, batch = "batch", reduction = "pca", dims = 10, per = 0.1, acceptance = TRUE, verbose = FALSE){

  md <- object[[]]
  if(!(reduction %in% Seurat::Reductions(object)))
    stop(paste0(reduction, " not found in the object's reductions."))

  if(!(batch %in% colnames(md)))
    stop(paste0(batch, " not found in the object's meta data."))

  data <- as.data.frame(Seurat::Embeddings(object = object, reduction = reduction)[,1:dims])
  meanBatch <- mean(table(md[[batch]]))

  scores <- lapply(per, function(p){
    k0 = floor(p*(meanBatch))
    score <- mean(kBET::kBET(df = data, batch = md[[batch]], do.pca = FALSE,
                             heuristic = FALSE, k0 = k0,
                             plot = FALSE)$stats$kBET.observed)
    return(score)
  })

  scores <- unlist(scores)
  scores <- mean(scores)

  if(acceptance)
    scores <- 1-scores

  return(scores)
}

#' RunSilhouette
#'
#' @param object A seurat object to correct batch effects.
#' @param batch Batch labels.
#' @param reduction Reduction to use.
#' @param dims Number of dimensions to use.
#'
#' @return Silhouette width score.
#' @export
RunSilhouette <- function(object = NULL, batch = "celltype", reduction = "pca", dims = 10){

  md <- object[[]]

  if(!(reduction %in% Seurat::Reductions(object)))
    stop(paste0(reduction, " not found in the object's reductions."))

  if(!(batch %in% colnames(md)))
    stop(paste0(batch, " not found in the meta data."))

  batch <- factor(md[[batch]])

  pcaData <- as.matrix(Seurat::Embeddings(object = object, reduction = reduction)[,1:dims])
  pcaData <- list(x = pcaData)

  score <- kBET::batch_sil(pca.data = pcaData, batch = batch, nPCs = dims)

  return(score)
}
