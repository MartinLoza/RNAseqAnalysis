#' GetTopGenes
#'
#' @param deg Data frame with DGE results.
#' @param pval P-value threshold.
#' @param fold Fold chagne threshold.
#' @param type Type of selected fold change. Available options are: "double", negatives and positive, "positive", and "negative".
#' @param ntop Number of top fold changes.
#' @param group If top. Group label to divide the data frame.
#'
#' @return Data frame with the top genes.
#' @export
GetTopGenes <- function(deg = NULL, pval = 0.05, fold = 1.5, type = "double", ntop = NULL, group = NULL){

  FCthreshold <- log2(fold)

  if(type ==  "double"){
    topGenes <- deg %>% filter(p_val_adj < pval & abs(avg_log2FC) > FCthreshold)
  }else if(type == "positive"){
    topGenes <- deg %>% filter(p_val_adj < pval & avg_log2FC > FCthreshold)
  }else{
    topGenes <- deg %>% filter(p_val_adj < pval & avg_log2FC < FCthreshold)
  }

  if(!is.null(ntop)){
    if(!is.null(group)){
      warning("Top genes by clusters", call. = TRUE)
      group = "cluster"
    }

    topGenes <- topGenes %>% group_by(cluster) %>% top_n(n = ntop)
  }

  return(topGenes)
}
