
#' myPCAplot
#'
#' Function to plot PCA
#'
#' @param df Data frame containing the PCA dimensions to plot and the labels to use as colors. PCA data columns should be named as "Dim1" and "Dim2"
#' @param group_by Labels to use as colors.
#' @param text_size Plot text size
#' @param point_size Plot point size
#' @param alpha Point transparency
#' @param legend_position The position of the legend. It can take values as "top", bottom", "left", "right", or "none.
#' @param legend_point_size Legend point size.
#' @param ... Arguments passed to other methods.
#'
#' @return A ggplot object containing the PCA plot.
#' @export
myPCAplot <- function(df = NULL, group_by = NULL, text_size = 20,
                      point_size = 2, alpha = 1.0, legend_position = "right",
                      legend_point_size = NULL, ...){


  p <- ggplot(df, aes_string(x = "Dim1", y = "Dim2", color = group_by)) +
    geom_point(size = point_size, alpha = alpha) +
    theme_classic() +
    theme(text = element_text(size = text_size), legend.position = legend_position, ...)

  if(!is.null(legend_point_size)){
    p <- p + guides(colour = guide_legend(override.aes = list(size=legend_point_size)))
  }

  return(p)
}
