
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


#' myFeaturePlot
#'
#' Function to plot a feature
#'
#' @param df Data frame containing the PCA dimensions to plot and the feature to use as colors. PCA data columns should be named as "Dim1" and "Dim2"
#' @param feature Feature to use as color scale.
#' @param text_size Plot text size
#' @param point_size Plot point size
#' @param alpha Point transparency
#' @param low_color Scale low color
#' @param high_color Scale high color
#' @param legend_position The position of the legend. It can take values as "top", bottom", "left", "right", or "none.
#' @param legend_point_size Legend point size.
#' @param order Whether the points should be ordered.
#' @param ... Arguments passed to other methods.
#'
#' @return A ggplot object containing the feature plot.
#' @export
myFeaturePlot <- function(df = NULL, feature = NULL,
                          text_size = 20, point_size = 2,
                          alpha = 1.0, low_color = "gray",
                          high_color = "tomato", legend_position = "right",
                          legend_point_size = NULL, order = TRUE,
                          discrete = FALSE, palette = NULL,
                          ...){

  # Order the feature
  if(order){
    df <- df[order(df[feature]),]
  }

  # Colors palette
  if(is.null(palette)){
    if(length(unique(df$feature)) == 2){
      palette <- c(low_color, high_color)
    }else{
      palette <- 1:length(unique(df$feature))
    }
  }

  # Create the initial plot
  p <-  ggplot(df, aes_string(x = "Dim1", y = "Dim2", color = feature)) +
    geom_point(size = point_size, alpha = alpha) +
    theme_classic() +
    theme(text = element_text(size = text_size), legend.position = legend_position, ... ) +
    ggtitle("Feature plot")

  # Change the legend's point size
  if(!is.null(legend_point_size)){
    p <- p + guides(colour = guide_legend(override.aes = list(size=legend_point_size)))
  }

  # Change the points' color
  if(discrete){
    p <- p + scale_color_manual(values = palette)
  }else{
    p <- p + scale_color_gradient(low = low_color, high =  high_color)
  }

  return(p)
}

#' PlotGroups
#'
#' @param object Seurat object with calculated UMAP
#' @param group Labels to use as colors.
#' @param combine Whether the plots should be combined
#' @param ncol If plots are combined, the number of columns of the final plot.
#' @param text_size Size of text.
#' @param point_size Size of points.
#' @param alpha Transparency of points.
#' @param low_color Color for background points.
#' @param high_color Color for main points.
#' @param legend_position Positions of plot legends.
#' @param legend_point_size Size of legends' points.
#' @param order Whether the points should be ordered.
#' @param palette Color pallete to use.
#' @param ... Arguments passed to other methods.
#'
#' @return
#' @export
PlotGroups <- function(object = NULL, group = NULL,
                       combine = FALSE, ncol = 2,
                       text_size = 20, point_size = 2,
                       alpha = 1.0, low_color = "gray",
                       high_color = "tomato", legend_position = "right",
                       legend_point_size = NULL, order = TRUE,
                       palette = NULL, ...){

  if(is.null(group)){
    print("The input parameter 'group' needs to be specified.")
    stop(call. = TRUE)
  }

  if(is.null(object)){
    print("The input parameter 'object' needs to be specified.")
    stop(call. = TRUE)
  }

  md <- object[[]]
  umapData <- Embeddings(object = object, reduction = "umap")

  plots <- list()

  for(g in unique(md[,group])){

    feature <- rep(0, nrow(md))
    idx <- which(md[,group] == g)
    feature[idx] <- 1
    feature <- factor(feature)

    df <- data.frame("Dim1" = umapData[,1], "Dim2" = umapData[,2], feature = feature)

    p <- RNAseqAnalysis::myFeaturePlot(df = df, feature = "feature",
                                       discrete = TRUE, text_size = text_size, point_size = point_size,
                                       alpha = alpha, low_color = low_color,
                                       high_color = high_color, legend_position = legend_position,
                                       legend_point_size = legend_point_size, order = order,
                                       palette = palette, ...) +
      ggtitle(g)

    plots[[g]] <- p
  }

  if(combine){
    p <- plots[[1]]
    for(i in 2:length(plots)){
      p <- p + plots[[i]]
    }
    p <- p + plot_layout(ncol = ncol)
  }else{
    p <- plots
  }

  return(p)
}

#' TextSize
#'
#' @param size Text size
#'
#' @return A theme with the selected text size.
#' @export
TextSize <- function(size = 10 ){
  tmp_theme <- theme(text = element_text(size = size))
  return(tmp_theme)
}

#' LegendPosition
#'
#' @param position Legend position
#'
#' @return Theme with the legened in the selected position. Available positions are :"none", "right", "left", "top", and "bottom".
#' @export
LegendPosition <- function(position = "right"){
  tmp_theme <- theme(legend.position = position)
}

