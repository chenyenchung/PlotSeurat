#' Plotting Multiple Features on a Dimensionality Reduced Space
#'
#' \code{MultiFeaturePlot()} plots cells in a Seurat object into a t-SNE space
#' or other dimensionality reduced space for visualization. Each point will be
#' colored according to log-normalized expression level (from \code{@data} slot
#' of a Seurat object). The first gene will be colored in red, the second gene
#' will be colored in green, and the third gene will be colored in blue.
#'
#' Specifically, the gene expression is translated into color by first taking
#' the full range of expression from zero to highest, binning them into 256
#' bins, and use the binned expression level as 8-bit color code.
#'
#' @param object a Seurat object (v2 or v3)
#' @param features a character vector of features to plot
#' @param assay a character string indicating the assay to plot (SeuratV3 only)
#' @param pt.size a number indicating the size of points on the plot
#' @param reduction a character string indicating the reduction to embed the
#' cells on the plot (e.g., "tsne" or "pca")
#' @param cells.use a character vector of cell names indicating the cells to
#' include in the plot (try \code{\link[Seurat]{WhichCells}()})
#' @param title a character string for customized plot title. If not set, the
#' default is the feature names
#'
#' @return a ggplot object
#' @export
MultiFeaturePlot <- function (object, features, assay = "RNA", pt.size = 1,
                                reduction = "tsne", cells.use = NULL,
                                title = NULL) {
  # Set default arguments if not provided
  if (is.null(cells.use)) {
    cells.use <- row.names(methods::slot(object, name = "meta.data"))
  }

  if (is.null(title)) {
    title <- paste(features, collapse = " / ")
  }

  # Check feature upper limit
  if (length(features) > 3) {
    stop("At most 3 features can be plotted at the same.")
  }

  # Extract cell embedings in dimensionality reduction
  reduction <- tolower(reduction)
  if (utils::packageVersion("Seurat")[[1, 1]] == 2) {
    reductionobj <- methods::slot(object, name = "dr")[[reduction]]
  } else {
    reductionobj <- methods::slot(object, name = "reductions")[[reduction]]
  }
  cellemb <- data.frame(methods::slot(reductionobj, name = "cell.embeddings"))

  if (utils::packageVersion("Seurat")[[1, 1]] == 3) {
    # Get keys for the assay
    assayobj <- methods::slot(object, name = "assays")[[assay]]
    assaykey <- methods::slot(assayobj, name = "key")
  }



  # Extract expression matrix
  if (utils::packageVersion("Seurat")[[1,1]] == 2) {
    featuresexp <- Seurat::FetchData(object, vars.all = features,
                             cells = cells.use)
    featuresexp <- data.frame(featuresexp)
    featurekey <- features
  } else {
    featurekey <- paste0(assaykey, features)
    featuresexp <- Seurat::FetchData(object, vars = featurekey,
                             cells = cells.use)
  }

  # Binning the expression levels from 0 - 255
  featuresexp$color <- bin_expression(featuresdf = featuresexp)

  # Create plotting df
  plotdf <- merge(cellemb, featuresexp, by = "row.names")

  # Plotting
  dim1 <- ggplot2::sym(colnames(plotdf)[2])
  dim2 <- ggplot2::sym(colnames(plotdf)[3])

  base_plot <- ggplot2::ggplot(plotdf, ggplot2::aes(x = !!dim1, y = !!dim2)) +
    ggplot2::geom_point(color = plotdf$color, size = pt.size) +
    ggplot2::ggtitle(label = title)
  # Make legends for each gene
  legend <- make_legend(featurekey = featurekey, plotdf = plotdf)

  result <- cowplot::plot_grid(base_plot, legend,
                      rel_widths = c(5, 1), ncol = 2)
  print(result)
  return(result)
}