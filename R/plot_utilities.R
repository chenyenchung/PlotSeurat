#' Bin Expression to 256 Bins
#'
#' \code{bin_expression()} takes a data frame and bin every column to 256 bins,
#' and generate Hex-coded color code according to the composite of the columns.
#' It takes three or less columns.
#' @param featuresdf a data frame containing expression information
#'
#' @return a character vector containing hexcoded color with a length the same
#' as the row number of the input data frame
bin_expression <- function (featuresdf) {
  if (class(featuresdf) != "data.frame") {
    stop(paste("bin_expression only takes data frames as input.",
               "Check the type of data FetchData() returns."))
  }
  binls <- lapply(as.list(featuresdf), function(x) {
    unit <- max(x) / 255
    binexp <- round(x / unit)
    hexbinexp <- format(as.hexmode(binexp), width = 2)
    return(hexbinexp)
  })
  if (length(binls) == 2) {
    binls[[3]] <- rep("00", length(binls[[1]]))
  }
  if (length(binls) == 1) {
    binls[[2]] <- rep("00", length(binls[[1]]))
    binls[[3]] <- rep("00", length(binls[[1]]))
  }
  hexcolor <- paste0("#", binls[[1]], binls[[2]], binls[[3]])
  return(hexcolor)
}

#' Extract Figure Legend from Individual Plots for Later Integration to the
#' Main Plot
#'
#' Since \code{ggplot2} does not support superimposing multiple
#' \code{geom_point()} with different color codes,
#' \code{\link{MultiFeaturePlot}()} bypasses the limitation by generating one
#' custom color scale integrating multiple expression data. To get legends for
#' that, individual \code{geom_point()} ggplots are generated here, and
#' extracted for later integration with the main plot.
#' @param featurekey a character vector containing the features prefixed by
#' keys (if using Seuratv3) to be plotted
#' @param plotdf a data frame containing cell embedings, expression level, and
#' color codes genearted by \code{\link{bin_expression}()}
#'
#' @return a ggplot2 object
make_legend <- function (featurekey, plotdf) {
  dim1 <- ggplot2::sym(colnames(plotdf)[2])
  dim2 <- ggplot2::sym(colnames(plotdf)[3])
  num_key <- 1
  legends <- list()
  for (key in featurekey) {
    # Quoting column names to deal with ggplot2 NSE
    gene <- ggplot2::sym(key)

    # The first gene will use red gradient, then green, and then blue
    col <- c("00", "00", "00")
    col[num_key] <- "FF"
    highcol <- paste0("#", paste(col, collapse = ""))

    legendonly <- cowplot::get_legend(
      ggplot2::ggplot(plotdf, ggplot2::aes(x = !! dim1, y = !! dim2,
                                           color = !! gene)) +
        ggplot2::geom_point() +
        ggplot2::scale_color_gradient(low = "#000000", high = highcol)
    )
    legends[[num_key]] <- legendonly
    num_key <- num_key + 1
  }

  if (length(legends) == 1) {
    return(cowplot::plot_grid(legends[[1]], ncol = 1))
  } else if (length(legends) == 2) {
    return(cowplot::plot_grid(legends[[1]], legends[[2]], ncol = 1))
  }
  return(cowplot::plot_grid(legends[[1]], legends[[2]], legends[[3]], ncol = 1))

}