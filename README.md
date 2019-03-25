# PlotSeurat

The goal of PlotSeurat is to enable custom plots derived from what Seurat has 
to offer.

## Installation

You can install the  PlotSeurat from this repository with:

``` r
devtools::install_github(repo = "chenyenchung/PlotSeurat")
```

## Example

This is a basic example of plotting a featureplot with multiple features 
overlaid with red, green, and blue respectively:

``` r
MultiFeaturePlot(object, features = c("feature1", "feature2", "feature3"), 
                 title = "My Awesome Plot")
```

