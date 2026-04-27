<!-- badges: start -->
[![R-CMD-check](https://github.com/Daniel-1232/NicheDiv/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Daniel-1232/NicheDiv/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

# NicheDiv R package

Performs DAPC niche divergence tests

NicheDiv provides functions for environmental data extraction, occurrence thinning,
background filtering, transformation of skewed variables, trimming to analogous
environmental space, DAPC-based niche divergence analyses, and plotting of results.

## Installation
Install R package from GitHub by running:

```r
install.packages("remotes")
remotes::install_github("Daniel-1232/NicheDiv")
```

## Main functionality
- convert integer columns to numeric
- crop background points to buffered occurrence extents
- down-sample occurrence or background data
- thin occurrence records spatially
- transform skewed environmental variables
- remove low-information variables
- filter to analogous environmental variables
- trim to analogous environmental space
- run DAPC-based niche divergence analyses
- plot DAPC results and predictor contributions

## Minimal example
```r
library(NicheDiv)

# Example toy data
x <- data.frame(
  Longitude = c(-100, -99.5, -99),
  Latitude = c(40, 40.5, 41),
  bio1 = c(10L, 11L, 12L),
  bio12 = c(500, 520, 510)
)

# Convert integer columns to numeric
x2 <- convert.integer.to.numeric(x)

# Thin occurrence records
x_thin <- thin.occurrence(
  occurrence.data = x2,
  longitude.col = "Longitude",
  latitude.col = "Latitude",
  thinning.dist.km = 1,
  calc.Morans.I = FALSE,
  plot.Morans.I = FALSE,
  verbose = FALSE
)

x_thin
```

## Author
Daniel Schönberger
