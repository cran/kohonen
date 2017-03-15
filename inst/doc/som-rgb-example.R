## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ------------------------------------------------------------------------
library(kohonen)
library(grid)
cells <- 20
colors <- rbind(c(1,0,0), c(1,1,0), c(0,1,0), c(0,1,1), c(0,0,1), c(1,0,1))
init <- replicate(3, runif(cells * cells))
somgrid <- somgrid(cells, cells, "rectangular", neighbourhood.fct = "gaussian")
som.rgb <- supersom(
  data = colors,
  grid = somgrid,
  init = init,
  rlen = 1000,
  mode = "online",
  keep.data = FALSE)


## ---- echo=FALSE---------------------------------------------------------
  plot(som.rgb)
  colors <- (som.rgb$codes[[1]])
  map <- rgb(colors[,1], colors[,2], colors[,3])
  col <- matrix(map, nrow = som.rgb$grid$xdim, ncol = som.rgb$grid$ydim)
  grid.raster(col, interpolate=FALSE)

## ------------------------------------------------------------------------
cells <- 20
colors <- rbind(c(1,0,0), c(1,1,0), c(0,1,0), c(0,1,1), c(0,0,1), c(1,0,1))
init <- replicate(3, runif(cells * cells))
somgrid <- somgrid(cells, cells, "rectangular", neighbourhood.fct = "bubble")
som.rgb <- supersom(
  data = colors,
  grid = somgrid,
  init = init,
  rlen = 1000,
  mode="online",
  keep.data = FALSE)


## ---- echo=FALSE---------------------------------------------------------
  plot(som.rgb)
  colors <- (som.rgb$codes[[1]])
  map <- rgb(colors[,1], colors[,2], colors[,3])
  col <- matrix(map, nrow = som.rgb$grid$xdim, ncol = som.rgb$grid$ydim)
  grid.raster(col, interpolate=FALSE)

