\name{plot.kohonen}
\alias{plot.kohonen}
\alias{plot.somgrid}
\alias{plot.kohchanges}
\alias{plot.kohcounts}
\alias{plot.kohcodes}
\alias{plot.kohprop}
\alias{plot.kohquality}
\alias{plot.heatkey}
\alias{plot.kohmapping}
\alias{plot.kohUmatrix}
\title{Plot kohonen object}
\description{
  Plot self-organising map, obtained from function kohonen. Several types
  of plots are supported.
}
\usage{
\method{plot}{kohonen}(x, type = c("codes", "changes", "counts",
                 "dist.neighbours", "mapping", "property", "quality"),
     classif = NULL, labels = NULL, pchs = NULL, main = NULL,
     palette.name = heat.colors, ncolors, bgcol = NULL,
     zlim = NULL, heatkey = TRUE, property, contin,
     whatmap = NULL, codeRendering = NULL, keepMargins = FALSE,
     heatkeywidth = .2, \dots) 
}
\arguments{
  \item{x}{kohonen object.}
  \item{type}{type of plot. (Wow!)}
  \item{classif}{classification object, as returned by
    \code{predict.kohonen}, or vector of unit numbers.
    Only needed if \code{type} equals "mapping" and
    "counts".}
  \item{labels}{labels to plot when \code{type} equals
    "mapping".}
  \item{pchs}{symbols to plot when \code{type} equals
    "mapping".}
  \item{main}{title of the plot.}
  \item{palette.name}{colors to use as unit background for "codes",
    "counts", "prediction", "property", and "quality" plotting types.}
  \item{ncolors}{number of colors to use for the unit
    backgrounds. Default is 20 for continuous data, and the number of
    distinct values (if less than 20) for categorical data.}
  \item{bgcol}{optional argument to colour the unit backgrounds for the
    "mapping" and "codes" plotting type. Defaults to "gray" and
    "transparent" in both types, respectively.}
  \item{zlim}{optional range for color coding of unit backgrounds.}
  \item{heatkey}{whether or not to generate a heatkey at the left side
    of the plot in the "property" and "counts" plotting types.}
  \item{property}{values to use with the "property" plotting
    type.}
  \item{contin}{whether or not the data should be seen as discrete
    (i.e. classes) or continuous in nature. Only relevant for the colour
    keys of plots of supervised networks.}
  \item{whatmap}{For \code{supersom} maps and a "codes" plot: what maps
    to show.}
  \item{codeRendering}{How to show the codes. Possible choices:
    "segments", "stars" and "lines".}
  \item{keepMargins}{if \code{FALSE} (the default), restore the original
    graphical parameters after plotting the kohonen map. If \code{TRUE},
    one retains the map coordinate system so that one can add symbols to
    the plot, or map unit numbers using the \code{identify} function.}
  \item{heatkeywidth}{width of the colour key; the default of 0.2 should
    work in most cases but in some cases, e.g. when plotting multiple
    figures, it may need to be adjusted.}
  \item{\dots}{other graphical parameters, e.g. colours of labels, or
    plotting symbols, in the "mapping" plotting type.}
}
\details{Several different types of plots are supported:
  \describe{
    \item{"changes"}{shows the mean distance to the closest codebook vector
      during training.} 
    \item{"codes"}{shows the codebook vectors.}
    \item{"counts"}{shows the number of objects mapped to the
      individual units. Empty units are depicted in gray.}
    \item{"dist.neighbours"}{shows the sum of the distances to all
      immediate neighbours. This kind of visualisation is also known as
      a U-matrix plot. Units near a class boundary can be expected to
      have higher average distances to their neighbours. Only available
      for the "som" and "supersom" maps, for the moment.}
    \item{"mapping"}{shows where objects are mapped. It
      needs the "classif" argument, and a "labels" or "pchs"
      argument.}
    \item{"property"}{properties of each unit can be calculated and
      shown in colour code. It can be used to visualise the similarity
      of one particular object to all units in the map, to show the mean
      similarity of all units and the objects mapped to them,
      etcetera. The parameter \code{property} contains the numerical
      values. See examples below.}
    \item{"quality"}{shows the mean distance of objects mapped to a
      unit to the codebook vector of that unit. The smaller the
      distances, the better the objects are represented by the codebook
      vectors.}
  }
}

\author{Ron Wehrens}
\seealso{\code{\link{som}}, \code{\link{bdk}}, \code{\link{xyf}}}

\examples{
data(wines)
set.seed(7)

kohmap <- xyf(scale(wines), classvec2classmat(wine.classes),
              grid = somgrid(5, 5, "hexagonal"), rlen=100)
plot(kohmap, type="changes")
plot(kohmap, type="codes", main = c("Codes X", "Codes Y"))
plot(kohmap, type="counts")
plot(kohmap, type="quality")
plot(kohmap, type="mapping", 
     labels=wine.classes, col=wine.classes+1,
     main="mapping plot")
xyfpredictions <- classmat2classvec(predict(kohmap)$unit.predictions)
bgcols <- c("gray", "pink", "lightgreen")
plot(kohmap, type="mapping", col=wine.classes+1,
     pchs=wine.classes, bgcol=bgcols[as.integer(xyfpredictions)], 
     main="another mapping plot")

}
\keyword{classif}
