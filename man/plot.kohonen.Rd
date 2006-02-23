\name{plot.kohonen}
\alias{plot.kohonen}
\alias{plot.somgrid}
\alias{plot.kohchanges}
\alias{plot.kohcounts}
\alias{plot.kohcodes}
\alias{plot.kohprop}
\alias{plot.kohpred}
\alias{plot.heatkey}
\alias{plot.kohmapping}
\title{Plot kohonen object}
\description{
  Plot self-organising map, obtained from function kohonen. Several types
  of plots are supported.
}
\usage{
plot.kohonen(x, type = c("changes", "codes", "counts", "mapping",
                         "prediction", "property"),
             classif, labels=NULL, pchs=NULL, main=NULL,
             palette.name = heat.colors, ncolors, 
             zlim=NULL, property, heatkey=TRUE, contin, \dots) 
}
\arguments{
  \item{x}{kohonen object.}
  \item{type}{type of plot. (Wow!)}
  \item{classif}{classification object, as returned by
    \code{predict.kohonen}, or vector of unit numbers.
    Only needed if \code{type} equals "classifications" and
    "counts".}
  \item{labels}{labels to plot when \code{type} equals
    "classifications".}
  \item{pchs}{symbols to plot when \code{type} equals
    "classifications".}
  \item{main}{title of the plot.}
  \item{palette.name}{colors to use as unit background for "codes",
    "counts", "prediction", and "property" plotting types.}
  \item{ncolors}{number of colors to use. Default is 20 for continuous
    data, and the number of distinct values (if less than 20) for
    class data.}
  \item{zlim}{optional range for color coding of unit backgrounds.}
  \item{property}{values to use with the "property" plotting
    type.}
  \item{heatkey}{whether or not to generate a heatkey at the left side
    of the plot in the "property" and "counts" plotting types.}
  \item{contin}{whether or not the data should be seen as discrete
    (i.e. classes) or continuous in nature. Only relevant for the colour
    keys of plots of supervised networks.}
  \item{\dots}{other graphical parameters, e.g. colours of labels, or
    plotting symbols, in the "classifications" plotting type.}
}
\details{Several different types of plots are supported:
  \describe{
    \item{"changes"}{shows the changes in similarity (on a more or
      less arbitrary scale) during training. }
    \item{"codes"}{shows the unit vectors.}
    \item{"counts"}{shows the number of objects mapped to the
      individual units. Empty units are depicted in gray.}
    \item{"mapping"}{show where objects are mapped. It
      needs the "classif" argument, and a "labels" or "pchs"
      argument.}
    \item{"prediction"}{for supervised maps, the class or
      predicted value for each unit.}
    \item{"property"}{properties of each unit can be calculated and
      shown in colour code. It can be used to visualise the similarity
      of one particular object to all units in the map, to show the mean
      similarity of all units and the objects mapped to them,
      etcetera. The parameter \code{property} contains the numerical
      values. See examples below.}
  }
}

\author{Ron Wehrens}
\seealso{\code{\link{som}}, \code{\link{bdk}}, \code{\link{xyf}}}

\examples{
data(wines)
set.seed(7)

kohmap <- xyf(scale(wines), classvec2classmat(wine.classes),
              grid = somgrid(5, 5, "hexagonal"), rlen=100)
plot(kohmap, type="changes", main="plot of changes")
plot(kohmap, type="codes", main="codes plot")
plot(kohmap, type="counts", main="counts plot")
plot(kohmap, type="mapping", 
     labels=wine.classes, col=wine.classes,
     main="mapping plot")
plot(kohmap, type="mapping", 
     pchs=wine.classes, col=wine.classes,
     main="another mapping plot")

### Made-up labels :-)
plot(kohmap, type="prediction",
     labels=c("Chianti", "Barolo", "Montepulciano"),
     palette.name = rainbow,
     main="unit class prediction", cex=.8)

### Plot mean similarity of mapped objects to their unit codebook vector
mapping <- predict(kohmap, classif=kohmap$classif)
hits <- sort(unique(kohmap$classif))
similarities <- rep(NA, 25)
for (i in seq(along=hits))
  similarities[hits[i]] <- mean(mapping$distances[kohmap$classif == hits[i]])

plot(kohmap, type="property", property=similarities, main="property plot")
}
\keyword{classif}
