\name{supersom}
\alias{supersom}
\title{Super-organising maps}
\description{
  An extension of \code{xyf} maps to multiple data layers, possibly with
  different numbers of variables (though equal numbers of objects). NAs
  are allowed (see below). A weighted distance over all layers is
  calculated to determine the winning units during training.
}

\usage{
supersom(data, grid=somgrid(), rlen = 100, alpha = c(0.05, 0.01),
    radius = quantile(nhbrdist, 0.67), 
    toroidal = FALSE, whatmap = NULL, weights = 1,
    maxNA.fraction = .5, keep.data = TRUE)
}

\arguments{
  \item{data}{list of data matrices.}
  \item{grid}{a grid for the representatives: see \code{\link{somgrid}}.}
  \item{rlen}{the number of times the complete data set will be
    presented to the network.}
  \item{alpha}{learning rate, a vector of two numbers indicating the
    amount of change. Default is to decline linearly from 0.05 to 0.01
    over \code{rlen} updates.}
  \item{radius}{the starting radius of the neighbourhood to be used for
    each update: the decrease is linear over \code{rlen} updates and
    reaches the value of 1 after one-third of the iterations.}
  \item{toroidal}{if TRUE, the edges of the map are joined. Note
    that in a hexagonal toroidal map, the number of rows must be even.}
  \item{whatmap}{For \code{supersom} maps: what layers to use in the
    mapping.}
  \item{weights}{the weights given to individual layers. Default is 1/n,
    with n the number of layers.}
  \item{maxNA.fraction}{the maximal fraction of values that may be NA to
    prevent the row or column to be removed.}
  \item{keep.data}{save data in return value.}
}

\value{
  an object of class "kohonen" with components
  \item{data}{data matrix, only returned if \code{keep.data == TRUE}.}
  \item{na.rows}{indices of objects (rows) that are removed because at
    least one of the layers has to many NAs for these objects.}
  \item{unit.classif}{winning units for all data objects,
    only returned if \code{keep.data == TRUE}.}
  \item{distances}{distances of objects to their corresponding winning
    unit, only returned if \code{keep.data == TRUE}.}
  \item{grid}{the grid, an object of class \code{\link{somgrid}}.}
  \item{codes}{a list of matrices containing codebook vectors.}
  \item{changes}{matrix of mean average deviations from code vectors;
    every map corresponds with one column.}
  \item{toroidal}{whether a toroidal map is used.}
  \item{weights}{For \code{supersom} maps: weights of layers uses in the
    mapping.}
  \item{whatmap}{For \code{supersom} maps: what layers to use in the
    mapping.}
  \item{method}{type of map, here "supersom".}
}

\references{W.J. Melssen, R. Wehrens, and L.M.C. Buydens.
Chemom. Intell. Lab. Syst., 83, 99-113 (2006).}
\author{Ron Wehrens}

\seealso{\code{\link{somgrid}}, \code{\link{plot.kohonen}}}
\examples{
}

\keyword{classif}

