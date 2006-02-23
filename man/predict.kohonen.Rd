\name{predict.kohonen}
\alias{predict.kohonen}
\title{Classify using a trained Kohonen map}
\description{
  Map objects to a trained Kohonen map. The map may be obtained by
  supervised methods (\code{xyf}, \code{bdk}), but only the X-space is
  considered in prediction. If the classification is already known, the
  function returns the similarity for each object to the unit it is
  mapped to.
}
\usage{
predict.kohonen(object, ...)
}
\arguments{
  \item{object}{trained network.}
  \item{\dots}{optional arguments. Two arguments in particular are
    important: if a \code{data} argument is given (a data matrix with
    the same number of variables as the codebook vectors of the
    network), predictions for these objects are returned, both the units
    at which the objects are mapped and the distance to these units. If
    a classification vector is given (argument \code{classif}), the
    distances to the corresponding codebook vectors are calculated.} 
}
\value{
  Returns a list with components
  \item{unit.classif}{unit numbers to which objects in the data matrix are
    mapped.}
  \item{classif}{classes associated with the units (only for the
    supervised \code{bdk} and \code{xyf} networks).}
  \item{distances}{similarities, expressed as distances, of objects to
    codebook vectors.}
}
\author{Ron Wehrens}
\seealso{\code{\link{som}},\code{\link{xyf}},\code{\link{bdk}}}
\examples{
### Simulated example
library(MASS)

csize <- 15
c1 <- mvrnorm(csize, mu=c(5,3,4),
              Sigma=matrix(c(2, 0,0,0, 2, 1, 0, 1, 2),3,3))
c2 <- mvrnorm(csize, mu=c(5.5, 3.5, 4.5),
              Sigma=matrix(c(2, 0,0,0, 2, 1, 0, 1, 2),3,3))
c3 <- mvrnorm(csize, mu=c(0,0,0),
              Sigma=matrix(c(2, 0,0,0, 2, 1, 0, 1, 2),3,3))

X <- rbind(c1, c2, c3)
classes <- c(rep(1, csize), rep(2, csize), rep(3, csize))

bdkmod <- bdk(X, classvec2classmat(classes), somgrid(4, 4, "hexagonal"))

predict(bdkmod)
}
\keyword{classif}
