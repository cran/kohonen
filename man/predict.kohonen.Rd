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
predict.kohonen(object, newdata, trainX, trainY, unit.predictions,
                 threshold = 0, ...)
}
\arguments{
  \item{object}{Trained network.}
  \item{newdata}{Data matrix for which predictions are to be made. If
    not given, defaults to the training data (when available).}
  \item{trainX}{Training data for obtaining predictions for unsupervised
    maps; necessary for \code{som} maps trained with the keep.data =
    FALSE option.}
  \item{trainY}{Values for the dependent variable for the training data;
    only necessary for \code{som} maps.}
  \item{unit.predictions}{Possible override of the predictions for each
    unit.}
  \item{threshold}{Used in class predictions; see
    \code{\link{classmat2classvec}}.}
  \item{\dots}{Currently not used.} 
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
\seealso{\code{\link{som}},\code{\link{xyf}},\code{\link{bdk}},
  \code{\link{map}}}
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
