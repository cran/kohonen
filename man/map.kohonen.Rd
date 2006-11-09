\name{map.kohonen}
\alias{map.kohonen}
\alias{map}
\title{Map data to a supervised or unsupervised SOM}
\description{
  Map a data matrix onto a trained SOM.
}
\usage{
map.kohonen(x, newdata, ...)
}
\arguments{
  \item{x}{A trained supervised or unsupervised SOM obtained from
    functions \code{som}, \code{xyf}, or \code{bdk}.}
  \item{newdata}{Data matrix, with rows corresponding to objects.}
  \item{\dots}{Currently ignored.}
}
\value{
  A list with elements
  \item{unit.classif}{a vector of units that are closest to the objects
    in the data matrix.}
  \item{dists}{distances (currently only Euclidean distances) of the
    objects to the units.}
}
\author{Ron Wehrens}
\seealso{\code{\link{predict.kohonen}}}
\examples{
data(wines)
set.seed(7)

training <- sample(length(wine.classes), 120)
Xtraining <- scale(wines[training,])
somnet <- som(Xtraining, somgrid(5, 5, "hexagonal"))

mapping <- map(somnet,
               scale(wines[-training,],
                     center=attr(Xtraining, "scaled:center"),
                     scale=attr(Xtraining, "scaled:scale")))
}
\keyword{classif}
