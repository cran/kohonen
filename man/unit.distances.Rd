\name{unit.distances}
\alias{unit.distances}
\title{Calculate distances between units in a SOM}
\description{
Calculate distances between units in a SOM. 
}
\usage{
unit.distances(grid, toroidal)
}
\arguments{
  \item{grid}{an object of class \code{somgrid}.}
  \item{toroidal}{if true, edges of the map are joined so that the
    topology is that of a torus.}
}
\value{
  Returns a matrix containing distances.
}
\author{Ron Wehrens}
\examples{
library(kohonen)
data(wines)

kohmap <- som(wines, grid = somgrid(5, 5, "hexagonal"), rlen=100)

par(mfrow=c(1,2))
dists <- unit.distances(kohmap$grid, toroidal=FALSE)
plot(kohmap, type="property", property=dists[1,],
     main="Distances to unit 1", zlim=c(0,4))
dists <- unit.distances(kohmap$grid, toroidal=TRUE)
plot(kohmap, type="property", property=dists[1,],
     main="Distances to unit 1 (toroidal)", zlim=c(0,4))
}
\keyword{classif}
