\name{xyf}
\alias{xyf}
\title{Supervised version of Kohonen's self-organising maps}
\description{
  Supervised version of self-organising maps for mapping
  high-dimensional spectra or patterns to 2D. The name stands for X-Y
  fused SOMs. One vector for each object is created by concatenating X
  and Y, and a SOM is trained in the usual way, with one exception: the
  distance of an object to a unit is the sum of separate distances for X
  and Y spaces. Prediction is done only using the X-space.
}

\usage{
xyf(data, Y, grid=somgrid(), rlen = 100, alpha = c(0.05, 0.01),
    radius = quantile(nhbrdist, 0.67), 
    xweight = 0.5, toroidal = FALSE, keep.data = TRUE)
}

\arguments{
  \item{data}{a matrix, with each row representing an object.}
  \item{Y}{property that is to be modelled. In case of classification, Y
    is a matrix of zeros, with exactly one '1' in each row indicating
    the class. For prediction of continuous properties, Y is a vector. A
    combination is possible, too, but one then should take care of
    appropriate scaling.}
  \item{grid}{a grid for the representatives: see \code{somgrid}.}
  \item{rlen}{the number of times the complete data set will be
    presented to the network.}
  \item{alpha}{learning rate, a vector of two numbers indicating the
    amount of change. Default is to decline linearly from 0.05 to 0.01
    over \code{rlen} updates.}
  \item{radius}{the starting radius of the neighbourhood to be used for
    each update: the decrease is linear over \code{rlen} updates and
    reaches the value of 1 after one-third of the iterations.}
  \item{xweight}{the weight given to the X map in the
    calculation of distances for updating Y. Default is 0.5.}
  \item{toroidal}{if TRUE, the edges of the map are joined. Note
    that in a hexagonal toroidal map, the number of rows must be even.}
  \item{keep.data}{save data in return value.}
}

\value{
  an object of class "kohonen" with components
  \item{grid}{the grid, an object of class "somgrid".}
  \item{changes}{vector of mean average deviations from code vectors}
  \item{codes}{a matrix of code vectors.}
  \item{codeYs}{a matrix of Y values associated with each unit.}
  \item{toroidal}{whether a toroidal map is used.}
  \item{data}{data matrix, only if \code{keep.data} is TRUE.}
  \item{Y}{Y, only if \code{keep.data} is TRUE.}
}

\references{W.J. Melssen, R. Wehrens, and L.M.C. Buydens.
Chemom. Intell. Lab. Syst., in press.}
\author{Ron Wehrens}

\seealso{\code{\link{som}}, \code{\link{bdk}},
  \code{\link{plot.kohonen}}, \code{\link{predict.kohonen}}}
\examples{
### Simulated example
library(MASS)

csize <- 150
c1 <- mvrnorm(csize, mu=c(5,3,4),
              Sigma=matrix(c(2, 0,0,0, 2, 1, 0, 1, 2),3,3))
c2 <- mvrnorm(csize, mu=c(5.5, 3.5, 4.5),
              Sigma=matrix(c(2, 0,0,0, 2, 1, 0, 1, 2),3,3))
c3 <- mvrnorm(csize, mu=c(0,0,0),
              Sigma=matrix(c(2, 0,0,0, 2, 1, 0, 1, 2),3,3))

X <- rbind(c1, c2, c3)
classes <- c(rep(1, csize), rep(2, csize), rep(3, csize))

xyfmap <- xyf(X, classvec2classmat(classes), somgrid(8, 8, "hexagonal"))
plot(xyfmap, "prediction")

### Wine example
data(wines)
set.seed(7)

training <- sample(length(wine.classes), 120)
Xtraining <- scale(wines[training,])
xyf.wines <- xyf(Xtraining, classvec2classmat(wine.classes[training]),
  grid = somgrid(5, 5, "hexagonal"), rlen=100)
Xtest <- scale(wines[-training,],
               center = attr(Xtraining, "scaled:center"),
               scale = attr(Xtraining, "scaled:scale"))
xyf.prediction <- predict(xyf.wines, data=Xtest)
xyfclass.prediction <-
  classmat2classvec(xyf.wines$codeYs[xyf.prediction$unit.classif,])
table(wine.classes[-training], xyfclass.prediction)
}

\keyword{classif}
