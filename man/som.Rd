\name{som}
\alias{som}
\title{Kohonen's self-organising maps}
\description{
  Self-organising maps for mapping high-dimensional spectra or patterns
  to 2D; the Euclidean distance is used. Modelled after the SOM function in
  package \code{class}.
}

\usage{
som(data, grid=somgrid(), rlen = 100, alpha = c(0.05, 0.01),
    radius = quantile(nhbrdist, 0.67), init,
    toroidal = FALSE, FineTune = TRUE, keep.data = TRUE)
}

\arguments{
  \item{data}{a matrix, with each row representing an object.}
  \item{grid}{a grid for the representatives: see 'somgrid'.}
  \item{rlen}{the number of times the complete data set will be
    presented to the network.}
  \item{alpha}{learning rate, a vector of two numbers indicating the
    amount of change. Default is to decline linearly from 0.05 to 0.01
    over \code{rlen} updates.}
  \item{radius}{the initial radius of the neighbourhood to be used for
    each update: it decreases linearly to 1 over \code{rlen}
    updates. The default is to start with a value that covers 2/3 of all
    units.} 
  \item{init}{the initial representatives, represented as a matrix.  If
    missing, chosen (without replacement) randomly from 'data'.}
  \item{toroidal}{if TRUE, the edges of the map are joined. Note
    that in a hexagonal toroidal map, the number of rows must be even.}
  \item{FineTune}{apply kmeans for fine-tuning the codebook vectors.}
  \item{keep.data}{save data in return value.}
}

\value{
  an object of class "kohonen" with components
  \item{grid}{the grid, an object of class "somgrid".}
  \item{changes}{vector of mean average deviations from code vectors.}
  \item{codes}{a matrix of code vectors.}
  \item{classif}{winning units for all data objects.}
  \item{toroidal}{whether a toroidal map is used.}
  \item{data}{data matrix, only if \code{keep.data} is TRUE.}
  \item{type}{the type of som, here "som"}
}

\references{W.J. Melssen, R. Wehrens, and L.M.C. Buydens.
Chemom. Intell. Lab. Syst., 83, 99-113 (2006)..}
\author{Ron Wehrens}

\seealso{\code{\link{xyf}}, \code{\link{bdk}}, \code{\link{plot.kohonen}}}
\examples{
data(wines)
set.seed(7)

training <- sample(length(wine.classes), 120)
Xtraining <- scale(wines[training,])
som.wines <- som(Xtraining, grid = somgrid(5, 5, "hexagonal"))

Xtest <- scale(wines[-training,],
               center = attr(Xtraining, "scaled:center"),
               scale = attr(Xtraining, "scaled:scale"))

som.prediction <-
  predict(som.wines, newdata = Xtest,
          trainX = Xtraining,
          trainY = classvec2classmat(wine.classes[training]))
table(wine.classes[-training], som.prediction$classif)
}
\keyword{classif}
