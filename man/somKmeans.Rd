\name{somKmeans}
\alias{somKmeans}
\title{Refinement of codebook vectors by kmeans}
\description{
  Do kmeans of the codebook vectors, which is equivalent to the
  SOM algorithm when the size of the neighbourhood is smaller than
  1. Not meant to be called directly. 
}
\usage{
somKmeans(codes, data, max.iter = 20, verbose = FALSE)
}
\arguments{
  \item{codes}{codebook matrix for the X space: one row per unit.}
  \item{data}{training data, usually.}
  \item{max.iter}{default is 20, but convergence usually is much
    faster. If more iterations are needed, a warning is printed.}
  \item{verbose}{print more info to the screen.}
}
\value{
  Returns a list with components
  \item{codes}{new codebook matrix.}
  \item{classif}{classification vector.}
  \item{niter}{number of iterations needed for kmeans.}
}
\author{Ron Wehrens}
\keyword{classif}
