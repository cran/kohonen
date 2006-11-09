\name{summary.kohonen}
\alias{summary.kohonen}
\alias{print.kohonen}
\title{Summary and print Methods for kohonen objects}
\description{
  Summary and print methods for \code{kohonen}}
\usage{
\method{summary}{kohonen}(object, \dots)
\method{print}{kohonen}(x, \dots)
}
\arguments{
  \item{x, object}{a \code{kohonen} object}
  \item{\dots}{Not used.}
}
\author{Ron Wehrens}
\seealso{\code{\link{som}}, \code{\link{xyf}}, \code{\link{bdk}}}
\examples{
data(wines)
xyf.wines <- xyf(scale(wines), classvec2classmat(wine.classes),
                 grid = somgrid(5, 5, "hexagonal"))
xyf.wines
summary(xyf.wines)
}
\keyword{classif}
