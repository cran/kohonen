\name{nirTshift}
\alias{nirTshift}
\alias{nirXtr}
\alias{nirYtr}
\alias{nirXtest}
\alias{nirYtest}
\docType{data}
\title{Near-infrared data with temperature effects}
\description{Near-infrared spectra of ternary mixtures of ethanol, water
  and iso-propanol, measured at five different temperatures (30, 40,
  ..., 70 degrees Centigrade). A training set of 65 mixtures (13 spectra
  per temperature) is available, as well as a test set of 30 spectra (5
  per temperature). Spectra are in variables \code{nirXtrain} and
  \code{nirXtest}; concentrations of the three compounds in \code{nirYtrain}
  and \code{nirYtest}.} 
\references{F. Wülfert , W.Th. Kok, A.K. Smilde: Anal. Chem. 1998, 1761-1767}
\examples{
  data(nirTshift)

  ## First model concentrations, all at once
  set.seed(3)
  nirnet <- xyf(data = nirXtr, Y = nirYtr, xweight=.75,
                grid = somgrid(6, 6, "hexagonal"), rlen=500)
  par(mfrow = c(1,3))
  plot(nirnet, "predict")

  ## Concentrate on water, compound 2: 
  par(mfrow = c(1,2))
  set.seed(3)
  nirnet <- xyf(data = nirXtr, Y = nirYtr[,2], xweight=.75,
                grid = somgrid(6, 6, "hexagonal"), rlen=500)
  plot(nirnet, "predict", main="Prediction of water content")
  ## Plot temperatures as circles
  temps <- rep(seq(30, 70, by=10), times=rep(13,5))
  symbols(nirnet$grid$pts[nirnet$unit.classif,] +
          matrix(rnorm(nrow(nirYtr)*2, sd=.1), ncol=2),
	  circles = (temps - 20)/250,
	  inches = FALSE, add = TRUE)

  ## Model temperatures	  
  set.seed(3)
  nirnet2 <- xyf(data = nirXtr, Y = classvec2classmat(temps), xweight=.3,
                 grid = somgrid(6, 6, "hexagonal"), rlen=500)
  plot(nirnet2, "predict", palette.name = rainbow,
       main="Prediction of temperatures")
  ## Plot concentrations of water as circles
  symbols(nirnet2$grid$pts[nirnet2$unit.classif,] +
          matrix(rnorm(nrow(nirYtr)*2, sd=.1), ncol=2),
	  circles = 0.05 + 0.4 * nirYtr[,2],
	  inches = FALSE, add = TRUE)
}
\keyword{datasets}
