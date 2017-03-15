## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ------------------------------------------------------------------------
code <- ' 
  #include <Rcpp.h>
  	     
  typedef double (*DistanceFunctionPtr)(double *, double *, int, int);
  
  /*
    * My own custom euclidean distance function.
  */
    double myDistanceFunction(double *data, double *codes, int n, int nNA) {
      if (nNA == n) {
        return NA_REAL;
      }
      double tmp, d = 0.0;
      int i;
      for (i = 0; i < n; i++) {
        if (!std::isnan(data[i])) {
          tmp = data[i] - codes[i];
          d += tmp * tmp;
        }
      }
      if (nNA > 0) {
        d *= n / (n - nNA);
      }
      d = sqrt(d);
      return d;
    }
  
  /*
    * This function creates the external pointer that can be passed to the Kohonen
  * package.
  */
    // [[Rcpp::export]]
  Rcpp::ExpressionVector createMyDistanceFunctionXPtr(int n = 1){
    Rcpp::ExpressionVector distanceFunctions(n);
    for (int i = 0; i < n; i++) {
      distanceFunctions[i] = Rcpp::XPtr<DistanceFunctionPtr>(new DistanceFunctionPtr(myDistanceFunction));
    }
    return distanceFunctions;
  }
'

## ------------------------------------------------------------------------
library(Rcpp)
library(kohonen)
sourceCpp(code = code)
data(yeast)
whatmap <- 3:6
distance.fcts <- createMyDistanceFunctionXPtr(length(yeast))
supersom.yeast <- supersom(yeast,
                           whatmap = whatmap,
                           grid = somgrid(10, 10, "rectangular"),
                           dist.fcts = "euclidean",
			   maxNA.fraction = 0.5
                  )

## ---- echo=FALSE---------------------------------------------------------
plot(supersom.yeast, type = "quality")

## ---- echo=FALSE---------------------------------------------------------
plot(supersom.yeast, type = "counts")

## ---- echo=FALSE---------------------------------------------------------
plot(supersom.yeast, type = "changes")

