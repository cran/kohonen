"map" <- function(x, ...)
{
  UseMethod("map")
}


"map.kohonen" <- function(x, newdata, ...)
{
  codes <- x$codes
  ncodes <- nrow(codes)
  np <- ncol(codes)
  nd <- nrow(newdata)
  
  res <- .C("VR_knn1",
            as.integer(ncodes),
            as.integer(nd),
            as.integer(np),
            as.double(codes),
            as.integer(1:ncodes),
            as.double(newdata),
            res = integer(nd),
            integer(ncodes+1),
            as.integer(ncodes),
            dists = double(nd),
            PACKAGE = "class")
  
  list(unit.classif = res$res, distances = sqrt(res$dists))
}
