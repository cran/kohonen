somKmeans <- function(codes, data, max.iter=20, verbose=FALSE)
{
  nd <- nrow(data)
  np <- ncol(data)
  ncodes <- nrow(codes)
  
  newclassif <- rep(0, nd)
  res <- .C("VR_knn1",
            as.integer(ncodes),
            as.integer(nd),
            as.integer(np),
            as.double(codes),
            as.integer(1:ncodes),
            as.double(data),
            res = integer(nd),
            integer(ncodes+1),
            as.integer(ncodes),
            double(nd),
            PACKAGE = "class")
  
  for (i in 1:max.iter) {
    oldclassif <- res$res
    
    for (j in 1:ncodes) {
      hits <- which(oldclassif == j)
      if (length(hits) > 0)
        codes[j,] <- colMeans(data[hits,,drop=FALSE])
    }

    res <- .C("VR_knn1",
              as.integer(ncodes),
              as.integer(nd),
              as.integer(np),
              as.double(codes),
              as.integer(1:ncodes),
              as.double(data),
              res = integer(nd),
              integer(ncodes+1),
              as.integer(ncodes),
              double(nd),
              PACKAGE = "class")
    newclassif <- res$res

    if (verbose)
      cat("\nIteration ", i, ": agreement ", sum(newclassif == oldclassif),
          " out of ", nd, " cases", sep="")

    if (all(newclassif == oldclassif)) break
  }

  if (i == max.iter)
    warning("Maximum number of iterations reached in kohonenKmeans")

  list(codes=codes, classif=newclassif, niter=i)
}
