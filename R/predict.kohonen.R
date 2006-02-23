"predict.kohonen" <- function(object, ...)
{
  aux <- list(...)
  classif <- aux$classif
  data <- aux$data
  if (is.null(data)) data <- object$data
  
  codes <- object$codes
  ncodes <- nrow(codes)
  np <- ncol(codes)
  nd <- nrow(data)
  
  if (is.null(classif)) {
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
              dists = double(nd),
              PACKAGE = "class")

    classif <- res$res
    dists <- sqrt(res$dists)
  } else {
    res <- .C("kohonenAssign",
              as.double(data),
              as.double(codes),
              as.integer(ncodes),
              as.integer(nd),
              as.integer(np),
              as.integer(classif-1),
              dists = double(nd),
              PACKAGE = "kohonen")
    dists <- res$dists
  }

  if (!is.null(object$codeYs)) {
    threshold <- aux$threshold
    if (is.null(threshold)) threshold <- 0
    classes <- classmat2classvec(object$codeYs, threshold=threshold)

    list(classif = classes[classif], unit.classif = classif,
         distances = dists)
  } else {
    list(unit.classif = classif, distances = dists)
  }
}
