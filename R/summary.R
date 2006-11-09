summary.kohonen <- function(object, ...)
{
  cat(object$method, " map of size ",
      object$grid$xdim, "x", object$grid$ydim,
      " with a ", object$grid$topo,
      if (object$toroidal) "toroidal", " topology.", sep="")
  if (!is.null(object$data)) {
    cat("\nTraining data included; dimension is",
        nrow(object$data), "by", ncol(object$data))
    if (!is.null(object$Y)) {
      cat("\nDimension of Y:", nrow(object$Y), "by", ncol(object$Y))
    }
    cat("\nMean distance to the closest unit in the map:",
        mean(object$distances))
  }
  if (!is.null(object$predict.type)) {
    cat("\nPrediction type:",
        ifelse(object$predict.type == "class", "classification", "regression"))
  }
  cat("\n")
}

print.kohonen <- function(x, ...)
{
  cat(x$method, " map of size ", x$grid$xdim, "x", x$grid$ydim,
      " with a ", x$grid$topo, if (x$toroidal) "toroidal",
      " topology.", sep="")
  if (!is.null(x$data))
    cat("\nTraining data included.")
  cat("\n")
}

