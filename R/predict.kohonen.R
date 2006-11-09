"predict.kohonen" <- function(object, newdata, trainX, trainY,
                              unit.predictions, threshold = 0, ...)
{
  mapping <- NULL
  if (missing(newdata)) {
    if (!is.null(object$data)) {
      newdata <- object$data
      mapping <- object$unit.classif # perhaps NULL
    } else {
      stop("No data given with which to predict")
    }
  }

  ## Find the mapping for the new data
  if (is.null(mapping))
    mapping <- map(object, newdata)$unit.classif

  ## Find the value for each unit. For unsupervised maps, we should
  ## have trainX and trainY values.
  if (missing(unit.predictions)) {
    ## this argument may override the predictions per unit that are
    ## the result of supervised training
    if (is.null(object$predict.type)) { # unsupervised mapping
      if (missing(trainY))
        stop("For unsupervised maps, argument trainY is required")
      if (is.vector(trainY))
        trainY <- matrix(trainY, ncol = 1)
      nY <- ncol(trainY)
      
      trainingMapping <- NULL
      if (missing(trainX) & !is.null(object$data)) {
        trainX <- object$data
        trainingMapping <- object$unit.classif
      }
      if (nrow(trainX) != nrow(trainY))
        stop("Number of rows in trainY not equal to number of rows in trainX")

      ## Find mapping for training data
      if (is.null(trainingMapping))
        trainingMapping <- map(object, trainX)$unit.classif

      ## Find unit.predictions for training data
      unit.predictions <- matrix(NA, nrow(object$codes), nY)
      huhn <- aggregate(trainY, by = list(cl = trainingMapping),
                        mean)
      unit.predictions[sort(as.numeric(levels(huhn[, 1]))),] <-
        as.matrix(huhn[, -1])

      ## Prediction for empty units
      nas <- which(apply(unit.predictions, 1, function(x) all(is.na(x))))
      nhbrdist <- unit.distances(object$grid, !is.null(object$toroidal))
      for (i in seq(along = nas)) {
        unit.predictions[nas[i], ] <-
          mean(unit.predictions[nhbrdist[i,] == 1, ], na.rm = TRUE)
      }
      
      colnames(unit.predictions) <- colnames(trainY)
    } else { ## For supervised maps, things are easier.
      unit.predictions <- object$codeYs
    }
  }

  ## FIXME: what about continuous prediction???
  classes <- classmat2classvec(unit.predictions, threshold=threshold)
  
  list(classif = classes[mapping], unit.classif = mapping)
}
