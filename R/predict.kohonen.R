"predict.kohonen" <- function(object, newdata, trainX, trainY,
                              unit.predictions = NULL, threshold = 0,
                              whatmap = NULL, weights = 1, ...)
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
    mapping <- map(object, newdata, whatmap, weights)$unit.classif
  ### FIXME: is hetvolgende OK?
  ## whatmap <- mapping$whatmap
  
  ## Find the value for each unit. For unsupervised maps, we should
  ## have trainX and trainY values. Argument unit.predictions may
  ## override the predictions per unit.
  if (is.null(unit.predictions)) {
    if (object$method %in% c("xyf", "bdk")) {
      unit.predictions <- object$codes$Y
    } else {
      ## If whatmap is given, and not all layers present in the map are
      ## included, the excluded layers are interpreted as being the ones
      ## for which a prediction is wanted.
      if (object$method == "supersom" & !is.null(whatmap)) {
        whatmap <- check.whatmap(object, whatmap)
        if (length(whatmap) < length(object$data))
          unit.predictions <- object$codes[-whatmap]    
      } else {
        if (missing(trainY))
          stop("For unsupervised forms of mapping, trainY is required")
        if (is.list(trainY))
          stop("Prediction for trainY lists not implemented")
        
        if (is.vector(trainY))
          trainY <- matrix(trainY, ncol = 1)
        nY <- ncol(trainY)
        
        trainingMapping <- NULL
        if (missing(trainX) & !is.null(object$data)) {
          trainX <- object$data
          trainingMapping <- object$unit.classif
        }

        nX <- ifelse(is.list(trainX),
                     nrow(trainX[[1]]),
                     nrow(trainX))
        if (nX != nrow(trainY))
          stop("Unequal number of rows in trainX and trainY")

        ## Find mapping for training data
        if (is.null(trainingMapping))
          trainingMapping <- map(object, trainX)$unit.classif

        ## Find unit.predictions for training data
        unit.predictions <- matrix(NA, nrow(object$grid$pts), nY)
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
      }
    }
  }
      

  if (!is.null(object$contin) && !object$contin) {
    prediction <- # xyf or bdk for categorical variables
      classmat2classvec(unit.predictions, threshold=threshold)[mapping]
  } else {
    if (is.list(unit.predictions)) { # supersom
      prediction <- sapply(unit.predictions,
                           function(x) x[mapping,])
    } else {
      prediction <- unit.predictions[mapping,]
    }
  }
  
  list(prediction = prediction, unit.classif = mapping,
       unit.predictions = unit.predictions)
}
