"map" <- function(x, ...)
{
  UseMethod("map")
}

## Map data to a trained supersom network using any of the layers, or any
## subset of them.

## Aug 2, 2016: Completely remodelled this function, now looks much
## simpler. Now only the winning units and corresponding distances
## are returned from C, which leads to much less memory consumption
## and hopefully faster processing. Especially important for large data.
## With the new structure a kohonen object will always have weights
## and whatmap elements
"map.kohonen" <- function(x,
                          newdata,
                          whatmap=NULL,
                          weights=NULL,
                          maxNA.fraction = 0L,
                          c.interface = c("c", "call", "rcpp"),
                          ...)
{
  ## ##########################################################################
  ## Get relevant info from kohonen object
  codes <- x$codes
  if (missing(newdata) & !is.null(x$data)) newdata <- x$data
  if (is.matrix(newdata)) newdata <- list(newdata) ## service to the user
  if (is.null(weights)) {
    weights <- x$weights
    useTrainingWeights <- TRUE
  } else {
    useTrainingWeights <- FALSE
  }
  if (is.null(whatmap)) whatmap <- x$whatmap
  dist.fcts <- x$dist.fcts
  dist.types <- x$dist.types

  ## ##########################################################################
  ## Check whatmap
  whatmap <- check.whatmap(newdata, whatmap)
  if (length(newdata) != length(codes) &&
      length(newdata) != length(whatmap))
    stop("Incorrect number of data layers in newdata")

  if (useTrainingWeights & any(weights[whatmap] < 1e-8))
    warning("Mapping new data using data layers not involved in training")

  ## ##########################################################################
  ## Apply whatmap
  newdata <- newdata[whatmap]
  dist.fcts <- dist.fcts[whatmap]
  dist.types <- dist.types[whatmap]
  codes <- codes[whatmap]
  weights.orig <- weights
  weights <- weights[whatmap]
  if (length(whatmap) == 1)
    weights <- 1
  if (sum(weights >= 1e-8) == 0)
    stop("Only weights of zero given")

  ## ##########################################################################
  ## Final preparations
  nvars <- sapply(codes, ncol)
  ncodes <- nrow(codes[[1]])

  if (any(factor.idx <- sapply(newdata, is.factor)))
    newdata[factor.idx] <- lapply(newdata[factor.idx], classvec2classmat)

  if (maxNA.fraction == 0L) {
    if (any(sapply(newdata, function(dt) any(is.na(dt)))))
      stop("maxNA.fraction equals 0, so no NAs allowed in the data!")
  } else {
    newdata <- check.data.na(newdata, maxNA.fraction, na.rm = FALSE)
  }

  ## Check that all matrices in the data set have the same number of objects
  nobjects <- unique(sapply(newdata, nrow))
  if (length(nobjects) > 1)
    stop("unequal numbers of objects in data list")

  if (maxNA.fraction > 0L) {
    nNA <- t(sapply(newdata,
                    function(x) apply(x, 1, function(y) sum(is.na(y)))))
  } else {
    nNA <- matrix(0, length(newdata), nobjects)
  }

  newdata <- matrix(unlist(newdata), ncol=nobjects, byrow=TRUE)
  codes <- matrix(unlist(codes), ncol=ncodes, byrow=TRUE)

  nmat <- length(whatmap)

  ## ##########################################################################
  ## Prepare C output variables
  winners <- rep(0L, nobjects)
  unit.distances <- rep(0, nobjects)

  ## ##########################################################################
  ## Go!
  c.interface <- match.arg(c.interface)
  switch(c.interface,
    c = {
      res <- .C("MapKohonen",
                data = as.double(newdata),
                nd = as.integer(nobjects),                 ## scalar
                nvar = as.integer(nvars),                   ## vector
                nNA = as.integer(nNA),
                codes = as.double(codes),
                ncodes = as.integer(ncodes),               ## scalar
                nmat = as.integer(nmat),                   ## scalar
                weights = as.double(weights),              ## vector
                distancefcts = dist.types,
                winners = as.integer(winners),             ## vector
                unitdistances = as.double(unit.distances), ## vector
                NAOK = TRUE,
                PACKAGE = "parkohonen")
    },
    call = {
      res <- .Call("CallMapKohonen",
                   data = as.double(newdata),
                   nd = as.integer(nobjects),
                   nNA = as.integer(nNA),
                   nvar = as.integer(nvars),
                   weights = as.double(weights),
                   codes = as.double(codes),
                   ncodes = as.integer(ncodes),
                   distancefcts = dist.types,
                   NAOK = TRUE,
                   PACKAGE = "parkohonen")
    },
    rcpp = {
      res <- RcppMap(data = newdata,
                     numVars = nvars,
                     numNAs = nNA,
                     codes = codes,
                     weights = weights,
                     distanceFunctions = dist.fcts)
    })

  list(unit.classif = res$winners + 1, ## C starts counting at 0...
       distances = res$unitdistances,
       whatmap = whatmap,
       weights = weights.orig)
}
