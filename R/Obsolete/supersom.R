supersom <- function(data,
                     grid = somgrid(),
                     rlen = 100,
                     alpha = c(0.05, 0.01),
                     radius = quantile(nhbrdist, 0.67),
                     whatmap = NULL,
                     weights = 1,
                     maxNA.fraction = 0L,
                     keep.data = TRUE,
                     dist.fcts = NULL,
                     mode = c("online", "batch", "pbatch", "sugar"),
                     c.interface = c("c", "call", "rcpp"),
                     init)
{
  ## ##########################################################################
  ## Check data
  if (!is.list(data) & (is.vector(data) | is.factor(data) | is.matrix(data)))
    data <- list(data)
  ## if there are any vectors, convert them to one-column matrices
  vctr.idx <- sapply(data, is.vector)
  if (any(vctr.idx))
    for (i in which(vctr.idx))
      data[[i]] <- matrix(data[[i]], ncol = 1)
  
  orig.data <- data

  ## Check whether data is a list of data matrices or factors
  if (!is.list(data) | !all(sapply(data, is.matrix) | sapply(data, is.factor)))
    stop("data should be a list of data matrices or factors")

  ## Convert vectors to one-column matrices in layers
  if (any (vector.idx <- sapply(data, is.vector))) {
    warning("Converting vectors to one-column matrices in layers",
            paste(vector.idx, sep = ", "),
            "if you meant factors please present the data as such")
    data[vector.idx] <- sapply(data, function(x) matrix(x, ncol = 1))
  }
  if (any(factor.idx <- sapply(data, is.factor)))
    data[factor.idx] <- lapply(data[factor.idx], classvec2classmat)

  ## Check whether data is numeric
  if (!all(sapply(data, is.numeric)))
    stop("Argument data should be numeric")

  ## ##########################################################################
  ## Check weights
  if (length(weights) == 1) {
    weights <- rep(weights, length(data))
  } else if (length(weights) != length(data)) {
    stop("Mismatch between the number of data layers and the number of weights")
  }

  ## ##########################################################################
  ## Check distance functions
  prefabDists <- c("sumofsquares", "euclidean", "manhattan", "tanimoto")
  defaultDist <- 1      ## use sumofsquares by default, for non-class vars
  defaultClassDist <- 4 ## use tanimoto for class variables

  if (is.null(dist.fcts)) {
    ## default is to use tanimoto distances for classes, and sumofsquares
    ## distances, the first element of prefabDists, for anything
    ## else. We use all data layers here, and later select the ones
    ## relevant using the whatmap argument
    factorLayers <- sapply(data,
                           function(dt)
                             is.factor(dt) ||
                             (is.matrix(dt) &&
                              min(dt) >= 0 && max(dt) <= 1 &&
                              all(abs(rowSums(dt) - 1) < 1e-8 )))
    dist.fcts <- as.list(rep(prefabDists[defaultDist], length(data)))
    if (any(factorLayers))
      dist.fcts[[which(factorLayers)]] <- prefabDists[defaultClassDist]
  } else {
    if (length(dist.fcts) == 1)
      dist.fcts <- as.list(rep(dist.fcts, length(data)))
  }

  if (length(dist.fcts) != length(data)) {
    stop(c("The number of distance functions should be one, ",
           "or equal to the number of data layers"))
  }

  ## we should also allow mixtures of text and explicit function pointers...
  dist.classes <- sapply(dist.fcts, class)
  dist.names <- rep("User-defined", length(dist.fcts))
  dist.types <- rep(NA, length(dist.fcts))
  if (!all(dist.classes %in% c("character", "externalptr")))
    stop(c("Wrong distance function specification: ",
           "should be specified as vectors of strings ",
           "or explicit function pointers"))
  if (any(char.idx <- dist.classes == "character")) {
    ## first convert to a factor, and then to a pointer...
    dist.names[char.idx] <- dist.classes[char.idx]
    dist.types[char.idx] <- as.integer(factor(dist.fcts[char.idx],
                                              levels = c(prefabDists)))
    dist.fcts[char.idx] <-
      sapply(lapply(dist.fcts[char.idx], factor, levels = prefabDists),
             function(type) CreateDistanceFunctionXPtr(type, maxNA.fraction > 0L))
  }
  orig.dist.fcts <- dist.fcts
  orig.dist.types <- dist.types

  ## ##########################################################################
  ## Check radius update parameters
  nhbrdist <- unit.distances(grid)
  if (length(radius) == 1)
    radius <- abs(radius) * c(1, -1)

  ## ##########################################################################
  ## Check whatmap
  whatmap <- check.whatmap(data, whatmap)
  whatmap <- whatmap[weights[whatmap] != 0]
  nmat <- length(whatmap)
  if (nmat == 0) {
    stop("No data layers left - check whatmap and weights arguments")
  }

  ## ##########################################################################
  ## Apply whatmap
  data <- data[whatmap]
  dist.fcts <- dist.fcts[whatmap]
  dist.types <- dist.types[whatmap]
  weights <- weights[whatmap]

  ## Check whatmapped weights
  if (abs(sum(weights)) < .Machine$double.eps) {
    stop("weights sum to zero")
  } else {
    weights <- weights / sum(weights)
  }
  orig.weights <- rep(0, length(orig.data))
  orig.weights[whatmap] <- weights

  ## ##########################################################################
  ## Check the number of NAs, remove rows or columns with too many
  if (maxNA.fraction == 0L) {
    if (any(sapply(data, function(dt) any(is.na(dt)))))
      stop("maxNA.fraction equals 0, so no NAs allowed in the data!")
  } else {
    data <- check.data.na(data, maxNA.fraction, na.rm = TRUE)
  }

  nvar <- sapply(data, ncol)
  nobjects <- unique(sapply(data, nrow))
  if (length(nobjects) > 1)
    stop("unequal numbers of objects in data list")
  
  if (maxNA.fraction > 0L) {
    nNA <- t(sapply(data, function(x) apply(x, 1, function(y) sum(is.na(y)))))
  } else {
    nNA <- matrix(0, length(data), nobjects)
  }

  data.matrix <- matrix(unlist(data), ncol = nobjects, byrow = TRUE)

  ## ##########################################################################
  ## Get or create initial codebooks
  ncodes <- nrow(grid$pts)
  if (missing(init)) {
    starters <- sample(1:nobjects, ncodes, replace = FALSE)
    init.matrix <- data.matrix[,starters]
  } else {
    ## Check length and dimensions
    if (!is.list(init) & (is.vector(init) | is.factor(init) | is.matrix(init)))
      init <- list(init)
    if (length(init) != nmat)
      stop("Incorrect number of initialization matrices")
    if (!all(sapply(init, nrow) == ncodes))
      stop("Incorrect number of objects in initalization matrices")
    if (!all(sapply(init, ncol) == nvar)) {
      if (maxNA.fraction == 0L) {
        stop("Incorrect number of variables in initialization matrices")
      } else {
        stop("Incorrect number of variables in initialization matrices, ",
             "maybe due to the removal of columns because of NAs")
      }
    }
    init.matrix <- matrix(unlist(init), ncol = ncodes, byrow = TRUE)
  }

  ## Skip this if no NAs are allowed - that has already been checked
  ## note that these matrices contain objects in the columns, and
  ## variables in the rows
  if (maxNA.fraction > 0L) {
    ## Fill in any NAs with random draws (nonNA) from the
    ## corresponding columns in data
    nastarters <- which(is.na(init.matrix), arr.ind = TRUE)
    if (nrow(nastarters) > 0) {
      for (i in unique(nastarters[,1])) { ## rows
        idx <- which(nastarters[,1] == i) ## columns
        imputers <- sample(data.matrix[i, !is.na(data.matrix[i,])],
                           length(idx),
                           replace = TRUE)
        init.matrix[i, nastarters[idx,2]] <- imputers
      }
    }
  }

  ## ##########################################################################
  ## Go!
  c.interface <- match.arg(c.interface)
  switch(c.interface,
    c = {
      changes <- rep(0, rlen * nmat)
      mode <- match.arg(mode)
      call.func <- switch(mode,
                          online = "Supersom",
                          batch = "BatchSupersom",
                          pbatch = "ParallelBatchSupersom")
      res <- .C(call.func,
                data = as.double(data.matrix),
                codes = as.double(init.matrix),
                nhbrdist = as.double(nhbrdist),
                alpha = as.double(alpha),
                radii = as.double(radius),
                weights = as.double(weights),
                changes = as.double(changes),
                nobjects = as.integer(nobjects),
                nmat = as.integer(nmat),
                nvar = as.integer(nvar),
                nNA = as.integer(nNA),
                ncodes = as.integer(ncodes),
                rlen = as.integer(rlen),
                distancefcts = dist.types,
                neighbourhoodfct = as.integer(grid$neighbourhood.fct),
                NAOK = TRUE,
                PACKAGE = "parkohonen")
        changes <- matrix(res$changes, ncol = nmat, byrow = TRUE)
        colnames(changes) <- names(data)
        mycodes <- t(matrix(res$codes, nrow = ncodes, byrow = TRUE))
    },
    call = {
      mode <- match.arg(mode)
      call.func <- switch(mode,
                          online = "CallSupersom")
      res <- .Call(call.func,
                   data = as.double(data.matrix),
                   codes = as.double(init.matrix),
                   nhbrdist = as.double(nhbrdist),
                   alpha = as.double(alpha),
                   radii = as.double(radius),
                   weights = as.double(weights),
                   numObjects = as.integer(nobjects),
                   nmat = as.integer(nmat),
                   nvar = as.integer(nvar),
                   nNA = as.integer(nNA),
                   ncodes = as.integer(ncodes),
                   rlen = as.integer(rlen),
                   distancefcts = dist.types,
                   neighbourhoodfct = as.integer(grid$neighbourhood.fct),
                   NAOK = TRUE,
                   PACKAGE = "parkohonen")
      changes <- matrix(res$changes, ncol = nmat, byrow = TRUE)
      colnames(changes) <- names(data)
      mycodes <- t(matrix(res$codes, nrow = ncodes, byrow = TRUE))
    },
    rcpp = {
      mode <- match.arg(mode)
      call.func <-
        switch(mode,
               online = match.fun(parkohonen:::RcppSupersom),
               sugar = match.fun(parkohonen:::RcppSugarSupersom),
               batch = match.fun(parkohonen:::RcppBatchSupersom),
               pbatch = match.fun(parkohonen:::RcppParallelBatchSupersom))
          res <- call.func(data = data.matrix,
                           codes = init.matrix,
                           numVars = nvar,
                           weights = weights,
                           distanceFunctions = dist.fcts,
                           numNAs = nNA,
                           neighbourhoodDistances = nhbrdist,
                           neighbourhoodFct = as.integer(grid$neighbourhood.fct),
                           alphas = alpha,
                           radii = radius,
                           numEpochs = rlen)
        changes <- matrix(res$changes, ncol = nmat, byrow = TRUE)
        colnames(changes) <- names(data)
        mycodes <- res$codes
    })

  ## ##########################################################################
  ## Format the codes
  layerID <- rep(1:nmat, nvar)
  mycodes2 <- split(as.data.frame(mycodes), layerID)
  mycodes3 <- lapply(mycodes2, function(x) t(as.matrix(x)))
  codes <- vector(length(orig.data), mode = "list")
  names(codes) <- names(orig.data)
  codes[whatmap] <- mycodes3
  for (ii in seq(along = whatmap))
    colnames(codes[[ whatmap[ii] ]]) <- colnames(data[[ii]])

  ## ##########################################################################
  ## Prepare results
  if (keep.data) {
    mapping <- map.kohonen(list(codes = codes,
                                dist.fcts = orig.dist.fcts,
                                dist.types = orig.dist.types),
                           newdata = orig.data,
                           whatmap = whatmap,
                           weights = orig.weights,
                           maxNA.fraction = maxNA.fraction,
                           c.interface = c.interface)
    structure(list(data = orig.data,
                   unit.classif = mapping$unit.classif,
                   distances = mapping$distances,
                   grid = grid,
                   codes = codes,
                   changes = changes,
                   alpha = alpha,
                   radius = radius,
                   weights = orig.weights,
                   whatmap = whatmap,
                   dist.fcts = orig.dist.fcts,
                   dist.types = orig.dist.types,
                   dist.names = dist.names),
              class = "kohonen")
  } else {
    structure(list(grid = grid,
                   codes = codes,
                   changes = changes,
                   alpha = alpha,
                   radius = radius,
                   weights = orig.weights,
                   whatmap = whatmap,
                   dist.fcts = orig.dist.fcts,
                   dist.names = dist.names,
                   dist.types = dist.types),
              class = "kohonen")
  }
}
