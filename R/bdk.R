### Version 2.0: codeYs are no longer returned, but codes are a list
### with two elements: X and Y.

"bdk" <- function(data, Y, grid = somgrid(), rlen = 100,
                  alpha = c(0.05, 0.01),
                  radius = quantile(nhbrdist, 0.67),
                  xweight = 0.75, contin = !(all(rowSums(Y) == 1)),
                  toroidal = FALSE, keep.data = TRUE)
{
  data <- as.matrix(data)
  
  nd <- nrow(data)
  nx <- ncol(data)
  
  if (is.vector(Y)) Y <- matrix(Y, ncol=1)
  ny <- ncol(Y)

  ng <- nrow(grid$pts)
  xdists <- ydists <- rep(0, ng)  
  
  starters <- sample(1:nd, ng, replace = FALSE)
  init <- data[starters, , drop = FALSE]
  codes <- init
  if (!contin) {
    ## rescale to .25 - .75 in order to make class transitions easier
    codeYs <- 0.5 + 0.5*(Y[starters,] - 0.5)
  } else {
    codeYs <- Y[starters,]
  }
  
  nhbrdist <- unit.distances(grid, toroidal)
  changes <- rep(0, rlen*2)

  if (contin) {
    res <- .C("BDK_Eucl",
              data = as.double(data),
              Ys = as.double(Y),
              codes = as.double(codes),
              codeYs = as.double(codeYs),
              nhbrdist = as.double(nhbrdist),
              alpha = as.double(alpha),
              radius = as.double(radius),
              xweight = as.double(xweight),
              changes = as.double(changes),
              xdists = as.double(xdists),
              ydists = as.double(ydists),
              n = as.integer(nd),
              px = as.integer(nx),
              py = as.integer(ny),
              ncodes = as.integer(ng),
              rlen = as.integer(rlen),
              PACKAGE = "kohonen")
  } else {
    res <- .C("BDK_Tani",
              data = as.double(data),
              Ys = as.double(Y),
              codes = as.double(codes),
              codeYs = as.double(codeYs),
              nhbrdist = as.double(nhbrdist),
              alpha = as.double(alpha),
              radius = as.double(radius),
              xweight = as.double(xweight),
              changes = as.double(changes),
              xdists = as.double(xdists),
              ydists = as.double(ydists),
              n = as.integer(nd),
              px = as.integer(nx),
              py = as.integer(ny),
              ncodes = as.integer(ng),
              rlen = as.integer(rlen),
              PACKAGE = "kohonen")
  }
  
  changes <- matrix(res$changes, ncol=2)
  codes <- list(X = matrix(res$codes, nrow(init), ncol(init)),
                Y = matrix(res$codeYs, ng, ny))
  colnames(codes$Y) <- colnames(Y)

  if (keep.data) {
    mapping <- map.kohonen(list(codes = codes), newdata = data, whatmap = 1)
    
    structure(list(data = data, Y = Y, contin = contin,
                   grid = grid, codes = codes,
                   changes = changes, toroidal = toroidal,
                   unit.classif = mapping$unit.classif,
                   distances = mapping$distances, method="bdk"),
              class = "kohonen")
  } else {
    structure(list(contin = contin,
                   grid = grid, codes = codes,
                   changes = changes, toroidal = toroidal, method="bdk"),
              class = "kohonen")
  }
}

