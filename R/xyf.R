"xyf" <- function(data, Y, grid = somgrid(), rlen = 100,
                  alpha = c(0.05, 0.01),
                  radius = quantile(nhbrdist, 0.67),
                  xweight = 0.5,
                  toroidal = FALSE, keep.data = TRUE)
{
  data <- as.matrix(data)
  
  nd <- nrow(data)
  nx <- ncol(data)

  if (is.vector(Y)) Y <- matrix(Y, ncol=1)
  ny <- ncol(Y)
  predict.type <- ifelse (all(rowSums(Y) == 1),
                          "class",
                          "continuous")

  ng <- nrow(grid$pts)
  xdists <- ydists <- rep(0, ng)  

  starters <- sample(1:nd, ng, replace = FALSE)
  
  init <- data[starters, , drop = FALSE]
  codeYs <- matrix(runif(ng*ny), ng, ny)
  codeYs <- sweep(codeYs, 1, rowSums(codeYs), FUN="/")
  codes <- init

  nhbrdist <- unit.distances(grid, toroidal)
  changes <- rep(0, rlen*2)

  if (predict.type == "continuous") {
    res <- .C("XYF_Eucl",
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
    res <- .C("XYF_Tani",
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
  codes <- res$codes
  codeYs <- res$codeYs
  dim(codes) <- dim(init)
  dim(codeYs) <- c(ng, ny)
  colnames(codeYs) <- colnames(Y)

  classif <- predict.kohonen(list(data = data, codes = codes))$unit.classif
  
  if (keep.data) {
    structure(list(data = data, Y = Y, predict.type = predict.type,
                   grid = grid, codes = codes, codeYs = codeYs,
                   changes = changes, toroidal = toroidal, classif = classif),
              class = "kohonen")
  } else {
    structure(list(predict.type = predict.type,
                   grid = grid, codes = codes, codeYs = codeYs,
                   changes = changes, toroidal = toroidal, classif = classif), 
              class = "kohonen")
  }
}

