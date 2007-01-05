"som" <- function (data, grid = somgrid(),
                   rlen = 100, alpha = c(0.05, 0.01),
                   radius = quantile(nhbrdist, 0.67),
                   init, toroidal = FALSE, keep.data = TRUE)
{
  data <- as.matrix(data)
  nd <- nrow(data)
  
  ng <- nrow(grid$pts)
  if (missing(init))
    init <- data[sample(1:nd, ng, replace = FALSE), , drop = FALSE]
  codes <- init
  nhbrdist <- unit.distances(grid, toroidal)

  changes <- rep(0, rlen)
  
  res <- .C("SOM_online",
            data = as.double(data),
            codes = as.double(codes),
            nhbrdist = as.double(nhbrdist),
            alpha = as.double(alpha),
            radii = as.double(radius),
            changes = as.double(changes),
            n = as.integer(nrow(data)),
            p = as.integer(ncol(data)),
            ncodes = as.integer(nrow(init)),
            rlen = as.integer(rlen),
            PACKAGE = "kohonen")

  changes <- matrix(res$changes, ncol=1)
  codes <- res$codes
  dim(codes) <- dim(init)
  colnames(codes) <- colnames(init)

  if (keep.data) {
    mapping <- map.kohonen(list(codes = codes), newdata = data)
    
    structure(list(data = data, grid = grid, codes = codes,
                   changes = changes, toroidal = toroidal,
                   unit.classif = mapping$unit.classif,
                   distances = mapping$distances, method="som"),
              class = "kohonen")
  } else {
    structure(list(grid = grid, codes = codes, changes = changes, 
                   toroidal = toroidal, method="som"),
              class = "kohonen")
  }
}
