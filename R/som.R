"som" <- function (data, grid = somgrid(),
                   rlen = 100, alpha = c(0.05, 0.01),
                   radius = quantile(nhbrdist, 0.67) * c(1, -1),
                   init, toroidal = FALSE,
                   n.hood, keep.data = TRUE)
{
  if (!is.numeric(data))
    stop("Argument data should be numeric")
  data <- as.matrix(data)
  
  nd <- nrow(data)
  ng <- nrow(grid$pts)
  if (missing(init)) {
    init <- data[sample(1:nd, ng, replace = FALSE), , drop = FALSE]
  } else {
    init <- as.matrix(init)
    if (nrow(init) != ng |
        ncol(init) != ncol(data) |
        !is.numeric(init))
      stop("incorrect init matrix supplied")
  }
  codes <- init

  if (missing(n.hood)) {
    n.hood <- switch(grid$topo,
                     hexagonal = "circular",
                     rectangular = "square")
  } else {
    n.hood <- match.arg(n.hood, c("circular", "square"))
  }
  grid$n.hood <- n.hood
  nhbrdist <- unit.distances(grid, toroidal)
  
  if (length(radius) == 1)
    radius <- sort(radius * c(1, -1), decreasing = TRUE)

  changes <- rep(0, rlen)
  res <- .C("SOM_online",
            data = as.double(data),
            codes = as.double(codes),
            code_snapshots = rep(as.double(codes), rlen*nrow(data)),
            code_snapshot_datum =as.integer(rep(1, rlen*nrow(data))),
            nhbrdist = as.double(nhbrdist),
            alpha = as.double(alpha),
            radii = as.double(radius),
            changes = as.double(changes),
            n = as.integer(nrow(data)),
            p = as.integer(ncol(data)),
            ncodes = as.integer(nrow(init)),
            rlen = as.integer(rlen),
            PACKAGE = "itsakohonen")

  changes <- matrix(res$changes, ncol=1)
  codes <- res$codes
  dim(codes) <- dim(init)
  colnames(codes) <- colnames(init)
  code_snapshots <- res$code_snapshots
  # itsakettle, first c(codes matrix size, number of iter size) gives
  # a matrix with columns of each code snapshot. Then 
  # c(codes rows, codes col, number of iter size)
  dim(code_snapshots) <- c(dim(init), rlen*nrow(data))
  
  code_snapshots_datum <- res$code_snapshot_datum + 1
  
  if (keep.data) {
    mapping <- map.kohonen(list(codes = codes), newdata = data)
    
    structure(list(data = data, grid = grid, codes = codes,
                   changes = changes, alpha = alpha,
                   radius = radius, toroidal = toroidal,
                   unit.classif = mapping$unit.classif,
                   code_snapshots=code_snapshots,
                   code_snapshots_datum=code_snapshots_datum,
                   distances = mapping$distances, method="som"),
              class = "kohonen")
  } else {
    structure(list(grid = grid, codes = codes, changes = changes,
                   alpha = alpha, radius = radius,
                   code_snapshots=code_snapshots,
                   code_snapshots_datum=code_snapshots_datum,
                   toroidal = toroidal, method="som"),
              class = "kohonen")
  }
}
