### Calculate distances in a Kohonen map. Crude and
### slow implementation, but hey. When training a map takes 20 hours,
### who cares about a few seconds?

unit.distances <- function(grid, toroidal)
{
  if (!toroidal) {
    if (grid$topo == "hexagonal")
      return(as.matrix(stats::dist(grid$pts)))
    if (grid$topo == "rectangular")
      return(as.matrix(stats::dist(grid$pts, method="maximum")))
  }

  np <- nrow(grid$pts)
  maxdiffx <- grid$xdim/2
  maxdiffy <- max(grid$pts[,2])/2
  
  result <- matrix(0, np, np)
  for (i in 1:(np-1)) {
    for (j in (i+1):np) {
      diffs <- abs(grid$pts[j,] - grid$pts[i,])
      if (diffs[1] > maxdiffx)
        diffs[1] <- 2*maxdiffx - diffs[1]
      if (diffs[2] > maxdiffy)
        diffs[2] <- 2*maxdiffy - diffs[2]

      result[i,j] <- sum(diffs^2)
    }
  }

  sqrt(result + t(result))
}
