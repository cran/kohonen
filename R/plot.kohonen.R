### Parameter 'main' is mentioned explicitly since in most cases the
### default leaves open ugly space above the plot. We explicitly put
### it 1.2 units above the top row using 'text', if this is within the
### par("usr") range. Else, we use the standard 'title' command.

"plot.kohonen" <- function (x,
                            type = c("changes", "codes", "counts",
                            "mapping", "prediction", "property"),
                            classif, labels = NULL, pchs = NULL,
                            main = NULL,
                            palette.name = heat.colors, ncolors,
                            zlim = NULL, property, heatkey = TRUE,
                            contin, bgcol=NULL, ...)
{
  type <- match.arg(type)

  switch(type,
         mapping = plot.kohmapping(x, classif, main,
           labels, pchs, bgcol, ...),
         prediction = plot.kohpred(x, property, main,
           palette.name, ncolors, zlim, heatkey, labels, contin, ...),
         property = plot.kohprop(x, property, main,
           palette.name, ncolors, zlim, heatkey, contin, ...),
         codes = plot.kohcodes(x, main, ...),
         counts = plot.kohcounts(x, classif, main,
           palette.name, ncolors, zlim, heatkey, ...),
         changes = plot.kohchanges(x, main, ...))

  invisible()
}


### Overwrite the original plot.somgrid in the class library since
### that leaves open an ugly space at the top of the plot in case of
### hexagonal grids

plot.somgrid <- function(x, type = "p", xlim, ylim, ...)
{
  ## Following two lines leave equal amounts of space on both
  ## sides of the plot if no xlim or ylim are given
  if (missing(xlim)) xlim <- c(0, max(x$pts[,1]) + min(x$pts[,1]))
  if (missing(ylim)) ylim <-  c(max(x$pts[,2]) + min(x$pts[,2]), 0)
  MASS::eqscplot(xlim, ylim,
                 axes = FALSE, type = "n", xlab = "", ylab = "", ...)
  invisible()
}


plot.kohpred <- function(x, Y, main, palette.name, ncolors, zlim, heatkey,
                         labels, contin, ...)
{
  if (is.null(x$predict.type))
    stop("Prediction plot only available for supervised SOMs: XYF and BDK\nUse property plot...")
  
  if (missing(Y)) {
    if (is.null(x$codeYs))
      stop("No predictions available")

    if (x$predict.type == "continuous") {
      Y <- x$codeYs
      
      if (ncol(Y) > 1) {       # make separate plots for each variable, use
                               # drop=FALSE to keep column names
        for (i in 1:ncol(Y))
          plot.kohpred(x, Y[,i,drop=FALSE], main, palette.name, ncolors,
                       zlim, heatkey, labels, contin, ...)

        return()
      }
    } else { # classification
      Y <- classmat2classvec(x$codeYs)
    }
  }

  if (missing(contin))
    contin <- (x$predict.type == "continuous")
  if (!contin) Y <- as.factor(Y)
  
  if (missing(ncolors))
    ncolors <- ifelse(contin, 20, min(nlevels(factor(Y)), 20))
  bgcol <- palette.name(ncolors)

  margins <- rep(0.6, 4)
  if (heatkey) 
    margins[2] <- margins[2] + 4
  if (is.null(main))
    main <- colnames(Y)
  
  margins[3] <- margins[3] + 2

  par(mar = margins)

  plot(x$grid, ...)
  title.y <- max(x$grid$pts[,2]) + 1.2
  if (title.y > par("usr")[4] - .2){
    title(main)
  } else {
    text(mean(range(x$grid$pts[,1])),
         title.y,
         main, adj = .5, cex = par("cex.main"),
         font = par("font.main"))
  }
  
  if (is.null(zlim)) {
    if (contin) {
      zlim <- range(Y, finite = TRUE)
    } else {
      zlim <- c(1, nlevels(Y))
    }
  }

  symbols(x$grid$pts[, 1], x$grid$pts[, 2],
          circles = rep(0.5, nrow(x$grid$pts)), inches = FALSE,
          add = TRUE, fg = "black",
          bg = bgcol[as.integer(cut(as.numeric(Y), seq(zlim[1], zlim[2],
            length = ncolors + 1), include.lowest = TRUE))])
  
  if (is.null(labels))
    labels <- levels(Y)
  
  if (heatkey) {
    plot.heatkey(x, zlim, bgcol, labels, contin = contin, ...)
  }
}


plot.kohmapping <- function(x, classif, main, labels, pchs, bgcol, ...)
{
  ifelse(is.null(main),
         par(mar = c(0.6, 0.6, 0.6, 0.6)),
         par(mar = c(0.6, 0.6, 2.6, 0.6)))
    
  if (missing(classif) & !is.null(x$unit.classif)) {
    classif <- x$unit.classif
  } else {
    if (!is.null(classif$unit.classif))
      classif <- classif$unit.classif
  }
  
  plot(x$grid, ...)
  title.y <- max(x$grid$pts[,2]) + 1.2
  if (title.y > par("usr")[4] - .2){
    title(main)
  } else {
    text(mean(range(x$grid$pts[,1])),
         title.y,
         main, adj = .5, cex = par("cex.main"),
         font = par("font.main"))
  }

  if (is.null(bgcol)) bgcol <- "gray"
  symbols(x$grid$pts[, 1], x$grid$pts[, 2],
          circles = rep(0.5, nrow(x$grid$pts)),
          inches = FALSE, add = TRUE, bg = bgcol)
  if (is.null(labels) & !is.null(pchs))
    points(x$grid$pts[classif, 1] + rnorm(length(classif), 0, 0.12),
           x$grid$pts[classif, 2] + rnorm(length(classif), 0, 0.12),
           pch = pchs, ...)
  if (!is.null(labels))
    text(x$grid$pts[classif, 1] + rnorm(length(classif), 0, 0.12),
         x$grid$pts[classif, 2] + rnorm(length(classif), 0, 0.12),
         labels, ...)
}


plot.kohprop <- function(x, property, main, palette.name, ncolors,
                         zlim, heatkey, contin, ...)
{
  margins <- rep(0.6, 4)
  if (heatkey) 
    margins[2] <- margins[2] + 4
  if (!is.null(main)) 
    margins[3] <- margins[3] + 2
  par(mar = margins)
  
  plot(x$grid, ...)
  title.y <- max(x$grid$pts[,2]) + 1.2
  if (title.y > par("usr")[4] - .2){
    title(main)
  } else {
    text(mean(range(x$grid$pts[,1])),
         title.y,
         main, adj = .5, cex = par("cex.main"),
         font = par("font.main"))
  }
  
  if (is.null(zlim))
    zlim <- range(property, finite = TRUE)

  if (missing(ncolors)) 
    ncolors <- min(length(unique(property)), 20)
  bgcol <- palette.name(ncolors)

  bgcolors <- rep("gray", nrow(x$grid$pts))
  showcolors <- as.integer(cut(property,
                               seq(zlim[1], zlim[2],
                                   length = ncolors + 1),
                               include.lowest = TRUE))
  bgcolors[!is.na(showcolors)] <- bgcol[showcolors[!is.na(showcolors)]]

  symbols(x$grid$pts[, 1], x$grid$pts[, 2],
          circles = rep(0.5, nrow(x$grid$pts)), inches = FALSE,
          add = TRUE, fg = "black", bg = bgcolors)

  ## if contin, a pretty labelling of z colors will be used; if not,
  ## all colours will have their own label. The latter only if the
  ## number of categories is smaller than 10, unless explicitly
  ## given.
  if (missing(contin))
    contin <- !(length(unique(property)) < 10)
  
  if (heatkey)
    plot.heatkey(x, zlim, bgcol, labels = NULL, contin = contin, ...)
}


plot.kohchanges <- function(x, main, ...)
{
  if (is.matrix(x$changes)) { # for supervised networks
    par(mar=c(5.1, 4.1, 4.1, 4.1)) # axis scale to the right as well
    
    ## scale so that both have the same max value; assume only
    ## positive values.
    huhn <- x$changes
    huhn[,2] <- max(x$changes[,1]) * huhn[,2] / max(x$changes[,2])
    ticks <- pretty(x$changes[,2], length(axTicks(2)))
    
    matplot(huhn, type = "l", lty = 1, col=c(1,2), main = main, 
            ylab = "Mean similarity change", xlab = "Iteration", ...)
    axis(4, col.axis=2, at=ticks * max(x$changes[,1]) / max(x$changes[,2]),
         labels=ticks)
    legend("topright", legend = c("X update", "Y update"),
           lty=c(1,1), col=c(1,2), bty="n") 
  } else {
    plot(x$changes, type = "l", ylab = "Mean change", main = main,
         xlab = "Iteration", ...)
  }
}


plot.kohcounts <- function(x, classif, main, palette.name, ncolors,
                           zlim, heatkey, ...) {
  margins <- rep(0.6, 4)
  if (heatkey) 
    margins[2] <- margins[2] + 4
  if (!is.null(main)) 
    margins[3] <- margins[3] + 2
  par(mar = margins)

  if (missing(classif) & !is.null(x$unit.classif)) {
    classif <- x$unit.classif
  } else {
    if (!is.null(classif$unit.classif))
      classif <- classif$unit.classif
  }
  
  plot(x$grid, ...)
  title.y <- max(x$grid$pts[,2]) + 1.2
  if (title.y > par("usr")[4] - .2){
    title(main)
  } else {
    text(mean(range(x$grid$pts[,1])),
         title.y,
         main, adj = .5, cex = par("cex.main"),
         font = par("font.main"))
  }
    
  bgcolors <- rep("gray", nrow(x$grid$pts))
  hits <- as.integer(names(table(classif)))
  counts <- as.integer(table(classif))
  if (is.null(zlim))
    zlim <- c(0, max(counts))

  if (missing(ncolors)) ncolors <- min(20, max(counts))
  bgcol <- palette.name(ncolors)
  showcolors <- as.integer(cut(counts,
                               seq(zlim[1], zlim[2],
                                   length = ncolors + 1),
                               include.lowest = TRUE))
  bgcolors[hits] <- bgcol[showcolors]
  
  symbols(x$grid$pts[, 1], x$grid$pts[, 2],
          circles = rep(0.5, nrow(x$grid$pts)), inches = FALSE,
          add = TRUE, fg = "black", bg = bgcolors)
  
  if (heatkey) plot.heatkey(x, zlim, bgcol, contin = TRUE, ...)
}


plot.kohcodes <- function(x, main, ...)
{
  ifelse(is.null(main),
         par(mar = c(0.6, 0.6, 0.6, 0.6)),
         par(mar = c(0.6, 0.6, 2.6, 0.6)))
  
  plot(x$grid, ...)
  title.y <- max(x$grid$pts[,2]) + 1.2
  if (title.y > par("usr")[4] - .2){
    title(main)
  } else {
    text(mean(range(x$grid$pts[,1])),
         title.y,
         main, adj = .5, cex = par("cex.main"),
         font = par("font.main"))
  }
  
  symbols(x$grid$pts[, 1], x$grid$pts[, 2],
          circles = rep(0.5, nrow(x$grid$pts)), inches = FALSE,
          add = TRUE, bg = "white")

  codes <- sweep(x$codes, 2, apply(x$codes, 2, min))
  nvars <- ncol(codes)
  
  if (nvars < 15) {
    stars(codes, location = x$grid$pts, labels = NULL, len = 0.4,
          add=TRUE, col.segments=rainbow(nvars), draw.segments=TRUE)
  } else {
    for (i in 1:nrow(x$grid$pts)) {
      lines(seq(x$grid$pts[i, 1] - 0.4,
                x$grid$pts[i, 1] + 0.4,
                length = ncol(codes)),
            x$grid$pts[i, 2] - 0.2 + codes[i, ] * 0.5/max(codes[i, ]),
            col = "red")
    }
  }
}


plot.heatkey <- function (x, zlim, bgcol, labels, contin, ...)
{
  ncolors <- length(bgcol)
  
  yrange <- range(x$grid$pts[, 2])
  smallestx <- min(x$grid$pts[,1])
  ## A width of .2 looks OK on my screen
  xleft <- c(smallestx - 1.2, smallestx - 1)
  yleft <- seq(yrange[1] - 0.5,
               yrange[2] + 0.5,
               length = ncolors + 1)
  rect(xleft[1], yleft[1:ncolors],
       xleft[2], yleft[2:(ncolors + 1)],
       border = "black", col = bgcol,
       xpd = TRUE)

  cex <- list(...)$cex

  if (contin) {
    zvals <- pretty(zlim)
    zvals <- zvals[zvals <= max(zlim) & zvals >= min(zlim)]
    yvals <- yrange[1] - .5 + (diff(yrange) + 1)*(zvals - zlim[1])/diff(zlim)

    text(xleft[2] - 1.3*diff(xleft),
         yvals,
         formatC(zvals),
         xpd=TRUE, adj=1, cex=cex)
  } else {
    if (is.null(labels))
      labels <- 1:ncolors
    
    text(xleft[2] - 1.3 * diff(xleft),
         yleft[-1] - 0.5*diff(yleft[1:2]),
         sort(labels),
         xpd = TRUE, adj=1, cex=cex)
  }
}
