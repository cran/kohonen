### Parameter 'main' is mentioned explicitly since in most cases the
### default leaves open ugly space above the plot. We explicitly put
### it 1.2 units above the top row using 'text', if this is within the
### par("usr") range. Else, we use the standard 'title' command.

### Addition 04/11/07: keepMargins, codeRendering and whatmap arguments.
### Decisions to take: * what to do with the "prediction" plot (not
###                      for supersom). Let's eliminate it for the moment.
###                    * what to do with the whatmap argument for
###                      non-supersom maps (simply ignore?)
###                    * Do we need a separate supersom class?

### Adapted for version 2.0: April 11, 2007
### Checked: nope

### Added quality plot, August 30 2007.

### Added default titles, August 31 2007.

"plot.kohonen" <- function (x,
                            type = c("codes", "changes", "counts",
                              "mapping", "property", "quality"),
                            classif = NULL,
                            labels = NULL, pchs = NULL, main = NULL,
                            palette.name = heat.colors, ncolors,
                            bgcol=NULL, zlim = NULL, heatkey = TRUE,
                            property, contin, whatmap = NULL,
                            codeRendering = NULL, keepMargins = FALSE,
                            ...)
{
  type <- match.arg(type)

  switch(type,
         mapping = plot.kohmapping(x, classif, main,
           labels, pchs, bgcol, keepMargins, ...),
         property = plot.kohprop(x, property, main,
           palette.name, ncolors, zlim, heatkey, contin, keepMargins, ...),
         codes = plot.kohcodes(x, main, bgcol, whatmap, codeRendering,
           keepMargins, ...),
         quality = plot.kohquality(x, classif, main,
           palette.name, ncolors, zlim, heatkey,keepMargins, ...),
         counts = plot.kohcounts(x, classif, main,
           palette.name, ncolors, zlim, heatkey, keepMargins, ...),
         changes = plot.kohchanges(x, main, keepMargins, ...))

  invisible()
}


### Overwrite the original plot.somgrid in the class library since
### that leaves open an ugly space at the top of the plot in case of
### hexagonal grids

### Unchanged in version 2.0

plot.somgrid <- function(x, xlim, ylim, ...)
{
  ## Following two lines leave equal amounts of space on both
  ## sides of the plot if no xlim or ylim are given
  if (missing(xlim)) xlim <- c(0, max(x$pts[,1]) + min(x$pts[,1]))
  if (missing(ylim)) ylim <-  c(max(x$pts[,2]) + min(x$pts[,2]), 0)
  MASS::eqscplot(xlim, ylim, axes = FALSE,
                 type = "n", xlab = "", ylab = "", ...)

  invisible()
}


### Adapted for version 2.0: April 11.
### Checked: nope

plot.kohmapping <- function(x, classif, main, labels, pchs, bgcol,
                            keepMargins, ...)
{
  if (is.null(main)) main <- "Mapping plot"
  
  margins <- rep(0.6, 4)
  if (main != "") margins[3] <- margins[3] + 2
  if (!keepMargins) {
    opar <- par("mar")
    on.exit(par(mar = opar))
  }
  par(mar=margins)
    
  if (is.null(classif) & !is.null(x$unit.classif)) {
    classif <- x$unit.classif
  } else {
    if (!is.null(classif$unit.classif))
      classif <- classif$unit.classif
  }
  if (is.null(classif))
    stop("No mapping available")
  
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

  if (is.null(bgcol)) bgcol <- "transparent"
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


### Adapted for version 2.0: April 11.
### Checked: nope

plot.kohprop <- function(x, property, main, palette.name, ncolors,
                         zlim, heatkey, contin, keepMargins, ...)
{
  if (is.null(main)) main <- "Property plot"
  
  margins <- rep(0.6, 4)
  if (heatkey) margins[2] <- margins[2] + 4
  if (main != "") margins[3] <- margins[3] + 2
  if (!keepMargins) {
    opar <- par("mar")
    on.exit(par(mar = opar))
  }
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
    ncolors <- min(length(unique(property[!is.na(property)])), 20)
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
  
  if (heatkey) {
    if (length(unique(property)) < 10 & !contin) {
      plot.heatkey(x, zlim, bgcol, labels = levels(as.factor(property)),
                   contin = contin, ...)
    } else {
      plot.heatkey(x, zlim, bgcol, labels = NULL, contin = contin, ...)
    }
  }
}


### Adapted for version 2.0, April 11.
### Whatmap argument left out because of trouble with indexing:
### changes is a matrix with a number of columns that is equal to the
### layers in whatmap; so there is no safe way to distinguish between,
### say, a whatmap of c(1,2,4) and 2:4 in a six-layer map where only
### the first four layers are used. Anyway.
### Checked: April 13.

plot.kohchanges <- function(x, main, keepMargins, ...)
{
  if (is.null(main)) main <- "Training progress"
  
  nmaps <- ncol(x$changes)
  
  ## check whether a legend is necessary and what names should be used
  if (nmaps > 1) {
    if (!is.null(colnames(x$changes))) {
      varnames <- colnames(x$changes)
    } else {
      varnames <- paste("Matrix", 1:ncol(x$changes))
    }
  }

  ## prepare a second y-axis in case of two maps
  if (nmaps == 2) {
    if (!keepMargins) {
      opar <- par("mar")
      on.exit(par(mar = opar))
    }
    par(mar=c(5.1, 4.1, 4.1, 4.1)) # axis scale to the right as well

    ## scale so that both have the same max value; assume only
    ## positive values.
    huhn <- x$changes
    huhn[,2] <- max(x$changes[,1]) * huhn[,2] / max(x$changes[,2])
    ticks <- pretty(x$changes[,2], length(axTicks(2)))
  } else {
    huhn <- x$changes
  }

  ## plot the plot!
  matplot(huhn, type = "l", lty = 1, main = main, 
          ylab = "Mean distance to closest unit", xlab = "Iteration", ...)
  abline(h=0, col="gray")

  ## plot the second axis
  if (nmaps == 2)
    axis(4, col.axis=2, at=ticks * max(x$changes[,1]) / max(x$changes[,2]),
         labels=ticks)

  ## plot the legend
  if (nmaps > 1)
    legend("topright", legend = varnames, lty=1, col = 1:nmaps, bty="n") 
}


### Adapted for version 2.0: April 11.
### Checked: April 13.

plot.kohcounts <- function(x, classif, main, palette.name, ncolors,
                           zlim, heatkey, keepMargins, ...)
{
  if (is.null(main)) main <- "Counts plot"
  
  if (is.null(classif) & !is.null(x$unit.classif)) {
    classif <- x$unit.classif
  } else {
    if (!is.null(classif$unit.classif))
      classif <- classif$unit.classif
  }
  if (is.null(classif))
    stop("No mapping available")

  counts <- rep(NA, nrow(x$grid$pts))
  huhn <- table(classif)
  counts[as.integer(names(huhn))] <- huhn
  
  plot.kohprop(x, property = counts, main = main,
               palette.name = palette.name, ncolors = ncolors,
               zlim = zlim, heatkey = heatkey, contin = FALSE,
               keepMargins = keepMargins, ...)
}

### Newly written for version 2.0, August 30 2007.
### Revised as a property plot: August 31 2007.

plot.kohquality <- function(x, classif, main, palette.name, ncolors,
                            zlim, heatkey, keepMargins, ...)
{
  if (is.null(main)) main <- "Similarity plot"

  distances <- NULL
  if (is.null(classif) & !is.null(x$unit.classif)) {
    classif <- x$unit.classif
    distances <- x$distances
  } else {
    if (!is.null(classif$unit.classif) & !is.null(classif$distances)) {
      classif <- classif$unit.classif
      distances <- classif$distances
    }
  }
  if (is.null(distances))
    stop("No mapping or mapping distances available")

  similarities <- rep(NA, nrow(x$grid$pts))
  hits <- as.integer(names(table(classif)))
  similarities[hits] <- sapply(split(distances, classif), mean)
                      
  plot.kohprop(x, property = similarities, main = main,
               palette.name = palette.name, ncolors = ncolors,
               zlim = zlim, heatkey = heatkey, contin = TRUE,
               keepMargins = keepMargins, ...)
}


### Adapted for version 2.0: April 11.
### Checked: April 13.
### New elements: whatmap, codeRendering, keepMargins, legend

plot.kohcodes <- function(x, main, bgcol, whatmap, codeRendering,
                          keepMargins, ...)
{
  if (!keepMargins) {
    opar <- par(c("mar", "ask"))
    on.exit(par(opar))
  }

  whatmap <- check.whatmap(x, whatmap)
  nmaps <- length(whatmap)
  
  ## check if x$codes is a list; if so, call this function for every
  ## list element separately.
  if (is.list(x$codes)) {
    if (prod(par("mfrow")) < nmaps) par(ask = TRUE)

    for (i in 1:nmaps) {
      ## make a new object that only has one set of codebook vectors
      huhn <- x
      huhn$codes <- huhn$codes[[whatmap[i]]]

      ## allow a different title for every plot
      if (length(main) == length(x$codes)) {
        main.title <- main[whatmap[i]]
      } else {
        if (length(main) == nmaps) {
          main.title <- main[i]
        } else {
          if (length(main) == 1) {
            main.title <- main
          } else {
            if (is.null(main)) {
              if (!is.null(names(x$codes))) {
                main.title <- names(x$codes)[whatmap[i]]
              } else {
                main.title <- "Codes plot"
              }
            }
          }
        }
      }

      ## allow a different codeRendering for every plot
      if (length(codeRendering) == length(x$codes)) {
        cR <- codeRendering[whatmap[i]]
      } else {
        if (length(codeRendering) == nmaps) {
          cR <- codeRendering[i]
        } else {
          cR <- codeRendering
        }
      }

      plot.kohcodes(huhn, main = main.title, bgcol=bgcol, whatmap = NULL,
                    codeRendering = cR, keepMargins = TRUE, ...)
    }
  } else {
    codes <- x$codes
    nvars <- ncol(codes)
    
    if (is.null(codeRendering)) { ## defaults
      if (nvars < 15) {
        codeRendering <- "segments"
        maxlegendcols <- 3
      } else {
        codeRendering <- "lines"
      }
    }
    
    margins <- rep(0.6, 4)  # no text annotation anywhere
    if (!is.null(main))
      margins[3] <- margins[3] + 2
    par(mar = margins)
    
    if (codeRendering == "segments" & # we need space for the legend here...
        nvars < 15 &
        !is.null(colnames(codes))) {
      plot(x$grid, 
           ylim = c(max(x$grid$pts[,2]) + min(x$grid$pts[,2]), -2))
      current.plot <- par("mfg")
      plot.width <- diff(par("usr")[1:2])

      cex <- 1 # First see if the legend fits
      leg.result <- legend(x = mean(x$grid$pts[,1]), xjust = 0.5,
                           y = 0, yjust = 1,
                           legend = colnames(codes),
                           cex=cex, plot=FALSE,
                           ncol = min(maxlegendcols, nvars),
                           fill = rainbow(nvars))
      while (leg.result$rect$w > plot.width) {
        cex <- cex*0.9 # if too large, decrease text size
        leg.result <- legend(x = mean(x$grid$pts[,1]), xjust = 0.5,
                             y = 0, yjust = 1,
                             legend = colnames(codes),
                             cex=cex, plot=FALSE,
                             ncol = min(maxlegendcols, nvars),
                             fill = rainbow(nvars))
      } # until it fits!

      leg.result <- legend(x = mean(x$grid$pts[,1]), xjust = 0.5,
                           y = 0, yjust = 1, cex=cex,
                           legend = colnames(codes), plot=FALSE,
                           ncol = min(maxlegendcols, nvars),
                           fill = rainbow(nvars), ...)

      par(mfg = current.plot)
      plot(x$grid, 
           ylim = c(max(x$grid$pts[,2]) + min(x$grid$pts[,2]),
             -leg.result$rect$h))

      legend(x = mean(x$grid$pts[,1]), xjust = 0.5,
             y = 0, yjust = 1, cex=cex, plot = TRUE,
             legend = colnames(codes),
             ncol = min(maxlegendcols, nvars),
             fill = rainbow(nvars), ...)
    } else {
      plot(x$grid, ...)
    }
    
    title.y <- max(x$grid$pts[,2]) + 1.2
    if (title.y > par("usr")[4] - .2){
      title(main)
    } else {
      text(mean(range(x$grid$pts[,1])),
           title.y,
           main, adj = .5, cex = par("cex.main"),
           font = par("font.main"))
    }
    
    if (is.null(bgcol)) bgcol <- "transparent"
    symbols(x$grid$pts[, 1], x$grid$pts[, 2],
            circles = rep(0.5, nrow(x$grid$pts)), inches = FALSE,
            add = TRUE, bg = bgcol)
    
    if (codeRendering == "lines") {
      yrange <- range(codes)
      codes <- codes - mean(yrange)
    } else {
      codemins <- apply(codes, 2, min)
      codes <- sweep(codes, 2, codemins)
    }

    switch(codeRendering,
           segments = {
             stars(codes, location = x$grid$pts,
                   labels = NULL, len = 0.4,
                   add=TRUE, col.segments=rainbow(nvars),
                   draw.segments=TRUE)
           },             
           lines = {
             for (i in 1:nrow(x$grid$pts)) { # draw baseline
               if (yrange[1]<0 & yrange[2] > 0) {
                 lines(seq(x$grid$pts[i, 1] - 0.4,
                           x$grid$pts[i, 1] + 0.4,
                           length = 2),
                       rep(x$grid$pts[i, 2], 2),
                       col = "gray")
               }
               lines(seq(x$grid$pts[i, 1] - 0.4,
                         x$grid$pts[i, 1] + 0.4,
                         length = ncol(codes)),
                     x$grid$pts[i, 2] + codes[i, ] * 0.8/diff(yrange),
                     col = "red")
             }
           },
           stars = stars(codes, location = x$grid$pts,
             labels = NULL, len = 0.4, add=TRUE)
           )
    
  }

  invisible()
}


### Unchanged in version 2.0

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
