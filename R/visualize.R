# visualize.R -- visualization functions for Ckmeans.1d.dp
#
# Joe Song
# Created: Oct 1, 2016

ahist <- function(
  x, k = c(1,9), plot = TRUE, xlab = deparse(substitute(x)),
  main = paste("Adaptive histogram of", deparse(substitute(x))),
  col = NULL, lwd = graphics::par("lwd"), col.stick = "gray", lwd.stick = 1,
  ...)
  # adaptive histogram
{
  xs <- sort(x)
  r <- Ckmeans.1d.dp(xs, k=k)
  kopt <- length(r$size)
  breaks <- vector("double", length=kopt+1)
  i <- r$size[1]

  if(kopt > 1) {
    for(q in 2:kopt) {
      breaks[q] <- (xs[i] + xs[i+1])/2
      i <- i + r$size[q]
    }
  }

  breaks[1] <- xs[1]
  breaks[kopt+1] <- xs[length(xs)]

  h <- graphics::hist(x, breaks=breaks, plot=FALSE,
                      warn.unused=FALSE, ...)
  h$xname <- deparse(substitute(x))

  if(plot) {
    opar <- graphics::par(lwd=lwd)
    graphics::plot(h, main=main, xlab=xlab, col=col, ...)
    if(h$equidist) {
      graphics::segments(x, -max(h$count)/10, x, 0,
                         col=col.stick, lwd=lwd.stick)
    } else {
      graphics::segments(x, -max(h$density)/10, x, 0,
                         col=col.stick, lwd=lwd.stick)
    }
    graphics::par(opar)
    invisible(h)
  } else {
    return(h)
  }
}

plot.Ckmeans.1d.dp <-
  function(x, xlab=NULL, ylab=NULL, main=NULL,
           sub=NULL, col.clusters=NULL, ...)
  {
    ck <- x
    if(is.null(xlab)) xlab <- ck$xname
    if(is.null(ylab)) ylab <- ifelse(ck$yname=="1", "Weight", ck$yname)
    if(is.null(main)) main <- paste("Optimal k-means clustering of", ck$xname)
    if(is.null(sub)) sub=paste("n =", length(ck$cluster))
    if(is.null(col.clusters)) col.clusters <- seq_along(x$size)

    if(exists(ck$xname, mode="numeric")) {
      x <- get(ck$xname, mode="numeric")
    } else {
      x <- eval(parse(text=ck$xname))
    }

    if(exists(ck$yname, mode="numeric")) {
      y <- get(ck$yname, mode="numeric")
    } else {
      y <- eval(parse(text=ck$yname))
    }

    if(length(y) == 1) {
      y <- rep(y, length(x))
    }

    graphics::plot(x, y, type="p",
                   xlab=xlab, ylab=ylab, main=main, sub=sub,
                   col=col.clusters[ck$cluster],
                   ...)

    ks <- seq_along(ck$size)
    sapply(ks, function(q) {
      graphics::segments(x[ck$cluster == q], 0,
                         x[ck$cluster == q], y[ck$cluster == q],
                         col=col.clusters[q], ...) } )

    invisible(ck)
  }
