# ahist.R -- Adaptive histograms using univariate k-means
#
# Joe Song
# Created: Oct 16, 2016. Extracted from visualize.R
# Modified:
#   Oct 12, 2016. Revised ahist().
#   December 6, 2016.
#     1. Added weighted option to ahist() function.
#     2. Broke ahist() to multiple functions

discontinuous.breaks <- function(xs, size)
{ # xs: sorted data
  # size: a vector for the number of points in each cluster

  k <- length(size)

  breaks <- vector("double", length=2*k)

  i <- 1
  for(q in 1:k) {

    # estimate the lower and upper boundaries
    #   cluster q

    j <- i+size[q]-1

    if(size[q] == 1 || stats::var(xs[i:j]) == 0) {
      if(k == 1) {
        # there exists a single cluster of identical elements
        breaks <- 1
        break
      } else {
        # distance to the left (lower) cluster:
        dl <- ifelse(q==1, Inf, xs[i] - xs[i-1])
        # distance to the right (upper) cluster:
        du <- ifelse(q==k, Inf, xs[j+1] - xs[j])
        d <- min(c(dl, du)) / 6
        lb <- xs[i] - d
        ub <- xs[j] + d
      }
    } else {
      # Max difference between consecutive points
      #   within cluster q:
      delta <- max(diff(xs[i:j]))

      lb <- xs[i] - delta
      lmid <- ifelse(q == 1, -Inf, (xs[i-1] + xs[i]) / 2.0)
      lb <- max(c(lb, lmid))

      ub <- xs[j] + delta
      umid <- ifelse(q == k, Inf, (xs[j] + xs[j+1]) / 2.0)
      ub <- min(c(ub, umid))
    }

    breaks[2*q-1] <- lb
    breaks[2*q] <- ub

    i <- i + size[q]
  }

  breaks <- unique(breaks)

  return(breaks)
}

midpoint.breaks <- function(xs, size)
{
  # xs: sorted data
  # size: a vector for the number of points in each cluster

  k <- length(size)

  breaks <- vector("double", length=k+1)

  if(k > 1) {
    i <- size[1]
    for(q in 2:k) {
      breaks[q] <- (xs[i] + xs[i+1])/2
      i <- i + size[q]
    }
  }

  breaks[1] <- xs[1]
  breaks[k+1] <- xs[length(xs)]
  return(breaks)
}

plot.ahist <- function(
  h, x, weight, xlab, wlab, main, col, border, lwd,
  col.stick, lwd.stick, add.sticks,
  skip.empty.bin.color, ...)
{
  opar <- graphics::par(lwd=lwd, mar=c(5,4,4,4)+0.1)
  if(skip.empty.bin.color
     && !is.null(col) && length(col) > 1
     && min(h$counts) == 0) {
    # adjust color for empty bins
    colors <- vector(length = length(h$counts))
    j <- 1
    for(i in seq_along(colors)) {
      colors[i] <- col[(j-1) %% length(col) + 1]
      if(h$counts[i] > 0) {
        j <- j + 1
      }
    }
  } else {
    colors <- col
  }

  graphics::plot(h, main=main, xlab=xlab, col=colors,
                 border=border, ...)

  if(add.sticks) {
    if(length(weight) == 1) {
      # equal weights
      if(h$equidist) {
        graphics::segments(x, -max(h$counts)/10, x, 0,
                           col=col.stick, lwd=lwd.stick)
      } else {
        graphics::segments(x, -max(h$density)/10, x, 0,
                           col=col.stick, lwd=lwd.stick)
      }
    } else {
      # unequal weights

      graphics::mtext(wlab, side = 4, line = 2.75,
                      cex=graphics::par("cex") * graphics::par("cex.lab"), ...)

      graphics::par(new=T)

      graphics::plot(x, weight, type="h", axes=FALSE,
                     xlab=NA, ylab=NA,
                     col=col.stick, lwd=lwd.stick)
      graphics::axis(4)

      mar <- graphics::par("mar")
      mar.legend <- c(mar[1], mar[2]-4, mar[3]-2, mar[4]-2)
      mar.legend <- ifelse(mar.legend >= 0, mar.legend, 0)
      graphics::par(new=T, mar=mar.legend)
      graphics::plot.new()

      ylab <- ifelse(h$equidist, "Frequency", "Density")

      if(0) {
        graphics::legend("bottomleft", c(ylab, wlab), horiz = TRUE,
                         pch=c(22, NA), lty=c(NA, 1),
                         bg="transparent",
                         col=c(border, col.stick),
                         pt.bg=c(colors[1], NA),
                         lwd=c(NA, lwd.stick)
        )
      } else {
        graphics::legend("topleft", c(ylab), bty="n",
                         pch=c(22), lty=c(NA),
                         bg="transparent",
                         col=c(border),
                         pt.bg=c(colors[1]), pt.cex=2, pt.lwd = lwd,
                         lwd=c(NA)
        )
        graphics::legend("topright", c(wlab), bty="n",
                         pch=c(NA), lty=c(1),
                         bg="transparent",
                         col=c(col.stick),
                         pt.bg=c(NA),
                         lwd=c(ifelse(lwd.stick >= 1, lwd.stick, 1))
        )

      }
    }
  }

  graphics::par(opar)
  invisible(h)
}

ahist <- function(
  x, k = c(1,9), breaks=NULL, data=NULL, weight=1,
  plot = TRUE, xlab = deparse(substitute(x)),
  wlab = deparse(substitute(weight)),
  main = NULL, col = "lavender", border = graphics::par("fg"),
  lwd = graphics::par("lwd"),
  col.stick = "gray", lwd.stick = 1, add.sticks=TRUE,
  style = c("discontinuous", "midpoints"),
  skip.empty.bin.color=TRUE,
  ...)
  # adaptive histogram
{
  # x or data may be unsorted, which the return value of the function
  #   must be based on

  style <- match.arg(style)

  xname <- deparse(substitute(x))

  if(is.null(main)) {
    main <- paste(ifelse(is.null(breaks),
                         "Adaptive histogram", "Histogram"),
                  "of", xname)
  }

  if(mode(x) == "numeric") {

    if(length(x) == 0) {
      warning(paste0("\'", xname, "\'", " is empty!"))
      return()
    }

    o <- order(x)
    xs <- x[o]     # xs <- sort(x)

    if(length(weight) > 1) {
      r <- Ckmeans.1d.dp(xs, k=k, y=weight[o])
    } else {
      r <- Ckmeans.1d.dp(xs, k=k)
    }

  } else if(class(x) == "Ckmeans.1d.dp") {
    r <- x
    xname <- r$xname
    x <- data
    if(is.null(data)) {
      stop("data argument must be specified when the class of x is Ckmeans.1d.dp!\n")
    }

    o <- order(data)
    xs <- data[o] # xs <- sort(data)
  } else {
    warning(paste0("\'", xname, "\'", " must be numeric or of class Ckmeans.1d.dp!"))
    return()
  }

  if(is.null(breaks)) {

    if(length(weight) > 1) { # r$size is weighted counts of each cluster.
      # compute actual count of points in each cluster:
      counts <- table(r$cluster[o])
    } else { # in case the cluster is not numbered increasingly with xs
      counts <- r$size
    }

    if(style=="midpoints") {

      breaks <- midpoint.breaks(xs, counts)

    } else if(style=="discontinuous") {

      breaks <- discontinuous.breaks(xs, counts)

    }
  }

  h <- graphics::hist(x, breaks=breaks, plot=FALSE,
                      warn.unused=FALSE, ...)

  if(length(weight) > 1) { # compute weighted count and density

    h$counts <- rep(0, length(breaks)-1)

    table <- stats::aggregate(weight,
                              by=list(bin=cut(x, breaks, labels=FALSE)),
                              FUN=sum) # weighted counts

    h$counts[table[,"bin"]] <- table[,2]

    h$density <- h$counts / diff(breaks) / sum(h$counts) # weighted density
  }

  h$xname <- xname

  if(plot) {
    plot.ahist(h, x, weight, xlab, wlab, main, col, border,
               lwd, col.stick, lwd.stick, add.sticks,
               skip.empty.bin.color, ...)
  } else {
    return(h)
  }
}
