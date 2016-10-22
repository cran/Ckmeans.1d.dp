# ahist.R -- Adaptive histograms using univariate k-means
#
# Joe Song
# Created: Oct 16, 2016. Extracted from visualize.R
# Modified:
#   Oct 12, 2016. Revised ahist().

ahist <- function(
  x, k = c(1,9), breaks=NULL, data=NULL, plot = TRUE,
  xlab = deparse(substitute(x)), main = NULL,
  col = NULL, lwd = graphics::par("lwd"),
  col.stick = "gray", lwd.stick = 1, add.sticks=TRUE,
  style = c("discontinuous", "midpoints"),
  skip.empty.bin.color=TRUE,
  ...)
  # adaptive histogram
{
  style <- match.arg(style)

  if(is.null(main)) {
    main <- paste(ifelse(is.null(breaks),
                         "Adaptive histogram", "Histogram"),
                  "of", deparse(substitute(x)))
  }

  if(mode(x) == "numeric") {
    xname <- deparse(substitute(x))
    xs <- sort(x)
    r <- Ckmeans.1d.dp(xs, k=k)
  } else if(class(x) == "Ckmeans.1d.dp") {
    r <- x
    xname <- r$xname
    x <- data
    if(is.null(data)) {
      cat("ERROR: data must be specified if x is of class Ckmeans.1d.dp!\n")
    }
    xs <- sort(data)
  } else {
    warning("ERROR: unrecognized input!")
    return()
  }

  if(is.null(breaks)) {

    kopt <- length(r$size)

    if(style=="midpoints") {
      breaks <- vector("double", length=kopt+1)

      if(kopt > 1) {
        i <- r$size[1]
        for(q in 2:kopt) {
          breaks[q] <- (xs[i] + xs[i+1])/2
          i <- i + r$size[q]
        }
      }

      breaks[1] <- xs[1]
      breaks[kopt+1] <- xs[length(xs)]

    } else if(style=="discontinuous") {

      breaks <- vector("double", length=2*kopt)

      i <- 1
      for(q in 1:kopt) {

        # estimate the lower and upper boundaries
        #   cluster q

        j <- i+r$size[q]-1

        if(r$size[q] == 1 || stats::var(xs[i:j]) == 0) {
          if(kopt == 1) {
            # there exists a single cluster of identical elements
            breaks <- 1
            break
          } else {
            # distance to the left (lower) cluster:
            dl <- ifelse(q==1, Inf, xs[i] - xs[i-1])
            # distance to the right (upper) cluster:
            du <- ifelse(q==kopt, Inf, xs[j+1] - xs[j])
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
          umid <- ifelse(q == kopt, Inf, (xs[j] + xs[j+1]) / 2.0)
          ub <- min(c(ub, umid))
        }

        breaks[2*q-1] <- lb
        breaks[2*q] <- ub

        i <- i + r$size[q]
      }

      breaks <- unique(breaks)
    }
  }

  h <- graphics::hist(x, breaks=breaks, plot=FALSE,
                      warn.unused=FALSE, ...)
  h$xname <- xname

  if(plot) {
    opar <- graphics::par(lwd=lwd)
    if(skip.empty.bin.color
       && !is.null(col) && length(col) > 1
       && min(h$count) == 0) {
      # adjust color for empty bins
      colors <- vector(length = length(h$count))
      j <- 1
      for(i in seq_along(colors)) {
        colors[i] <- col[(j-1) %% length(col) + 1]
        if(h$count[i] > 0) {
          j <- j + 1
        }
      }
    } else {
      colors <- col
    }

    graphics::plot(h, main=main, xlab=xlab, col=colors, ...)

    if(add.sticks) {
      if(h$equidist) {
        graphics::segments(x, -max(h$count)/10, x, 0,
                           col=col.stick, lwd=lwd.stick)
      } else {
        graphics::segments(x, -max(h$density)/10, x, 0,
                           col=col.stick, lwd=lwd.stick)
      }
    }

    graphics::par(opar)
    invisible(h)
  } else {
    return(h)
  }
}
