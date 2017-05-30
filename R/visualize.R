# visualize.R -- visualization functions for Ckmeans.1d.dp
#
# Joe Song
# Created: Oct 1, 2016
# Modified:
#   Oct 12, 2016. Revised ahist().
#   Oct 16, 2016. Moved ahist() function to a new R file ahist.R
#   May 29, 2017. Added two functions plot.Ckmedian.1d.dp() and
#     plot.Cksegs.1d.dp

plotBIC <-
  function(ck, xlab="Number of clusters k",
           ylab = "BIC/n", type="b",
           sub=paste("n =", length(ck$cluster)),
           main=paste("Bayesian information criterion",
                      "(normalized by sample size)", sep="\n"),
           ...)
  {
    if(length(ck$BIC)>0) {
      n <- length(ck$cluster)
      k <- as.numeric(sub("k=", "", names(ck$BIC)))
      graphics::plot(k, ck$BIC/n, xlab=xlab, ylab=ylab,
                     type=type, sub=sub, main=main, ...)

      if(length(ck$BIC)>1) {
        kstar <- length(ck$size)
        graphics::abline(v=kstar, lty="dashed")
        if(kstar <= mean(k)) {
          graphics::text(kstar, min(ck$BIC)/n,
                         paste0(" k*=", kstar), font=2, adj=0)
        } else {
          graphics::text(kstar, min(ck$BIC)/n,
                         paste0("k*=", kstar, " "), font=2, adj=1)
        }
      }
    }
    invisible(ck)
  }

plot.Ckmeans.1d.dp <-
  function(x, xlab=NULL, ylab=NULL, main=NULL,
           sub=NULL, col.clusters=NULL, ...)
    # function(ck, xlab=ck$xname,
    #          ylab=ifelse(ck$yname=="1", "Weight", ck$yname),
    #          main=paste("Optimal k-means clustering of", ck$xname),
    #          sub=paste("n =", length(ck$cluster)),
    #          col.clusters=seq_along(ck$size),
    #          ...)
  {
    ck <- x
    if(is.null(xlab)) xlab <- ck$xname
    if(is.null(ylab)) ylab <- ifelse(ck$yname=="1", "Weight", ck$yname)
    if(is.null(main)) main <- paste("Optimal",
                                    ifelse(ck$yname=="1", "univariate",
                                           "weighted univariate"),
                                    "k-means clustering of", ck$xname)
    if(is.null(sub)) sub=paste("n =", length(ck$cluster))
    if(is.null(col.clusters)) col.clusters <- seq_along(ck$size)

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


plot.Ckmedian.1d.dp <-function(x, xlab=NULL, ylab=NULL, main=NULL,
                               sub=NULL, col.clusters=NULL, ...)
{
  ck <- x
  if(is.null(main)) {
    main <- paste("Optimal",
                  ifelse(ck$yname=="1", "univariate",
                         "weighted univariate"),
                  "k-median clustering of", ck$xname)
  }

  plot.Ckmeans.1d.dp(ck, xlab, ylab, main=main, sub, col.clusters, ...)
}

plot.Cksegs.1d.dp <- function(x, xlab=NULL, ylab=NULL, main=NULL,
                              sub=NULL, col.clusters=NULL, ...)
{
  ck <- x
  if(is.null(xlab)) xlab <- ck$xname
  if(is.null(ylab)) ylab <- ifelse(ck$yname=="1", "Weight", ck$yname)
  if(is.null(main)) main <- paste("Optimal segmentation of", ck$xname, "by", ck$yname,
                                  "\n(Piece-wise constant approximation)")
  if(is.null(sub)) sub=paste("n =", length(ck$cluster))
  if(is.null(col.clusters)) col.clusters <- seq_along(ck$size)

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
    graphics::segments(min(x[ck$cluster == q]), mean(y[ck$cluster == q]),
                       max(x[ck$cluster == q]), mean(y[ck$cluster == q]),
                       col=col.clusters[q], ...) } )

  invisible(ck)
}
