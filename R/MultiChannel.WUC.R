MultiChannel.WUC <- function(x, y,  k=c(1,9)){
  result <- MCW_main(x = x, y = y, Kmin = min(k), Kmax = max(k),
                     estimate_k = "BIC", method = "linear")

  if(length(result$cluster) > 0) {
    k.opt <- max(result$cluster)
  } else {
    k.opt <- 0
  }

  if(nrow(y) == length(x) && all(colSums(y) != 0)) {
    totss <- sum(apply(y, 2, function(y.tmp){
      return(sum(y.tmp * (x - sum(as.numeric(x * y.tmp)) / sum(y.tmp))^2))
    }))
  } else {
    # totss <- sum(scale(x, scale=FALSE)^2) scale function is VERY SLOW!
    totss <- sum((x - sum(as.numeric(x)) / length(x))^2)
  }

  tot.withinss <- sum(result$withinss[1:k.opt])
  betweenss <- totss - tot.withinss
  BIC <- result$BIC
  names(BIC) <- paste0("k=", k)

  r <- structure(
    list(
      cluster = result$cluster, centers = result$centers[1:k.opt],
      withinss = result$withinss[1:k.opt], size = result$size[1:k.opt],
      totss = totss, tot.withinss = tot.withinss, betweenss = betweenss,
      BIC = BIC, xname = deparse(substitute(x)), yname = deparse(substitute(y))
    ),
    class = "MultiChannelClusters")

  return( r )
}

plot.MultiChannelClusters <- function(
  x, xlab=NULL, ylab=NULL, main=NULL,
  sub=NULL, col.clusters=NULL, ...)
{
  ck <- x
  if(is.null(xlab)) xlab <- ck$xname
  if(is.null(ylab)) ylab <- ifelse(ck$yname=="1", "Weight", ck$yname)
  if(is.null(main)) main <-
    "Optimal multichannel weighted clustering"
  # if(is.null(sub)) sub=paste("n =", length(ck$cluster))
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

  graphics::split.screen(c(ncol(y), 1))

  for(ch in seq_len(ncol(y))) {
    graphics::screen(ch)
    main.text <- paste0(
      ifelse(ch == 1, main, ""),
      "\nChannel ", ch, ",  weight = ", colSums(y)[ch])

    graphics::plot(
      x, y[, ch], type="h", xlab=xlab, ylab=ylab,
      main=main.text, # sub=sub,
      col=col.clusters[ck$cluster],
      ...)
  }
  graphics::close.screen(all.screens = TRUE)
  invisible(ck)
}
