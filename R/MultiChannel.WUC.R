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
    ))

  return( r )
}
