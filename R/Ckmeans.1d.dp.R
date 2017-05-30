## Ckmeans.1d.dp.R
##
##   Ckmeans.1d.dp()-----Interface funcition to call C++ version kmeans.1d.dp()
##
##   Haizhou Wang
##   Computer Science Department
##   New Mexico State University
##   hwang@cs.nmsu.edu
##
##   Created: Oct 20, 2009
##
## Modified:
#    May 17, 2016. MS
#    September 25, 2016. MS. Introduced function ahist()

##Ckmeans.1d.dp : function which implement optimal one-dimensional clustering
## x is one-dimensional input vector
## k indicates cluster level
cluster.1d.dp <- function( x, k, y, method, estimate.k, criterion, xname, yname )
{
  if(is.null(k)) {
    k <- 1: min( 9, length(x) )
  } else {
    k <- as.integer(ceiling(k))
  }

  if(length(k) > 1) {
    k.min <- min(k)
    k.max <- max(k)
  } else {
    k.min <- k
    k.max <- k
  }

  ##Check to see if k is less than 0.
  ##If k is less than 0, stop and exit.
  if( k.max <= 0 ) {
    stop ("Can NOT classify vector into 0 or less cluster\n")
  }

  n.unique <- length(unique(x))
  ##Check to see if cluster level bigger than the unique number of the input vector
  ##If k is bigger than the unique number of the input vector,
  ##force k set to the number of unique number of input.

  if(n.unique == 0) {
    warning(paste("Input vector", deparse(substitute(x)), "is empty!\n"))
  } else {
    if(n.unique < k.min) {
      warning("Min number of clusters is greater than the unique number of elements in\n",
              "the input vector, both k.min and k.max are set to the number of\n",
              "unique number of input values.\n")
      k.min <- n.unique
      k.max <- n.unique
    } else if (n.unique >= k.min && n.unique < k.max) {
      warning("Max number of clusters is greater than the unique number of\n",
              "elements in the input vector, and k.max is set to the number of\n",
              "unique number of input values.\n")
      k.max <- n.unique
    }
  }

  ##Form data which will be passed to external C++ function.
  clusters <- vector("integer", length(x))
  center <- vector("double", k.max)
  withinss <- vector("double", k.max)
  size <- vector("double", k.max) # vector("integer", k.max)
  BIC <- vector("double", k.max-k.min+1)

  #Call external C++ function
  if(criterion == "L2") {
    result <- .C("Ckmeans_1d_dp", PACKAGE="Ckmeans.1d.dp",
                 data=as.double(x), length=as.integer(length(x)),
                 weight=as.double(y), weight_length=as.integer(length(y)),
                 Kmin=as.integer(k.min), Kmax=as.integer(k.max),
                 cluster=as.integer(clusters), centers=as.double(center),
                 withinss=as.double(withinss), size=as.double(size), # size=as.integer(size),
                 BIC=as.double(BIC),
                 estimate.k=as.character(estimate.k),
                 method=as.character(method))
    class <- "Ckmeans.1d.dp"
  } else if (criterion == "L1") {
    result <- .C("Ckmedian_1d_dp", PACKAGE="Ckmeans.1d.dp",
                 data=as.double(x), length=as.integer(length(x)),
                 weight=as.double(y), weight_length=as.integer(length(y)),
                 Kmin=as.integer(k.min), Kmax=as.integer(k.max),
                 cluster=as.integer(clusters), centers=as.double(center),
                 withinss=as.double(withinss), size=as.double(size), # size=as.integer(size),
                 BIC=as.double(BIC),
                 estimate.k=as.character(estimate.k),
                 method=as.character(method))
    class <- "Ckmedian.1d.dp"
  } else if(criterion == "L2Y") {
    result <- .C("Cksegs_1d_dp", PACKAGE="Ckmeans.1d.dp",
                 data=as.double(x), length=as.integer(length(x)),
                 weight=as.double(y), weight_length=as.integer(length(y)),
                 Kmin=as.integer(k.min), Kmax=as.integer(k.max),
                 cluster=as.integer(clusters), centers=as.double(center),
                 withinss=as.double(withinss), size=as.double(size), # size=as.integer(size),
                 BIC=as.double(BIC),
                 estimate.k=as.character(estimate.k),
                 method=as.character(method))
    class <- "Cksegs.1d.dp"
  }

  if(length(result$cluster) > 0) {
    k.opt <- max(result$cluster) # length(unique(result$cluster))
  } else {
    k.opt <- 0
  }

  if(n.unique > 0) {
    if(k.min < k.max) {
      if (k.opt == k.min && k.min != 1) {
        warning("Min number of clusters used. Consider decreasing k!\n")
      } else if(k.opt == k.max && k.max != length(x)) {
        warning("Max number of clusters used. Consider increasing k!\n")
      }
    }
  }

  if(length(y) == length(x) && sum(y) != 0 && criterion != "L2Y") {
    totss <- sum(y * (x - sum(as.numeric(x * y)) / sum(y))^2)
  } else if(criterion == "L2Y") {
    totss <- sum((y - sum(as.numeric(y)) / length(y))^2)
  } else {
    # totss <- sum(scale(x, scale=FALSE)^2) scale function is VERY SLOW!
    totss <- sum((x - sum(as.numeric(x)) / length(x))^2)
  }

  tot.withinss <- sum(result$withinss[1:k.opt])
  betweenss <- totss - tot.withinss
  BIC <- result$BIC
  names(BIC) <- paste0("k=", k.min : k.max)

  r <- structure(
    list(
      cluster = result$cluster, centers = result$centers[1:k.opt],
      withinss = result$withinss[1:k.opt], size = result$size[1:k.opt],
      totss = totss, tot.withinss = tot.withinss, betweenss = betweenss,
      BIC = BIC, xname=xname, yname=yname
    ),
    class = class)

  return( r )
} ##end of cluster.1d.dp()


Ckmeans.1d.dp <- function( x, k=c(1,9), y=1,
                           method=c("linear", "loglinear", "quadratic"),
                           estimate.k=c("BIC", "BIC 3.4.12") )
{
  method <- match.arg(method)
  estimate.k <- match.arg(estimate.k)
  cluster.1d.dp(x, k, y, method, estimate.k, "L2",
                deparse(substitute(x)), deparse(substitute(y)))
}

Ckmedian.1d.dp <- function( x, k=c(1,9), y=1,
                            method=c("linear", "loglinear", "quadratic"),
                            estimate.k=c("BIC", "BIC 3.4.12") )
{
  method <- match.arg(method)
  estimate.k <- match.arg(estimate.k)

  if(length(y) > 1) {
    warning("Ckmedian.1d.dp() ignores weight vector y and uses equal weight instead!")
  }

  cluster.1d.dp(x, k, y=1, method, estimate.k, "L1",
                deparse(substitute(x)), deparse(substitute(y)))
}

Cksegs.1d.dp <- function(y, k=c(1,9), x=seq_along(y),
                         method=c("quadratic", "linear", "loglinear"),
                         estimate.k=c("BIC", "BIC 3.4.12") )
{
  method <- match.arg(method)
  estimate.k <- match.arg(estimate.k)
  cluster.1d.dp(x, k, y, method, estimate.k, "L2Y",
                deparse(substitute(x)), deparse(substitute(y)))
}
