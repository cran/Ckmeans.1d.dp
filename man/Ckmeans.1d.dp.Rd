\name{Ckmeans.1d.dp}
\alias{Ckmeans.1d.dp}
\title{Optimal and Fast Univariate \var{k}-Means Clustering}

\description{
Perform optimal and fast univariate \var{k}-means clustering by dynamic programming and divide-and-conquer.
}

\usage{
Ckmeans.1d.dp(x, k=c(1,9), y=1)
}

\arguments{
  \item{x}{a numeric vector of data to be clustered. All \code{NA} elements must be removed from \code{x} before calling this function.}

  \item{k}{either an exact integer number of clusters, or a vector of length two specifying the minimum and maximum numbers of clusters to be examined. The default is \code{c(1,9)}. When \code{k} is a range, the actual number of clusters is determined by Bayesian information criterion.}

  \item{y}{a value of 1 to specify equal weights or a numeric vector of unequal weights for each element. The default weight is one. It is highly recommended to use positive (instead of zero) weights to account for the influence of every element. The weights have a strong impact on the clustering result.}
}

\author{
	Joe Song and Haizhou Wang
}

\details{
In contrast to the heuristic \var{k}-means algorithms implemented in function \code{\link{kmeans}}, this function optimally assigns elements in numeric vector \code{x} into \code{k} clusters by dynamic programming (Wang and Song, 2011). It minimizes the total of within-cluster sums of squared distances (\var{withinss}) between each element and its corresponding cluster mean. When a range is provided for \code{k}, the exact number of clusters is determined by Bayesian information criterion. Different from the heuristic \var{k}-means algorithms whose results may be non-optimal or change from run to run, the result of Ckmeans.1d.dp is guaranteed to be optimal and reproducible, and its advantage in efficiency and accuracy over heuristic \var{k}-means methods becomes increasingly evident as \var{k} increases.
}

\value{
	An object of class \code{"Ckmeans.1d.dp"} with a \code{print} method. It is a list containing the following components:

	\item{cluster}{a vector of clusters assigned to each element in \code{x}. Each cluster is indexed by an integer from 1 to \code{k}.}

	\item{centers}{a numeric vector of the (weighted) means for each cluster.}

	\item{withinss}{a numeric vector of the (weighted) within-cluster sum of squares for each cluster.}

	\item{size}{a vector of the (weighted) number of elements in each cluster.}

	\item{totss}{total sum of squared distances between each element and the sample mean. This statistic is not dependent on the clustering result.}

	\item{tot.withinss}{total sum of within-cluster squared distances between each element and its cluster mean. This statistic is minimized given the number of clusters.}

	\item{betweenss}{sum of squared distances between each cluster mean and sample mean weighed by cluster size. This statistic is maximized given the number of clusters.}

	\item{xname}{a character string. The actual name of the \code{x} argument.}

	\item{yname}{a character string. The actual name of the \code{y} argument.}

}

\references{
Wang, H. and Song, M. (2011) Ckmeans.1d.dp: optimal \var{k}-means clustering in one dimension by dynamic programming. \emph{The R Journal} \bold{3}(2), 29--33. Retrieved from \url{https://journal.r-project.org/archive/2011-2/RJournal_2011-2_Wang+Song.pdf}
}

\seealso{
  The \code{\link{kmeans}} function in package \pkg{\link{stats}} that implements several heuristic \eqn{k}{k}-means algorithms.
}

\examples{
# Ex. 1 The number of clusters is provided.
# Generate data from a Gaussian mixture model of three components
x <- c(rnorm(50, sd=0.3), rnorm(50, mean=1, sd=0.3), rnorm(50, mean=2, sd=0.3))
# Divide x into 3 clusters
k <- 3
result <- Ckmeans.1d.dp(x, k)
plot(result)
plot(x, col=result$cluster, pch=result$cluster, cex=1.5,
     main="Optimal k-means clustering given k",
     sub=paste("Number of clusters given:", k))
abline(h=result$centers, col=1:k, lty="dashed", lwd=2)
legend("bottom", paste("Cluster", 1:k), col=1:k, pch=1:k, cex=1.5, bty="n")

# Ex. 2 The number of clusters is determined by Bayesian information criterion
# Generate data from a Gaussian mixture model of three components
x <- c(rnorm(50, mean=-1, sd=0.3), rnorm(50, mean=1, sd=1), rnorm(50, mean=2, sd=0.4))
# Divide x into k clusters, k automatically selected (default: 1~9)
result <- Ckmeans.1d.dp(x)
plot(result)
k <- max(result$cluster)
plot(x, col=result$cluster, pch=result$cluster, cex=1.5,
     main="Optimal k-means clustering with k estimated",
     sub=paste("Number of clusters is estimated to be", k))
abline(h=result$centers, col=1:k, lty="dashed", lwd=2)
legend("topleft", paste("Cluster", 1:k), col=1:k, pch=1:k, cex=1.5, bty="n")

# Ex. 3 Segmenting a time course using weighted k-means
n <- 160
t <- seq(0, 2*pi*2, length=n)
n1 <- 1:(n/2)
n2 <- (max(n1)+1):n
y1 <- abs(sin(1.5*t[n1]) + 0.1*rnorm(length(n1)))
y2 <- abs(sin(0.5*t[n2]) + 0.1*rnorm(length(n2)))
y <- c(y1, y2)

w <- y^8 # stress the peaks
res <- Ckmeans.1d.dp(t, k=c(1:10), w)
plot(res)
plot(t, w, main = "Time course segmentation",
     col=res$cluster, pch=res$cluster,
     xlab="Time t", ylab="Transformed intensity w",
     type="h")
abline(v=res$centers, col="chocolate", lty="dashed")
text(res$centers, max(w) * .95, cex=0.5, font=2,
     paste(round(res$size / sum(res$size) * 100, 2), "/100"))

}
