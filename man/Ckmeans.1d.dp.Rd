\name{Ckmeans.1d.dp}
\alias{Ckmeans.1d.dp}
\title{Optimal K-means Clustering in One-dimension by Dynamic Programming}

\description{
Perform optimal \var{k}-means clustering on one-dimensional data.
}

\usage{
Ckmeans.1d.dp(x, k=c(1,9))
}

\arguments{
  \item{x}{a numeric vector of data to be clustered. The vector cannot contain any \code{NA}.}
  \item{k}{either an exact integer number of clusters, or a vector of two integers specifying the minimum and maximum numbers of clusters. The default is \code{c(1,9)}. When a range is provided, the actual number of clusters will be determined by Bayesian information criterion.}
}

\author{
	Joe Song and Haizhou Wang
}

\details{
Distance-based \var{k}-means clustering assigns all elements in the input vector \code{x} into \code{k} clusters to minimize the sum of squares of within-cluster distances (\var{withinss}) from each element to its corresponding cluster center (mean). The algorithm was first described by Richard Bellman (1973). When a range is provided for \code{k}, the exact number of clusters will be determined within the range by Bayesian information criterion.  The Ckmeans.1d.dp algorithm groups 1-D data given by \code{x} into \code{k} cluster by dynamic programming (Wang and Song, 2011).  It guarantees the optimality of clustering -- the sum of \var{withinss} for each cluster is always the minimum.  In contrast, heuristic \var{k}-means algorithms may be inconsistent or non-optimal from run to run.

Since version 3.4.0, the runtime of this program has been greatly reduced with two changes. First, an upper bound for the sum of within cluster square distances is checked to reduce the dynamic programming search space. Second, the unnecessary calculation of \eqn{n-1}{n-1} elements in the dynamic programming matrix that do not influence the final result has been eliminated. This led to enormous reduction in runtime when the number of cluster is 2: assigning one million points into two clusters took half a second on iMac with 2.93 GHz Intel Core i7 processor. Overall, the practical runtime of the program has been greatly reduced. Although the asymptotic runtime on \eqn{n}{n} points is still \eqn{O(\max(k)n^2)}{O(max(k) n^2)} for \eqn{k>2}{k>2}, it has been improved to \eqn{O(n\lg n)}{O(n lg(n))} for \eqn{k=2}{k=2}. These improvements were not discussed in (Wang and Song, 2011).
}

\value{
	An object of class \code{"Ckmeans.1d.dp"} with a \code{print} method. It is a list containing the following components:
	\item{cluster}{a vector of clusters assigned to each element in \code{x}. Each cluster is indexed by an integer from 1 to \code{k}.}
	\item{centers}{a numeric vector of the means for each cluster.}
	\item{withinss}{a numeric vector of the within-cluster sum of squares for each cluster.}
	\item{size}{a vector of the number of points in each cluster.}
}

\references{
Wang, H. and Song, M. (2011) Ckmeans.1d.dp: optimal \var{k}-means clustering in one dimension by dynamic programming. \emph{The R Journal} \bold{3}(2), 29--33. Retrieved from \url{http://journal.r-project.org/archive/2011-2/RJournal_2011-2_Wang+Song.pdf}

Bellman, R. (1973) A note on cluster analysis and dynamic programming. \emph{Mathematical Biosciences} \bold{18}(3), 311--312.
}

\examples{
# Ex. 1 The number of clusters is provided.
# Generate data from a Gaussian mixture model of two components
x <- c(rnorm(50, sd=0.3), rnorm(50, mean=1, sd=0.3))
# Divide x into 2 clusters
k <- 2
result <- Ckmeans.1d.dp(x, k)
plot(x, col=result$cluster, pch=result$cluster, cex=1.5,
     main="Optimal k-means clustering given k",
     sub=paste("Number of clusters given:", k))
abline(h=result$centers, col=1:k, lty="dashed", lwd=2)
legend("bottom", paste("Cluster", 1:k), col=1:k, pch=1:k, cex=1.5, bty="n")

# Ex. 2 The number of clusters is determined by Bayesian information criterion
# Generate data from a Gaussian mixture model of two components
x <- c(rnorm(50, mean=-1, sd=0.3), rnorm(50, mean=1, sd=1))
# Divide x into k clusters, k automatically selected (default: 1~9)
result <- Ckmeans.1d.dp(x)
k <- max(result$cluster)
plot(x, col=result$cluster, pch=result$cluster, cex=1.5,
     main="Optimal k-means clustering with k estimated",
     sub=paste("Number of clusters is estimated to be", k))
abline(h=result$centers, col=1:k, lty="dashed", lwd=2)
legend("topleft", paste("Cluster", 1:k), col=1:k, pch=1:k, cex=1.5, bty="n")

}
