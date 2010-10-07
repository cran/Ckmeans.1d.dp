\name{Ckmeans.1d.dp}
\alias{Ckmeans.1d.dp}
\title{Ckmeans.1d.dp: Optimal k-means clustering in one-dimension by dynamic programming}
\description{
Perform optimal distance-based clustering on one-dimensional data.
}
\usage{
Ckmeans.1d.dp(x, k)
}
\arguments{
  \item{x}{ a one-dimensional array containing input data to be clustered}
  \item{k}{ the number of clusters}
}
\author{
	Joe Song <joemsong@cs.nmsu.edu> 
	and 
	Haizhou Wang <hwang@cs.nmsu.edu>
}
\details{
	Distance-based clustering is defined as assigning elements of the input data into k clusters, so that the sum of squares of within-cluster distances (withinss) from each element to its corresponding cluster mean is minimized.  The Ckmeans.1d.dp algorithm groups 1-D data given by x into k clusters.  By dynamic programming, Ckmeans.1.dp guarantees the optimality of clustering -- the sum of withinss for each cluster is always the minimum.  In contrast, the results of the heuristic k-means algorithm may be inconsistent or non-optimal from run to run.  
}
\value{
	An object of class "Ckmeans.1d.dp" which has a print method and is a list with components:
	\item{cluster}{ A vector of integers (from 1:k) indicating the cluster to which each point is allocated.}
	\item{centers}{ A vector of cluster centres.}
	\item{withinss}{ The within-cluster sum of squares for each cluster.}
	\item{size}{ The number of points in each cluster.}
}
\references{
	Wang, H. and Song, M. (2010) Ckmeans.1d.dp: optimal k-means clustering in one dimension by dynamic programming. \emph{The R Journal}. (In preparation)
}

\examples{

# a one-dimensional example with two mixed Gaussian model
x <- rnorm(50, mean=1, sd=0.3)
x <- append(x, rnorm(50, sd=0.3) )
result <- Ckmeans.1d.dp(x, 2)
plot(x, col = result$cluster)
abline(h=result$centers, col = 1:2, pch = 8, cex=2)

}