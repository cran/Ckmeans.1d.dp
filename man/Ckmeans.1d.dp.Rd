\name{Ckmeans.1d.dp}
\alias{Ckmeans.1d.dp}
\title{Optimal K-means Clustering in One-dimension by Dynamic Programming}
\description{
Perform optimal \var{k}-means clustering on one-dimensional data.
}
\usage{
Ckmeans.1d.dp(x, k)
}
\arguments{
  \item{x}{ a one-dimensional array containing input data to be clustered.}
  \item{k}{ the number of clusters.}
}
\author{
	Joe Song \email{joemsong@cs.nmsu.edu}
	and 
	Haizhou Wang \email{hwang@cs.nmsu.edu}
}
\details{
	Distance-based \var{k}-means clustering assigns all elements in the input data in vector \code{x} into \code{k} clusters to minimize the sum of squares of within-cluster distances (\var{withinss}) from each element to its corresponding cluster centre (mean).  The Ckmeans.1d.dp algorithm groups 1-D data given by \code{x} into \code{k} cluster by dynamic programming (Wang and Song, 2011).  It guarantees the optimality of clustering -- the sum of \var{withinss} for each cluster is always the minimum.  In contrast, heuristic \var{k}-means algorithms may be inconsistent or non-optimal from run to run.  
}
\value{
	An object of class "\code{Ckmeans.1d.dp}" which has a print method and is a list with components:
	\item{cluster}{ A vector of integers (\code{1:k}) indicating the cluster to which each point is allocated.}
	\item{centers}{ A vector of cluster centres.}
	\item{withinss}{ The within-cluster sum of squares for each cluster.}
	\item{size}{ The number of points in each cluster.}
}
\references{
	Wang, H. and Song, M. (2011) Ckmeans.1d.dp: optimal \var{k}-means clustering in one dimension by dynamic programming. \emph{The R Journal} \bold{3}(2), 29--33. Retrieved from \url{http://journal.r-project.org/archive/2011-2/RJournal_2011-2_Wang+Song.pdf}
}

\examples{

# Example: clustering data generated from a Gaussian mixture model of two components
x <- rnorm(50, mean=1, sd=0.3)
x <- append(x, rnorm(50, sd=0.3) )
result <- Ckmeans.1d.dp(x, 2)
plot(x, col = result$cluster)
abline(h=result$centers, col = 1:2, pch = 8, cex=2)

}