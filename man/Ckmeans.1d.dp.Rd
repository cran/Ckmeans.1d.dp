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
  \item{x}{a one-dimensional array containing input data to be clustered.}
  \item{k}{the number of clusters, or an array of required min and max numbers of clusters. The default is \code{c(1,9)}. When a rannge is provided, the number of clusters will be determined within the range by Bayesian information criterion.}
}

\author{
	Joe Song and Haizhou Wang
}

\details{
	Distance-based \var{k}-means clustering assigns all elements in the input vector \code{x} into \code{k} clusters to minimize the sum of squares of within-cluster distances (\var{withinss}) from each element to its corresponding cluster centre (mean).  When a rannge is provided for \var{k}, the exact number of clusters will be determined within the range by Bayesian information criterion.  The Ckmeans.1d.dp algorithm groups 1-D data given by \code{x} into \code{k} cluster by dynamic programming (Wang and Song, 2011).  It guarantees the optimality of clustering -- the sum of \var{withinss} for each cluster is always the minimum.  In contrast, heuristic \var{k}-means algorithms may be inconsistent or non-optimal from run to run.  The run time of the algorithm is \var{O(max(k)n^2)}.
}

\value{
	An object of class "\code{Ckmeans.1d.dp}" which has a print method and is a list with components:
	\item{cluster}{A vector of cluster indices assigned to each element in \code{x}. Each cluster is indexed by an integer from 1 to k.}
	\item{centers}{A vector of cluster centres.}
	\item{withinss}{The within-cluster sum of squares for each cluster.}
	\item{size}{The number of points in each cluster.}
}

\references{
Wang, H. and Song, M. (2011) Ckmeans.1d.dp: optimal \var{k}-means clustering in one dimension by dynamic programming. \emph{The R Journal} \bold{3}(2), 29--33. Retrieved from \url{http://journal.r-project.org/archive/2011-2/RJournal_2011-2_Wang+Song.pdf}
}

\examples{
# Ex. 1 The number of clusters is provided. 
# Generate data from a Gaussian mixture model of two components
x <- c(rnorm(50, sd=0.3), rnorm(50, mean=1, sd=0.3))
# Divide x into 2 clusters
result <- Ckmeans.1d.dp(x, 2)
plot(x, col=result$cluster)
abline(h=result$centers, col=1:length(result$centers), pch=8, cex=2)

# Ex. 2 The number of clusters is determined by Bayesian information criterion
# Generate data from a Gaussian mixture model of two components
x <- c(rnorm(50, mean=-1, sd=0.3), rnorm(50, mean=1, sd=1))
# Divide x into k clusters, k automatically selected (default: 1~9)
result <- Ckmeans.1d.dp(x)
plot(x, col=result$cluster)
abline(h=result$centers, col=1:length(result$centers), pch=8, cex=2)
}