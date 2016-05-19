\name{print.Ckmeans.1d.dp}
\alias{print.Ckmeans.1d.dp}
\title{Print Results from Ckmeans.1d.dp}
\description{
Print the result returned from Ckmeans.1d.dp
}
\usage{
\method{print}{Ckmeans.1d.dp}(x, \dots)
}
\arguments{
  \item{x}{ object returned by calling Ckmeans.1d.dp}
  \item{...}{ Ignored arguments }
}

\author{
	Joe Song and Haizhou Wang
}

\details{
  The maximum ratio of between-cluster sum of squares to total sum of squares is printed unless all input numbers are zero. The ratio is an indicator of clustering quality: 100\% for perfect clustering and 0\% for worst clustering.
}

\value{
	An object of class "\code{Ckmeans.1d.dp}" which has a print method and is a list with components:
	\item{cluster}{a vector of integers (\code{1:k}) indicating the cluster to which each point is allocated.}
	\item{centers}{a numeric vector of cluster means.}
	\item{withinss}{a numeric vector of within-cluster sums of squares for each cluster.}
	\item{size}{a vector of the number of points in each cluster.}
	\item{totss}{total sum of squares of the input numbers.}
	\item{tot.withinss}{total sum of within-cluster distance squares.}
	\item{betweenss}{between-cluster sum of squares, equal to the sum of squared cluster means weighed by cluster size.}
}

\references{
  Wang, H. and Song, M. (2011) Ckmeans.1d.dp: optimal \var{k}-means clustering in one dimension by dynamic programming. \emph{The R Journal} \bold{3}(2), 29--33. Retrieved from \url{http://journal.r-project.org/archive/2011-2/RJournal_2011-2_Wang+Song.pdf}
}

\examples{
# Example: clustering data generated from a Gaussian mixture model of two components
x <- rnorm(50, mean=-1, sd=0.3)
x <- append(x, rnorm(50, mean=1, sd=0.3) )
res <- Ckmeans.1d.dp(x)
print(res)
}
