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
  \item{x}{object returned by calling Ckmeans.1d.dp.}
  \item{...}{ignored arguments.}
}

\author{
	Joe Song and Haizhou Wang
}

\details{
The maximum ratio of between-cluster sum of squares to total sum of squares is printed unless all input elements are zero. The ratio is an indicator of maximum achievable clustering quality given the number of clusters: 100\% for a perfect clustering and 0\% for no cluster patterns.
}

\value{
	An object of class "\code{Ckmeans.1d.dp}" which has a \code{print} method and is a list with components:

	\item{cluster}{a vector of clusters assigned to each element in \code{x}. Each cluster is indexed by an integer from 1 to \code{k}.}

	\item{centers}{a numeric vector of the (weighted) means for each cluster.}

	\item{withinss}{a numeric vector of the (weighted) within-cluster sum of squares for each cluster.}

	\item{size}{a vector of the (weighted) number of elements in each cluster.}

	\item{totss}{total sum of squared distances between each element and the sample mean. This statistic is not dependent on the clustering result.}

	\item{tot.withinss}{total sum of within-cluster squared distances between each element and its cluster mean. This statistic is minimized given the number of clusters.}

	\item{betweenss}{sum of squared distances between each cluster mean and sample mean weighed by cluster size. This statistic is maximized given the number of clusters.}
}

\references{
  Wang, H. and Song, M. (2011) Ckmeans.1d.dp: optimal \var{k}-means clustering in one dimension by dynamic programming. \emph{The R Journal} \bold{3}(2), 29--33. Retrieved from \url{http://journal.r-project.org/archive/2011-2/RJournal_2011-2_Wang+Song.pdf}
}

\examples{
# Example: clustering data generated from a Gaussian
#          mixture model of two components
x <- c(rnorm(15, mean=-1, sd=0.3),
       rnorm(15, mean=1, sd=0.3))
res <- Ckmeans.1d.dp(x)
print(res)
}
