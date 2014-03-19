\name{print.Ckmeans.1d.dp}
\alias{print.Ckmeans.1d.dp}
\title{Print Results from Ckmeans.1d.dp}
\description{
Print the result returned by calling Ckmeans.1d.dp
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
x <- rnorm(50, mean=-1, sd=0.3)
x <- append(x, rnorm(50, mean=1, sd=0.3) )
res <- Ckmeans.1d.dp(x)
print(res)

}