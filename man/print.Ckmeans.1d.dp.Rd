\name{print.Ckmeans.1d.dp}
\alias{print.Ckmeans.1d.dp}
\title{print.Ckmeans.1d.dp}
\description{
Print the result returned by calling Ckmeans.1d.dp
}
\usage{
print.Ckmeans.1d.dp(x, ...)
}
\arguments{
  \item{x}{ object returned by calling Ckmeans.1d.dp}
  \item{...}{ Ignored arguments }
}
\author{
	Joe Song <joemsong@cs.nmsu.edu> 
	and 
	Haizhou Wang <hwang@cs.nmsu.edu>
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
Ckmeans.1d.dp(x, 2)
}