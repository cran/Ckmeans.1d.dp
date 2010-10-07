\name{Ckmeans.1d.dp}
\alias{Ckmeans.1d.dp}
\title{Optimal k-means clustering in one-dimension by dynamic programming}
\description{
Using dynamic programming to do one-dimensional clustering. Which will
guarantee optimality, i.e. the smallest "within-cluster sum of squares".
}
\usage{
Ckmeans.1d.dp(x, K)
}
\arguments{
  \item{x}{x is a one-dimensional input data, which will be clustered}
  \item{K}{K stands for the cluster level}
}
\author{
	Joe Song <joemsong@cs.nmsu.edu> 
	and 
	Haizhou Wang <hwang@cs.nmsu.edu>
}
\keyword{One-dimensional clustering}
\keyword{Dynamic programming}

\examples{

x <- rnorm(100,0,1)
Ckmeans.1d.dp(x,3)
}
