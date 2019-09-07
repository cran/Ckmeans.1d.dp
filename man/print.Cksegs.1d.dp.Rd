\name{print.Cksegs.1d.dp}
\alias{print.Cksegs.1d.dp}

\title{Print Optimal Univariate Segmentation Results}
\description{
Print optimal univariate segmentation results obtained from \code{Cksegs.1d.dp}.
}
\usage{
\method{print}{Cksegs.1d.dp}(x, \dots)
}
\arguments{
  \item{x}{object returned by calling \code{Cksegs.1d.dp}.}
  \item{...}{arguments passed to \code{\link{print}} function.}
}

\author{
	Joe Song
}

\details{
Function \code{print.Cksegs.1d.dp} prints the maximum ratio of between-cluster sum of squares to total sum of squares unless all input elements are zero. The ratio is an indicator of maximum achievable clustering quality given the number of clusters: 100\% for a perfect clustering and 0\% for no cluster patterns.
}

\value{
	An object of class "\code{Cksegs.1d.dp}" as defined in \code{\link{Cksegs.1d.dp}}.
}

\examples{
# Example: clustering data generated from a Gaussian
#          mixture model of two components
x <- c(rnorm(15, mean=-1, sd=0.3),
       rnorm(15, mean=1, sd=0.3))

y <- x^3
res <- Cksegs.1d.dp(y, x=x)
print(res, lwd=2)
}

\keyword{ univar }
\keyword{ cluster }
\keyword{ print }
