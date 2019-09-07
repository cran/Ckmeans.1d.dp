\name{print.Ckmeans.1d.dp}
\alias{print.Ckmeans.1d.dp}
\alias{print.Ckmedian.1d.dp}

\title{Print Optimal Univariate Clustering Results}
\description{
Print optimal univariate clustering results obtained from \code{Ckmeans.1d.dp} or \code{Ckmedian.1d.dp}.
}

\usage{
\method{print}{Ckmeans.1d.dp}(x, \dots)
\method{print}{Ckmedian.1d.dp}(x, \dots)
}

\arguments{
  \item{x}{object returned by calling \code{Ckmeans.1d.dp} or \code{Cksegs.1d.dp}.}
  \item{...}{arguments passed to \code{\link{print}} function.}
}

\author{
	Joe Song and Haizhou Wang
}

\details{
Function \code{print.Ckmeans.1d.dp} and \code{print.Ckmedian.1d.dp} prints the maximum ratio of between-cluster sum of squares to total sum of squares unless all input elements are zero. The ratio is an indicator of maximum achievable clustering quality given the number of clusters: 100\% for a perfect clustering and 0\% for no cluster patterns.
}

\value{
	An object of class "\code{Ckmeans.1d.dp}" or "\code{Ckmedian.1d.dp}" as defined in \code{\link{Ckmeans.1d.dp}}.
}

\examples{
# Example: clustering data generated from a Gaussian
#          mixture model of two components
x <- c(rnorm(15, mean=-1, sd=0.3),
       rnorm(15, mean=1, sd=0.3))
res <- Ckmeans.1d.dp(x)
print(res)

res <- Ckmedian.1d.dp(x)
print(res)

y <- (rnorm(length(x)))^2
res <- Ckmeans.1d.dp(x, y=y)
print(res)

res <- Ckmedian.1d.dp(x)
print(res)
}

\keyword{ univar }
\keyword{ cluster }
\keyword{ print }
