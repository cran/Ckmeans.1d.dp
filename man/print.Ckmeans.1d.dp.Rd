\name{print.Ckmeans.1d.dp}
\alias{print.Ckmeans.1d.dp}
\alias{print.Ckmedian.1d.dp}
\alias{print.Cksegs.1d.dp}

\title{Print Optimal Univariate Clustering Results}
\description{
Print optimal univariate clustering results obtained from \code{Ckmeans.1d.dp}.
}
\usage{
\method{print}{Ckmeans.1d.dp}(x, \dots)
\method{print}{Ckmedian.1d.dp}(x, \dots)
\method{print}{Cksegs.1d.dp}(x, \dots)
}
\arguments{
  \item{x}{object returned by calling Ckmeans.1d.dp.}
  \item{...}{arguments passed to \code{\link{print}} function.}
}

\author{
	Joe Song and Haizhou Wang
}

\details{
Function \code{print.Ckmeans.1d.dp} prints the maximum ratio of between-cluster sum of squares to total sum of squares unless all input elements are zero. The ratio is an indicator of maximum achievable clustering quality given the number of clusters: 100\% for a perfect clustering and 0\% for no cluster patterns.
}

\value{
	An object of class "\code{Ckmeans.1d.dp}", "\code{Ckmedian.1d.dp}", or "\code{Cksegs.1d.dp}" as defined in \code{\link{Ckmeans.1d.dp}}.
}

\references{
  Wang, H. and Song, M. (2011) Ckmeans.1d.dp: optimal \var{k}-means clustering in one dimension by dynamic programming. \emph{The R Journal} \bold{3}(2), 29--33. Retrieved from \url{https://journal.r-project.org/archive/2011-2/RJournal_2011-2_Wang+Song.pdf}
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

y <- x^3
res <- Cksegs.1d.dp(y, x=x)
print(res, lwd=2)
}

\keyword{ univar }
\keyword{ cluster }
\keyword{ print }
