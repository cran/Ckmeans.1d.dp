\name{plot.Ckmeans.1d.dp}
\alias{plot.Ckmeans.1d.dp}
\alias{plot.Ckmedian.1d.dp}
\alias{plot.Cksegs.1d.dp}
\title{Plot Optimal Univariate Clustering Results}
\description{
Plot optimal univariate clustering results returned from \code{Ckmeans.1d.dp}.
}
\usage{
\method{plot}{Ckmeans.1d.dp}(x, xlab=NULL, ylab=NULL, main=NULL,
      sub=NULL, col.clusters=NULL, \dots)
\method{plot}{Ckmedian.1d.dp}(x, xlab=NULL, ylab=NULL, main=NULL,
      sub=NULL, col.clusters=NULL, \dots)
\method{plot}{Cksegs.1d.dp}(x, xlab=NULL, ylab=NULL, main=NULL,
      sub=NULL, col.clusters=NULL, \dots)
}
%(x, xlab=ck$xname,
%     ylab=ifelse(x$yname=="1", "Weight", x$yname),
%     main="Optimal univariate clustering of", x$xname),
%     sub=paste("n =", length(x$cluster)),
%     col.clusters=seq_along(x$size),
%     \dots)

\arguments{
  \item{x}{an object of class as returned by \code{\link{Ckmeans.1d.dp}}, \code{\link{Ckmedian.1d.dp}}, or \code{\link{Cksegs.1d.dp}}.}

 \item{xlab}{a character string. The x-axis label for the plot.}

  \item{ylab}{a character string. The x-axis label for the plot.}

 \item{main}{a character string. The title for the plot.}

  \item{sub}{a character string. The subtitle for the plot.}

  \item{col.clusters}{a vector of colors, defined either by integers or by color names. If the length is shorter than the number of clusters, the colors will be reused.}

  \item{...}{arguments passed to \code{\link{plot}} function in package \pkg{graphics}.}
}

\author{
	Joe Song
}

\details{
The functions \code{plot.Ckmeans.1d.dp} and \code{plot.Ckmedian.1d.dp} visualize the input data as sticks whose heights are the weights. They use different colors to indicate clusters.

The function \code{plot.Cksegs.1d.dp} shows segments as horizontal lines from the univariate
segmentation results obtained from function \code{Cksegs.1d.dp}. It uses different colors to indicate segments.
}

\value{
	An object of class "\code{Ckmeans.1d.dp}", "\code{Ckmedian.1d.dp}", or "\code{Cksegs.1d.dp}" defined in \code{\link{Ckmeans.1d.dp}}.
}

\references{
  Wang, H. and Song, M. (2011) Ckmeans.1d.dp: optimal \var{k}-means clustering in one dimension by dynamic programming. \emph{The R Journal} \bold{3}(2), 29--33. Retrieved from \url{https://journal.r-project.org/archive/2011-2/RJournal_2011-2_Wang+Song.pdf}
}

\examples{
# Example: clustering data generated from a Gaussian
#          mixture model of three components
x <- c(rnorm(50, mean=-1, sd=0.3),
       rnorm(50, mean=1, sd=0.3),
       rnorm(50, mean=3, sd=0.3))

res <- Ckmeans.1d.dp(x)
plot(res)

y <- (rnorm(length(x)))^2
res <- Ckmeans.1d.dp(x, y=y)
plot(res)

res <- Ckmedian.1d.dp(x)
plot(res)

y <- x^3
res <- Cksegs.1d.dp(y, x=x)
plot(res, lwd=2)

}

\keyword{ univar }
\keyword{ cluster }
\keyword{ hplot }
\keyword{ distribution }

