\name{plot.Ckmeans.1d.dp}
\alias{plot.Ckmeans.1d.dp}
\title{Plot Optimal Univariate \var{k}-Means Clustering Results}
\description{
Plot optimal univariate clustering results returned from \code{Ckmeans.1d.dp}.
}
\usage{
\method{plot}{Ckmeans.1d.dp}(x, xlab=NULL, ylab=NULL, main=NULL,
      sub=NULL, col.clusters=NULL, \dots)
%(x, xlab=ck$xname,
%     ylab=ifelse(x$yname=="1", "Weight", x$yname),
%     main="Optimal k-means clustering of", x$xname),
%     sub=paste("n =", length(x$cluster)),
%     col.clusters=seq_along(x$size),
%     \dots)
}
\arguments{
  \item{x}{an object of class \code{Ckmeans.1d.dp} returned by \code{\link{Ckmeans.1d.dp}}.}

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
The function visualizes the input data as sticks whose heights are the weights. It uses different colors to indicate optimal \var{k}-means clusters.
}

\value{
	An object of class "\code{Ckmeans.1d.dp}" defined in \code{\link{Ckmeans.1d.dp}}.
}

\references{
  Wang, H. and Song, M. (2011) Ckmeans.1d.dp: optimal \var{k}-means clustering in one dimension by dynamic programming. \emph{The R Journal} \bold{3}(2), 29--33. Retrieved from \url{http://journal.r-project.org/archive/2011-2/RJournal_2011-2_Wang+Song.pdf}
}

\examples{
# Example: clustering data generated from a Gaussian mixture model of two components
x <- rnorm(50, mean=-1, sd=0.3)
x <- append(x, rnorm(50, mean=1, sd=0.3) )
res <- Ckmeans.1d.dp(x)
plot(res)

y <- (rnorm(length(x)))^2
res <- Ckmeans.1d.dp(x, y=y)
plot(res)
}
