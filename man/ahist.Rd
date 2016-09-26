\name{ahist}
\alias{ahist}
\title{Adaptive Histograms by Optimal Univariate k-Means Clustering}
\description{
Generate histograms that are adaptive to patterns in data. The number and widths of histogram bins are automatically adjusted to the data and thus the bins are unlikely of equal width. The breaks are computed from optimal univariate k-means clusters of the input data.
}
\usage{
ahist(x, k=c(1,9), plot = TRUE, xlab = deparse(substitute(x)),
  main = paste("Adaptive histogram of", deparse(substitute(x))),
  \dots)
}
\arguments{
  \item{x}{a numeric vector of data. All \code{NA} elements must be removed from \code{x} before calling this function. The function will run faster on sorted \code{x} (in non-decreasing order) than an unsorted input.}

  \item{k}{either an exact integer number of clusters, or a vector of length two specifying the minimum and maximum numbers of clusters to be examined. The default is \code{c(1,9)}. When \code{k} is a range, the actual number of clusters is determined by Bayesian information criterion.}

 \item{plot}{a logical. If \code{TRUE}, the histogram will be plotted.}

 \item{xlab}{a character string of the x-axis lab for the plot}

 \item{main}{a character string of the title for the plot}
 \item{...}{additional arguments passed to \code{\link{hist}} and \code{\link{plot.histogram}}.}
}

\author{
	Joe Song
}

\details{
The number of histogram bins is the optimal number of univariate k-means clusters estimated using Bayesian information criterion evaluated on Gaussian mixture models fitted to the input data in \code{x}. Breaks in the histogram are the midpoints between clusters identified by optimal univariate k-means (\code{\link{Ckmeans.1d.dp}}). The histogram is by default plotted using the \code{\link{plot.histogram}} method. The plot can be optionally disabled with the \code{plot=FALSE} argument.
}

\value{
An object of class \code{histogram} defined in \code{\link{hist}}. It has a S3 \code{plot} method \code{\link{plot.histogram}}.
}

\references{
Wang, H. and Song, M. (2011) Ckmeans.1d.dp: optimal \var{k}-means clustering in one dimension by dynamic programming. \emph{The R Journal} \bold{3}(2), 29--33. Retrieved from \url{http://journal.r-project.org/archive/2011-2/RJournal_2011-2_Wang+Song.pdf}
}

\seealso{
  \code{\link{hist}} in package \pkg{\link{graphics}}.
}

\examples{
# Example: plot an optimal histogram from a Gaussian mixture
#          model with three components
x <- c(rnorm(40, mean=-1, sd=0.3),
       rnorm(25, mean=1, sd=0.3),
       rnorm(60, mean=3, sd=0.2))
ahist(x, col="seagreen")
}
