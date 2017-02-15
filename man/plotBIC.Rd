\name{plotBIC}
\alias{plotBIC}
\title{Plot Bayesian Information Criterion as a Function of Number of Clusters}
\description{
Plot Bayesian information criterion (BIC) as a function of the number of clusters obtained from optimal univariate clustering results returned from \code{Ckmeans.1d.dp}. The BIC normalized by sample size (BIC/n) is shown.
}
\usage{
plotBIC(ck, xlab="Number of clusters k",
        ylab = "BIC/n", type="b",
        sub=paste("n =", length(ck$cluster)),
        main="Bayesian information criterion\n(normalized by sample size)",
        ...)
}

\arguments{
  \item{ck}{an object of class \code{Ckmeans.1d.dp} returned by \code{\link{Ckmeans.1d.dp}}.}

 \item{xlab}{a character string. The x-axis label for the plot.}

 \item{ylab}{a character string. The x-axis label for the plot.}

 \item{type}{the type of plot to be drawn. See \code{\link[graphics]{plot}}.}

 \item{main}{a character string. The title for the plot.}

 \item{sub}{a character string. The subtitle for the plot.}

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
  Wang, H. and Song, M. (2011) Ckmeans.1d.dp: optimal \var{k}-means clustering in one dimension by dynamic programming. \emph{The R Journal} \bold{3}(2), 29--33. Retrieved from \url{https://journal.r-project.org/archive/2011-2/RJournal_2011-2_Wang+Song.pdf}
}

\examples{
# Example: clustering data generated from a Gaussian mixture model of two components
x <- rnorm(50, mean=-1, sd=0.3)
x <- append(x, rnorm(50, mean=1, sd=0.3) )
res <- Ckmeans.1d.dp(x)
plotBIC(res)

y <- (rnorm(length(x)))^2
res <- Ckmeans.1d.dp(x, y=y)
plotBIC(res)
}
