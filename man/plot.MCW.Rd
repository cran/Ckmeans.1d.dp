\name{plot.MultiChannelClusters}
\alias{plot.MultiChannelClusters}
\title{Plot Multi-Channel Clustering Results}
\description{
Plot multi-channel clustering results returned from \code{MultiChannel.WUC}.
}
\usage{
\method{plot}{MultiChannelClusters}(x, xlab=NULL, ylab=NULL, main=NULL,
      sub=NULL, col.clusters=NULL, \dots)
}

\arguments{
  \item{x}{an object of class as returned by \code{\link{MultiChannel.WUC}}.}

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
The function visualizes the input data as sticks whose heights are the weights. Colors indicate clusters. See \code{\link{MultiChannel.WUC}} for examples.
}

\value{
	An object of class "\code{MultiChannelClusters}" defined in \code{\link{MultiChannel.WUC}}.
}

\keyword{ cluster }
\keyword{ hplot }
\keyword{ distribution }
