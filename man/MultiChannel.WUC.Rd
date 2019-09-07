\name{MultiChannel.WUC}
\alias{MultiChannel.WUC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Optimal Multi-channel Weighted Univariate Clustering
}
\description{
  Perform optimal multi-channel weighted univariate \eqn{k}{k}-means clustering in linear time.
}
\usage{
  MultiChannel.WUC(x, y,  k=c(1,9))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a numeric vector of data to be clustered. All \code{NA} elements must be removed from \code{x} before calling this function. The function will run faster on sorted \code{x} (in non-decreasing order) than an unsorted input.}

  \item{y}{a numeric matrix of non-negative weights for each element in \code{x}. Columns of the matrix are channels. It is highly recommended to use positive (instead of zero) weights to account for the influence of every element. Weights strongly influence clustering results. When the number of clusters \code{k} is given as a range, the weights should be linearly scaled to sum up to the observed sample size.}

  \item{k}{either an exact integer number of clusters, or a vector of length two specifying the minimum and maximum numbers of clusters to be examined. The default is \code{c(1,9)}. When \code{k} is a range, the actual number of clusters is determined by Bayesian information criterion (BIC).}
}
\details{
  \code{MultiChannel.WUC} minimizes the total weighted within-cluster sum of squared distance. It uses the SMAWK algorithm \insertCite{aggarwal1987geometric}{Ckmeans.1d.dp} with modified data structure to speed up the dynamic programming to linear runtime. The method selects an optimal \code{k} based on an approximate Gaussian mixture model using the BIC.
}
\value{
  A list object containing the following components:

	\item{cluster}{a vector of clusters assigned to each element in \code{x}. Each cluster is indexed by an integer from 1 to \code{k}.}

	\item{centers}{a numeric vector of the (weighted) means for each cluster.}

	\item{withinss}{a numeric vector of the (weighted) within-cluster sum of squares for each cluster.}

	\item{size}{a vector of the (weighted) number of elements in each cluster.}

	\item{totss}{total sum of (weighted) squared distances between each element and the sample mean. This statistic is not dependent on the clustering result.}

	\item{tot.withinss}{total sum of (weighted) within-cluster squared distances between each element and its cluster mean. This statistic is minimized given the number of clusters.}

	\item{betweenss}{sum of (weighted) squared distances between each cluster mean and sample mean. This statistic is maximized given the number of clusters.}

	\item{xname}{a character string. The actual name of the \code{x} argument.}

	\item{yname}{a character string. The actual name of the \code{y} argument.}
}
\references{
  \insertAllCited{}
}
\author{
  Hua Zhong and Mingzhou Song
}
\examples{
  X <- sort(sample(x = c(1:100), size = 20, replace = TRUE))
  Y <- matrix(sample(x = c(1:100), size = 40, replace = TRUE), ncol=2, nrow=length(X))

  res <- MultiChannel.WUC(x = X, y = Y, k = c(5:10))
}
\keyword{ univar }
\keyword{ cluster }
