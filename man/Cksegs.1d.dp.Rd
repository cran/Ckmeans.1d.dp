\name{Univariate Segmentation}
\alias{Univariate Segmentation}
\alias{Cksegs.1d.dp}

\title{Optimal Univariate Segmentation}

\description{
Perform optimal univariate \eqn{k}{k}-segmentation.
}

\usage{
Cksegs.1d.dp(y, k=c(1,9), x=seq_along(y),
             method=c("quadratic", "linear", "loglinear"),
             estimate.k=c("BIC", "BIC 3.4.12"))
}

\arguments{

  \item{y}{a numeric vector of y values. Values can be negative.}

  \item{k}{either an exact integer number of clusters, or a vector of length two specifying the minimum and maximum numbers of clusters to be examined. The default is \code{c(1,9)}. When \code{k} is a range, the actual number of clusters is determined by Bayesian information criterion.}

  \item{x}{an optional numeric vector of data to be clustered. All \code{NA} elements must be removed from \code{x} before calling this function. The function will run faster on sorted \code{x} (in non-decreasing order) than an unsorted input.}

  \item{method}{a character string to specify the speedup method to the original cubic runtime dynamic programming. The default is \code{"quadratic"}, which generates optimal results. The other options do not guarantee optimal solution and differ in runtime or memory usage. See Details.}

  \item{estimate.k}{a character string to specify the method to estimate optimal \code{k}. The default is \code{"BIC"}. See Details.}

}

\author{
	Joe Song
}

\details{

\code{Cksegs.1d.dp} minimizes within-cluster sum of squared distance on \code{y}. It offers optimal piece-wise constant approximation of \code{y} within clusters of \code{x}. Only \code{method="quadratic"} guarantees optimality. The "linear" and "loglinear" options are faster but not always optimal and are provided for comparison purposes.

The Bayesian information criterion (BIC) method to select optimal \code{k} is updated to deal with duplicates in the data. Otherwise, the estimated k would be the same with previous versions. Set \code{estimate.k="BIC"} to use the latest method; use \code{estimate.k="BIC 3.4.12"} to use the BIC method in version 3.4.12 or earlier to estimated \code{k} from the given range. This option is effective only when a range for \code{k} is provided.

\code{method} specifies one of three options to speed up the original dynamic programming taking a runtime cubic in sample size \var{n}. The default \code{"quadratic"} option, with a runtime of \eqn{O(kn^2)}{O(kn^2)}, guarantees optimality. The next two options do not guarantee optimality. The \code{"linear"} option, giving a total runtime of \eqn{O(n \lg n + kn)}{O(n lg n + kn)} or \eqn{O(kn)}{O(kn)} (if \code{x} is already sorted in ascending order) is the fastest option but uses the most memory (still \eqn{O(kn)}{O(kn)}); the \code{"loglinear"} option, with a runtime of \eqn{O(kn \lg n)}{O(kn lg n)}, is slightly slower but uses the least memory.
}

\value{
	An object of class "\code{Cksegs.1d.dp}". It is a list containing the following components:

	\item{cluster}{a vector of clusters assigned to each element in \code{x}. Each cluster is indexed by an integer from 1 to \code{k}.}

	\item{centers}{a numeric vector of the (weighted) means for each cluster.}

	\item{withinss}{a numeric vector of the (weighted) within-cluster sum of squares for each cluster.}

	\item{size}{a vector of the (weighted) number of elements in each cluster.}

	\item{totss}{total sum of (weighted) squared distances between each element and the sample mean. This statistic is not dependent on the clustering result.}

	\item{tot.withinss}{total sum of (weighted) within-cluster squared distances between each element and its cluster mean. This statistic is minimized given the number of clusters.}

	\item{betweenss}{sum of (weighted) squared distances between each cluster mean and sample mean. This statistic is maximized given the number of clusters.}

	\item{xname}{a character string. The actual name of the \code{x} argument.}

	\item{yname}{a character string. The actual name of the \code{y} argument.}

The class has a print and a plot method: \code{\link{print.Cksegs.1d.dp}} and \code{\link{plot.Cksegs.1d.dp}}.
}

\seealso{
\code{\link{plot.Cksegs.1d.dp}} and \code{\link{print.Cksegs.1d.dp}}.
}

\examples{
# Ex 1. Segmenting by y

y <- c(1,1,1,2,2,2,4,4,4,4)

res <- Cksegs.1d.dp(y, k=c(1:10))

main <- "k-segs giving 3 clusters\nsucceeded in finding segments"

opar <- par(mfrow=c(1,2))

plot(res, main=main, xlab="x")

res <- Ckmeans.1d.dp(x=seq_along(y), k=c(1:10), y)
main <- "Weighted k-means giving 1 cluster\nfailed to find segments"

plot(res, main=main, xlab="x")

par(opar)

# Ex 2. Segmenting by y

y <- c(1,1,1.1,1, 2,2.5,2, 4,5,4,4)
res <- Cksegs.1d.dp(y, k=c(1:10))
plot(res, xlab="x")

# Ex 3. Segmenting a sinusoidal curve by y
x <- 1:125
y <- sin(x * .2)
res.q <- Cksegs.1d.dp(y, k=8, x=x)
plot(res.q, lwd=3, xlab="x")

# Ex 4. Segmenting by y

y <- rep(c(1,-3,4,-2), each=20)
y <- y + 0.5*rnorm(length(y))
k <- 1:10
res.q <- Cksegs.1d.dp(y, k=k, method="quadratic")
main <- paste("Cksegs (method=\"quadratic\"):\ntot.withinss =",
              format(res.q$tot.withinss, digits=4), "BIC =",
              format(res.q$BIC[length(res.q$size)], digits=4),
              "\nGUARANTEE TO BE OPTIMAL")
plot(res.q, main=main, xlab="x")
res.l <- Cksegs.1d.dp(y, k=k, method="linear")
main <- paste("Cksegs (method=\"linear\"):\ntot.withinss =",
               format(res.l$tot.withinss, digits=4), "BIC =",
              format(res.l$BIC[length(res.l$size)], digits=4),
               "\nFAST BUT MAY NOT BE OPTIMAL")
plot(res.l, main=main, xlab="x")
res.g <- Cksegs.1d.dp(y, k=k, method="loglinear")
main <- paste("Cksegs (method=\"loglinear\"):\ntot.withinss =",
              format(res.g$tot.withinss, digits=4), "BIC =",
              format(res.g$BIC[length(res.g$size)], digits=4),
              "\nFAST BUT MAY NOT BE OPTIMAL")
plot(res.g, main=main, xlab="x")
}


\keyword{ univar }
\keyword{ cluster }

