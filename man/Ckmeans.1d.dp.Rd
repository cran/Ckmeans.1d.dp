\name{Ckmeans.1d.dp}
\alias{Ckmeans.1d.dp}
\title{Optimal and Fast Univariate Clustering}

\description{
Perform optimal univariate clustering in linear (fastest), loglinear, or quadratic (slowest) time.
}

\usage{
Ckmeans.1d.dp(x, k=c(1,9), y=1,
              method=c("linear", "loglinear", "quadratic"),
              estimate.k=c("BIC", "BIC 3.4.12"))
}

\arguments{
  \item{x}{a numeric vector of data to be clustered. All \code{NA} elements must be removed from \code{x} before calling this function. The function will run faster on sorted \code{x} (in non-decreasing order) than an unsorted input.}

  \item{k}{either an exact integer number of clusters, or a vector of length two specifying the minimum and maximum numbers of clusters to be examined. The default is \code{c(1,9)}. When \code{k} is a range, the actual number of clusters is determined by Bayesian information criterion.}

  \item{y}{a value of 1 to specify equal weights or a numeric vector of unequal weights for each element. The default weight is one. It is highly recommended to use positive (instead of zero) weights to account for the influence of every element. The weights have a strong impact on the clustering result.}

  \item{method}{a character string to specify the speedup method to the original cubic runtime dynamic programming. The default is \code{"linear"}. All methods generate the same optimal results but differ in runtime or memory usage. See Details.}

  \item{estimate.k}{a character string to specify the method to estimate optimal \code{k}. The default is \code{"BIC"}. See Details.}

}

\author{
	Joe Song and Haizhou Wang
}

\details{
In contrast to the heuristic \var{k}-means algorithms implemented in function \code{\link{kmeans}}, this function optimally assigns elements in numeric vector \code{x} into \code{k} clusters by dynamic programming (Wang and Song, 2011). It minimizes the total of within-cluster sums of squared distances (\var{withinss}) between each element and its corresponding cluster mean. When a range is provided for \code{k}, the exact number of clusters is determined by Bayesian information criterion. Different from the heuristic \var{k}-means algorithms whose results may be non-optimal or change from run to run, the result of Ckmeans.1d.dp is guaranteed to be optimal and reproducible, and its advantage in efficiency and accuracy over heuristic \var{k}-means methods is most pronounced at large \var{k}.

The Bayesian information criterion (BIC) method to select optimal \code{k} is updated to deal with duplicates in the data. Otherwise, the estimated k would be the same with previous versions. Set \code{estimate.k="BIC"} to use the lastest method; use \code{estimate.k="BIC 3.4.12"} to use the BIC method in version 3.4.12 or earlier to estimated \code{k} from the given range. This option is effective only when a range for \code{k} is provided.

\code{method} specifies one of three options to speed up the original dynamic programming taking a runtime cubic in sample size \var{n}. The default \code{"linear"} option, giving a total runtime of \eqn{O(n \lg n + kn)}{O(n lg n + kn)} or \eqn{O(kn)}{O(kn)} (if \code{x} is already sorted in ascending order) is the fastest option but uses the most memory (still \eqn{O(kn)}{O(kn)}); the \code{"loglinear"} option, with a runtime of \eqn{O(kn \lg n)}{O(kn lg n)}, is slightly slower but uses the least memory; the slowest \code{"quadratic"} option, with a runtime of \eqn{O(kn^2)}{O(kn^2)}, is provided for the purpose of testing on small data sets.

When the sample size \var{n} is too large to create two \eqn{k \times n}{k x n} dynamic programming matrices in memory, we recommend the heuristic solutions impelmented in the \code{\link{kmeans}} function in package \pkg{\link{stats}}.
}

\value{
	An object of class \code{"Ckmeans.1d.dp"} with a \code{\link{print.Ckmeans.1d.dp}} and a \code{\link{plot.Ckmeans.1d.dp}} method. It is a list containing the following components:

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
Wang, H. and Song, M. (2011) Ckmeans.1d.dp: optimal \var{k}-means clustering in one dimension by dynamic programming. \emph{The R Journal} \bold{3}(2), 29--33. Retrieved from \url{https://journal.r-project.org/archive/2011-2/RJournal_2011-2_Wang+Song.pdf}
}

\seealso{
\code{\link{ahist}}, \code{\link{plot.Ckmeans.1d.dp}}, \code{\link{print.Ckmeans.1d.dp}} in this package.

\code{\link[stats]{kmeans}} in package \pkg{\link{stats}} that implements several heuristic \eqn{k}{k}-means algorithms.
}

\examples{
# Ex. 1 The number of clusters is provided.
# Generate data from a Gaussian mixture model of three components
x <- c(rnorm(50, sd=0.2), rnorm(50, mean=1, sd=0.3), rnorm(100, mean=-1, sd=0.25))
# Divide x into 3 clusters
k <- 3
result <- Ckmeans.1d.dp(x, k)

plot(result)

plot(x, col=result$cluster, pch=result$cluster, cex=1.5,
     main="Optimal univariate clustering given k",
     sub=paste("Number of clusters given:", k))
abline(h=result$centers, col=1:k, lty="dashed", lwd=2)
legend("bottomleft", paste("Cluster", 1:k), col=1:k, pch=1:k, cex=1.5, bty="n")

# Ex. 2 The number of clusters is determined by Bayesian information criterion
# Generate data from a Gaussian mixture model of three components
x <- c(rnorm(50, mean=-3, sd=1), rnorm(50, mean=0, sd=.5), rnorm(50, mean=3, sd=1))
# Divide x into k clusters, k automatically selected (default: 1~9)
result <- Ckmeans.1d.dp(x)

plot(result)

k <- max(result$cluster)
plot(x, col=result$cluster, pch=result$cluster, cex=1.5,
     main="Optimal univariate clustering with k estimated",
     sub=paste("Number of clusters is estimated to be", k))
abline(h=result$centers, col=1:k, lty="dashed", lwd=2)
legend("topleft", paste("Cluster", 1:k), col=1:k, pch=1:k, cex=1.5, bty="n")

# Ex. 3 Segmenting a time course using optimal weighted univariate clustering
n <- 160
t <- seq(0, 2*pi*2, length=n)
n1 <- 1:(n/2)
n2 <- (max(n1)+1):n
y1 <- abs(sin(1.5*t[n1]) + 0.1*rnorm(length(n1)))
y2 <- abs(sin(0.5*t[n2]) + 0.1*rnorm(length(n2)))
y <- c(y1, y2)

w <- y^8 # stress the peaks
res <- Ckmeans.1d.dp(t, k=c(1:10), w)
plot(res)
plot(t, w, main = "Time course segmentation",
     col=res$cluster, pch=res$cluster,
     xlab="Time t", ylab="Transformed intensity w",
     type="h")
abline(v=res$centers, col="chocolate", lty="dashed")
text(res$centers, max(w) * .95, cex=0.5, font=2,
     paste(round(res$size / sum(res$size) * 100), "/ 100"))
}
