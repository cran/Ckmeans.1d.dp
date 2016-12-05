\name{Ckmeans.1d.dp-package}
\alias{Ckmeans.1d.dp-package}
\docType{package}
\title{
Optimal and Fast Univariate \var{k}-Means Clustering
}
\description{
The Ckmeans.1d.dp algorithm clusters univariate data given by a numeric vector \eqn{x}{x} into \eqn{k}{k} groups by dynamic programming (Wang and Song, 2011). It guarantees the optimality of clustering---the total of within-cluster sums of squares is always the minimum given the number of clusters \eqn{k}{k}. In contrast, heuristic \eqn{k}{k}-means algorithms may be non-optimal or inconsistent from run to run. Apart from the time for sorting \eqn{x}{x}, the algorithm takes a runtime linear in number of clusters and log-linear in sample size, comparable to heuristic \eqn{k}{k}-means. The space complexity is \eqn{O(kn)}{O(kn)}. It is practical for Ckmeans.1d.dp to cluster millions of sample points within seconds using a single processor on a not very recent desktop computer.

A major speed-up has been introduced since version 3.4.6 using a divide-and-conquer strategy in dynamic programming. The runtime is reduced from \eqn{O(kn^2)}{O(kn^2)} down to \eqn{O(kn\lg n)}{O(kn lg n)}, log-linear in sample size \eqn{n}{n}. This method is numerically stable and observed runtime is often close to linear in practice. These improvements, not discussed in (Wang and Song, 2011), will be described in detail in a future publication. The weighted optimal clustering option is still quadratic in \eqn{n}{n} and will be updated to loglinear in a future release.

Richard Bellman (1973) first described a general dynamic programming strategy for solving univariate clustering problems with additive optimality measures. The strategy, however, did not address any specific characteristics of the \eqn{k}{k}-means problem and its implied general algorithm will have a time complexity of \eqn{O(kn^3)}{O(k n^3)} and space complexity of \eqn{O(kn)}{O(kn)} on an input of \eqn{n}{n} points.

Utilizing the optimal clusters, a function can generate histograms adaptive to patterns in data.

As unequal non-negative weights are supported for each point, the algorithm can also segment a time course using the time points as input and the values at each time point as weight. For univariate data, the package provides a powerful alternative to heuristic clustering algorithms and also new functionality for weighted clustering / segmentation.
}

\details{
\tabular{ll}{
Package: \tab Ckmeans.1d.dp\cr
Type: \tab Package\cr
Version: \tab 3.4.6-5\cr
Initial version: \tab 1.0\cr
Initial date: \tab 2010-10-26\cr
License: \tab LGPL (>= 3) \cr
}

}

\seealso{
  The \code{\link{kmeans}} function in package \pkg{\link{stats}} that implements several heuristic \eqn{k}{k}-means algorithms.
}

\author{
Joe Song and Haizhou Wang
}

\references{
Wang, H. and Song, M. (2011) Ckmeans.1d.dp: optimal \var{k}-means clustering in one dimension by dynamic programming. \emph{The R Journal} \bold{3}(2), 29--33. Retrieved from \url{https://journal.r-project.org/archive/2011-2/RJournal_2011-2_Wang+Song.pdf}

Bellman, R. (1973) A note on cluster analysis and dynamic programming. \emph{Mathematical Biosciences} \bold{18}(3), 311--312.
}

\keyword{ package }
