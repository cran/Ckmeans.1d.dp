\name{Ckmeans.1d.dp-package}
\alias{Ckmeans.1d.dp-package}
\docType{package}
\title{
Optimal and Fast Univariate \var{k}-Means Clustering
}
\description{
The Ckmeans.1d.dp algorithm clusters 1-D data given by a numeric vector \eqn{x}{x} into \eqn{k}{k} groups by dynamic programming (Wang and Song, 2011). It guarantees the optimality of clustering -- the total of within-cluster sum of squares is always the minimum given the number of clusters. In contrast, heuristic \eqn{k}{k}-means algorithms may be non-optimal or inconsistent from run to run. Apart from sorting in log-linear time, clustering contributes only linear in number of clusters and log-linear in sample size, comparable to heuristic \eqn{k}{k}-means when no-restart is allowed. It is practical for Ckmeans.1d.dp to cluster millions of sample points within seconds using a single processor on a recent desktop computer.

A major speed-up has been introduced since version 3.4.6 using a divide-and-conquer strategy in dynamic programming. The runtime is reduced from \eqn{O(kn^2)}{O(kn^2)} down to \eqn{O(kn\lg n)}{O(kn lg n)}, log-linear in sample size \eqn{n}{n}. This method is numerically stable and observed runtime is often close to linear in practice. These improvements, not discussed in (Wang and Song, 2011), will be described in detail in a future publication.

Richard Bellman (1973) first described a general dynamic programming strategy for solving univariate clustering problems with additive optimality measures. The strategy, however, did not address any specific characteristics of the \eqn{k}{k}-means problem and its implied general algorithm will have a time complexity of \eqn{O(kn^3)}{O(k n^3)} and space complexity of \eqn{O(kn)}{O(kn)} on an input of \eqn{n}{n} points.

%Since version 3.4.0, the runtime of this dynamic programming algorithm has been greatly reduced after two changes to trim down the search space. First, an upper bound for the sum of within cluster square distances is checked to reduce the dynamic programming search space. Second, the unnecessary calculation of \eqn{n-1}{n-1} elements in the dynamic programming matrix that do not influence the final result has been eliminated. The time complexity runtime on an input of \eqn{n}{n} points is improved to \eqn{O( (k-2) n^2 + n \lg n)}{O((k-2) n^2 + n lg(n))}.

The space complexity is \eqn{O(kn)}{O(kn)}.
}

\details{
\tabular{ll}{
Package: \tab Ckmeans.1d.dp\cr
Type: \tab Package\cr
Version: \tab 3.4.6-1\cr
Initial version: \tab 1.0\cr
Initial date: \tab 2010-10-26\cr
License: \tab LGPL (>= 3) \cr
}

}
\seealso{
  The \code{stats::kmeans} function that implements various heuristic \eqn{k}{k}-means algorithms.
}

\author{
Joe Song and Haizhou Wang
}

\references{
Wang, H. and Song, M. (2011) Ckmeans.1d.dp: optimal \var{k}-means clustering in one dimension by dynamic programming. \emph{The R Journal} \bold{3}(2), 29--33. Retrieved from \url{http://journal.r-project.org/archive/2011-2/RJournal_2011-2_Wang+Song.pdf}

Bellman, R. (1973) A note on cluster analysis and dynamic programming. \emph{Mathematical Biosciences} \bold{18}(3), 311--312.
}

\keyword{ package }


