\name{Ckmeans.1d.dp-package}
\alias{Ckmeans.1d.dp-package}
\docType{package}
\title{
Optimal \var{k}-Means Clustering in One-dimension by Dynamic Programming
}
\description{
The Ckmeans.1d.dp algorithm clusters 1-D data given by a numeric vector \eqn{x}{x} into \eqn{k}{k} groups by dynamic programming (Wang and Song, 2011). It guarantees the optimality of clustering -- the sum of within-cluster sum of squares is always the minimum. In contrast, heuristic \eqn{k}{k}-means algorithms may be non-optimal or inconsistent from run to run.

Since version 3.4.0, the runtime of this dynamic programming algorithm has been greatly reduced after two changes to trim down the search space. These improvements, not discussed in (Wang and Song, 2011), will be described in future publications. While runtime was at least two times faster for all tested examples, the most notable improvement occurred when \eqn{k=2}{k=2} or \eqn{k}{k} is very large. First, an upper bound for the sum of within cluster square distances is checked to reduce the dynamic programming search space. Second, the unnecessary calculation of \eqn{n-1}{n-1} elements in the dynamic programming matrix that do not influence the final result has been eliminated. The asymptotic runtime on \eqn{n}{n} points is improved to \eqn{O( (k-1) n^2 + n \lg n)}{O((k-1) n^2 + n lg(n))}. When \eqn{k=2}{k=2}, it is \eqn{O(n\lg n)}{O(n lg(n))}, leading to enormous reduction in runtime: assigning one million points into two clusters took half a second, down from about a week, on iMac with a 2.93 GHz Intel Core i7 processor. Overall, the function now runs fast for all tested input, and the improvements are especially pronounced for \eqn{k=2}{k=2} and large \eqn{k}{k}. The space complexity is \eqn{O(kn)}{O(kn)}.

Richard Bellman (1973) first described a general dynamic programming strategy for solving univariate clustering problems with additive optimality measures. The strategy, however, did not address any specific characteristics of the \eqn{k}{k}-means problem and its implied general algorithm will have a time complexity of \eqn{O(kn^3)}{O(k n^3)} and space complexity of \eqn{O(kn)}{O(kn)}.
}

\details{
\tabular{ll}{
Package: \tab Ckmeans.1d.dp\cr
Type: \tab Package\cr
Version: \tab 3.4.0-1\cr
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


