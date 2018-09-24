\name{Ckmeans.1d.dp-package}
\alias{Ckmeans.1d.dp-package}
\docType{package}
\title{
Optimal and Fast Univariate Clustering
}
\description{
This package provides fast optimal univariate clustering and segementation by dynamic programming. Three types of problem including univariate \eqn{k}{k}-means, \eqn{k}{k}-median, and \eqn{k}{k}-segments are solved with guaranteed optimality and reproducibility. The core algorithm minimizes the sum of within-cluster distances using respective metrics. Its advantage over heuristic clustering algorithms in efficiency and accuracy is increasingly pronounced as the number of clusters k increases. Weighted k-means and unweighted k-segments algorithms can also optimally segment time series and perform peak calling. An auxiliary function generates histograms that are adaptive to patterns in data. In contrast to heuristic methods, this package provides a powerful set of tools for univariate data analysis with guaranteed optimality.

The Ckmeans.1d.dp algorithm clusters univariate data given by a numeric vector \eqn{x}{x} into \eqn{k}{k} groups by dynamic programming (Wang and Song, 2011). It guarantees the optimality of clustering---the total of within-cluster sums of squares is always the minimum given the number of clusters \eqn{k}{k}. In contrast, heuristic univariate clustering algorithms may be non-optimal or inconsistent from run to run. As unequal non-negative weights are supported for each point, the algorithm can also segment a time course using the time points as input and the values at each time point as weight. Utilizing the optimal clusters, a function can generate histograms adaptive to patterns in data. This method is numerically stable.

Apart from the time for sorting \eqn{x}{x}, the algorithm, now covering both unweighted and weighted univariate clustering problems, takes a runtime of \eqn{O(kn)}{O(kn)}, linear in both sample size \eqn{n}{n} and the number of clusters \eqn{k}{k}. The new approach was first introduced into version 3.4.9 (not publicly released) on July 16, 2016, using a new divide-and-conquer strategy integrating a previous theoretical result on matrix search (Aggarwal et al., 1987) and a novel in-place search space reduction method. Earlier versions since 3.4.6 used a divide-and-conquer strategy in dynamic programming to reduce the runtime from \eqn{O(kn^2)}{O(kn^2)} down to \eqn{O(kn\lg n)}{O(kn lg n)}. Both fast algorithms take runtime comparable to heuristic \eqn{k}{k}-means.

Richard Bellman (1973) first described a general dynamic programming strategy for solving univariate clustering problems with additive optimality measures. The strategy, however, did not address any specific characteristics of the \eqn{k}{k}-means problem and its implied general algorithm will have a time complexity of \eqn{O(kn^3)}{O(kn^3)} on an input of \eqn{n}{n} points.

The space complexity in all cases is \eqn{O(kn)}{O(kn)}.

This package provides a powerful alternative to heuristic clustering algorithms and also new functionality for weighted clustering, segmentation, and peak calling with guaranteed optimality. It is practical for Ckmeans.1d.dp to cluster millions of sample points with optional weights in seconds using a single processor on a not very recent desktop computer.

Third parties have ported various versions of this package to JavaScript, MATLAB, Python, Ruby, SAS, and Scala.
}

\details{
\tabular{ll}{
Package: \tab Ckmeans.1d.dp\cr
Type: \tab Package\cr
Version: \tab 4.2.2\cr
Initial release version: \tab 1.0\cr
Initial release date: \tab 2010-10-26\cr
License: \tab LGPL (>= 3)}

The most important function of this package is \code{\link{Ckmeans.1d.dp}} that implements optimal univariate clustering either weighted or unweighted. It also contains an adaptive histogram function \code{\link{ahist}}, plotting functions \code{\link{plot.Ckmeans.1d.dp}} and \code{\link{plotBIC}}, and a print function \code{\link{print.Ckmeans.1d.dp}}.
}

\seealso{
  The \code{\link{kmeans}} function in package \pkg{\link{stats}} that implements several heuristic \eqn{k}{k}-means algorithms.
}

\section{Disclaimer}{
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see \url{https://www.gnu.org/licenses/}.}

\author{
Joe Song and Haizhou Wang
}

\references{
Bellman, R. (1973) A note on cluster analysis and dynamic programming. \emph{Mathematical Biosciences} \bold{18}(3), 311--312.

Aggarwal, A., Klawe, M. M., Moran, S., Shor, P., and Wilber, R. (1987) Geometric applications of a matrix-searching algorithm. \emph{Algorithmica} \bold{2}(1-4), 195--208.

Wang, H. and Song, M. (2011) Ckmeans.1d.dp: optimal \var{k}-means clustering in one dimension by dynamic programming. \emph{The R Journal} \bold{3}(2), 29--33. Retrieved from \url{https://journal.r-project.org/archive/2011-2/RJournal_2011-2_Wang+Song.pdf}
}

\keyword{ package }
\keyword{ cluster }
\keyword{ univar }




