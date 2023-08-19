\name{Ckmeans.1d.dp-package}
\alias{Ckmeans.1d.dp-package}
\docType{package}
\title{
Optimal, Fast, and Reproducible Univariate Clustering
}
\description{
This package provides a powerful set of tools for univariate data analysis with guaranteed optimality, efficiency, and reproducibility.
Four problems including univariate \eqn{k}{k}-means, \eqn{k}{k}-median, \eqn{k}{k}-segments, and multi-channel weighted \eqn{k}{k}-means are solved with guaranteed optimality and reproducibility. The core algorithm minimizes the (weighted) sum of within-cluster distances using respective metrics. Its advantage over heuristic clustering in efficiency and accuracy is increasingly pronounced as the number of clusters \eqn{k}{k} increases. Weighted \eqn{k}{k}-means can also optimally segment time series to perform peak calling. An auxiliary function generates histograms that are adaptive to patterns in data.

The Ckmeans.1d.dp algorithm clusters (weighted) univariate data given by a numeric vector \eqn{x}{x} into \eqn{k}{k} groups by dynamic programming \insertCite{Wang2011Ckmeans,song2020wuc}{Ckmeans.1d.dp}. It guarantees the optimality of clustering---the total of within-cluster sums of squares is always the minimum given the number of clusters \eqn{k}{k}. In contrast, heuristic univariate clustering algorithms may be non-optimal or inconsistent from run to run. As non-negative weights are supported, the algorithm can partition a time course using time points as input and corresponding values as weight. Based on an optimal clustering, a function can generate histograms adaptive to patterns in data. The implementation of this algorithm is numerically stable.

A linear time solution. Excluding the time for sorting \eqn{x}{x}, the weighted univariate clustering algorithm takes a runtime of \eqn{O(kn)}{O(kn)} \insertCite{song2020wuc}{Ckmeans.1d.dp}, linear in both sample size \eqn{n}{n} and the number of clusters \eqn{k}{k}. The approach was first introduced into version 3.4.9 (not publicly released) on July 16, 2016, using a new divide-and-conquer strategy integrating a previous theoretical result on matrix search \insertCite{aggarwal1987geometric}{Ckmeans.1d.dp} and a novel in-place search space reduction method \insertCite{song2020wuc}{Ckmeans.1d.dp}.

A log-linear time solution. Since version 3.4.6, a divide-and-conquer strategy that is simple to code reduces the runtime from \eqn{O(kn^2)}{O(kn^2)} down to \eqn{O(kn\lg n)}{O(kn lg n)} \insertCite{song2020wuc}{Ckmeans.1d.dp}.

A quadratic time solution. Before version 3.4.6, Ckmeans.1d.dp uses an algorithm that runs in quadratic runtime \eqn{O(kn^2)}{O(kn^2)} \insertCite{Wang2011Ckmeans}{Ckmeans.1d.dp}.

A cubic time solution. \insertCite{Bellman1973;textual}{Ckmeans.1d.dp} first described a general dynamic programming strategy for solving univariate clustering problems with additive optimality measures. The strategy, however, did not address any specific characteristics of the \eqn{k}{k}-means problem and its implied general algorithm will have a time complexity of \eqn{O(kn^3)}{O(kn^3)} on an input of \eqn{n}{n} points.

The space complexity in all cases is \eqn{O(kn)}{O(kn)}.

This package provides a powerful alternative to heuristic clustering and also new functionality for weighted clustering, segmentation, and peak calling with guaranteed optimality. It is practical for Ckmeans.1d.dp to find a few clusters on millions of sample points with optional weights in seconds using a single core on a typical desktop computer.

Third parties have ported various versions of this package to JavaScript, MATLAB, Python, Ruby, SAS, and Scala.
}

\details{
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
Joe Song, Hua Zhong, and Haizhou Wang
}

\references{

  \insertAllCited{}

% Bellman, R. (1973) A note on cluster analysis and dynamic programming. \emph{Mathematical Biosciences} \bold{18}(3), 311--312.

% Aggarwal, A., Klawe, M. M., Moran, S., Shor, P., and Wilber, R. (1987) Geometric applications of a matrix-searching algorithm. \emph{Algorithmica} \bold{2}(1-4), 195--208.

% Wang, H. and Song, M. (2011) Ckmeans.1d.dp: optimal \var{k}-means clustering in one dimension by dynamic programming. \emph{The R Journal} \bold{3}(2), 29--33. Retrieved from \url{https://journal.r-project.org/archive/2011-2/RJournal_2011-2_Wang+Song.pdf}
}

\keyword{ package }
\keyword{ cluster }
\keyword{ univar }




