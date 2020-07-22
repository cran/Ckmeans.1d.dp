# NEWS

## Version 4.3.3

  2020-07-21

  1. Created version 4.3.3 from 4.3.2.1

## Version 4.3.2.1 (not publicly released)

  2020-07-18
  
  1. Updated REFERENCES.bib, CITATION.
  2. Updated README.md text and badges.
  3. Updated user manuals.
  
  2020-03-15
  
  1. Created version 4.3.2.1 from 4.3.2
  2. Updated REFERENCES.bib, CITATION, and manual Rd files.

## Version 4.3.2
  2020-03-13
  
  1. Changed the random data from uniform integer to standard normal
  to avoid examples with multile optimal solutions occuring to
  integers in test_MC_WUC.R.
  
  2020-01-20
  
  1. Updated a vignette to illustrate how to find boundaires between
  consecutive clusters.
  
  2020-01-03
  
  1. Internally, a specialized version for unweighted Euclidean (L2)
  distance based univariate clustering is added. In previous versions, 
  unweighted and weighted multiple metric clustering algorithms were
  implemented in a unified code framework, which is good for software
  engineering but carries unnecessary overhead. The new specialized
  version can speed up the unweighted L2 algorithm by about 20%. This
  is perhaps the most popular task, thus benefiting most users. There
  is no change in the user interface.
  
  2020-01-02
  
  1. Created version 4.3.2 from 4.3.0
  2. Revised CITATION file
  3. Revised Readme.md file
  
## Version 4.3.1 (not publicly released)
  
  2019-12-08
  
  1. Version 4.3.1 was created.

## Version 4.3.0

  2019-09-06
  
  1. Updated package documentation.
  2. Introduced NEWS.md instead the plain text NEWS.
  3. Changed the package title from **Optimal and Fast Univariate Clustering** to **Optimal, Fast, and Reproducible Univariate Clustering**
  4. Created README.md to introduce the package.
  
  2019-09-02
  
  1. Added the optimal multi-channel weighted univariate
  clustering function, called "MultiChannel.WUC" in short for now.
  Added related R document and testthat cases. The example in the
  "MultiChannel.WUC" R document illustrates how to run the function.
  2. Added source files: MCW_main.cpp, MCW_functions.cpp,
  MCW_functions.h, MCW_fill_SMAWK.cpp
  3. Added R file: MultiChannel.WUC.R
  4. Added Rd file: MultiChannel.WUC.Rd
  5. Added testthat file: test_MC_WUC.R
  6. Modified DESCRIPTION file.
  7. Added two imports in the NAMESPACE file.

## Version 4.2.3 (not publicly released)
  
  2018-09-26
  
  1. Removed unnecessary version requirement for Rcpp introduced in
  the preivous version.

  2018-09-24
  
  1. Version created to remove typos in DESCRIPTION.
  2. Updated vignette on weight scaling.

## Version 4.2.2
  
  2018-09-21
  
  1. Modified the package to use Rcpp interface instead of the
     old-style C interface.

## Version 4.2.1
  
  2017-06-10
  
  1. Added a new vignette "Tutorial: Linear weight scaling in cluster analysis".
  2. Re-organized manuals and updated documentation.

## Version 4.2.0
  
  2017-05-29
  
  1. Increased log likelihood calculation to long double
     precision in C++ function weighted_select_levels.cpp.
  2. Replaced std::accumulate() function by for-loop addition in
     C++ function weighted_select_levels.cpp. This resolved a
     numerical overflow issue when the weight values are large.
  3. Now R function plotBIC() automatically adjusts the "k*="
     text position, so that the text label is placed entirely
     within the BIC curve area and will not extend into the
     figure margin.
  4. The cluster size has been changed from integer to double to
     accomodate weighted cluster size in both R and C++ code.
  5. Force any weight vector to be equal weight in new R
     function Ckmedian.1d.dp.
  6. Introduced S3 methods print and plot for Ckmedian.1d.dp
     and Cksegs.1d.dp objects.

  2017-03-18
  
  1. Introduced Cksegs.1d.dp() function for k-segments clustering
     of y with or without x. Only method="quadratic" guarantees
     optimality.
  2. Expanded k-median clustering to work with all possible
     methods. Only unweighted solution guarantees optimality.

## Version 4.1.0
  
  2017-03-02
  
  1. Introduced function Ckmedian.1d.dp() for k-median unweighted
     clustering.

## Version 4.0.2

  2017-02-16
  
  1. Fixed symbol encoding used in NEWS.
  2. Updated documentation.

## Version 4.0.1
  
  2017-02-16
  
  1. Fixed a warning message in the use of 'R_registerRoutines' and
     'R_useDynamicSymbols'.
  2. Fixed a memory leak issue: invalid read of size 8.

## Version 4.0.0
  
  2017-02-11
  
  1. Removed some examples for future use.

## Version 3.4.15

  1. Minor changes.

## Version 3.4.14
  
  2017-01-03
  
  1. Minor changes in documentation files.

  2017-01-02
  
  1. Changed package title to "Optimal and Fast Univariate Clustering".

  2016-12-27
  
  1. When the input vector x is empty, function Ckmeans.1d.dp now generates a
     warning message instead of stops on error. Ckmeans.1d.dp.R is modified.
  2. Print out appropriate warning messages when input x does not have an
     appropriate type. ahist.R is modified.

## Version 3.4.13
  
  2016-12-19
  
  1. Revised the comparison function in sorting so that the code can be
     compiled by C++98, as requested by a user.

  2016-12-06
  
  1. Expanded the ahist() function to support weighted adaptive histogram

  2016-12-04
  
  1. Expanded the vignette of adaptive histograms to a tutorial.
  2. Expanded the vignette of optimal univariate k-means clustering to a tutorial.
  3. Update the time course example in Ckmeans.1d.dp function

  2016-10-22
  
  1. Added a vignette to visualize examples of adaptive histograms.
  2. Added a vignette to visualize examples of optimal univariate k-means
  clustering.

  2016-10-21
  
  1. Added an equal bin width histogram example to contrast with the adaptive
  histogram.

  2016-10-16
  
  1. Moved ahist() function from visualize.R to a new R file ahist.R.

  2016-10-15
  
  1. Added a breaks argument to ahist() so as use default graphics::hist() but
     with the capacity to add sticks to the histogram generated.
  2. Added a skip.empty.bin.color argument to ahist() to gain more control over
     colors of the histogram bars.

  2016-10-12
  
  1. Added a data argument to ahist() to provide raw data for visualization.
  2. Allow x to ahist() to be an object of the class "Ckmeans.1d.dp" to avoid
     recomputing the clustering if it has already been done. This requires the
     data for clustering to be provided via the data argument.

  2016-10-11
  
  1. Added an argument add.sticks=TRUE to ahist() to turn on or off the sticks
  just above the horizontal axis.
  2. Added an argument style to ahist() for different styles of adaptive
  histogram.

  2016-10-01
  
  1. Added a new function plot() to visualize the clusters.
  2. Added a new function plotBIC() to show the Bayesian information criterion
  as a function of number of clusters.

  2016-09-30
  
  1. Updated examples of ahist().
  2. Added sticks to ahist() to show the original input data.

  2016-09-27
  
  1. Fixed ahist() when there is only a single bin detected.
  2. Made ahist() run much faster than the previous version.
  3. Updated previous examples and added more examples to illustrate
     the use of ahist() better.

  2016-09-25
  
  1. Introduced a new function ahist() to generate adaptive histograms
     corresponding to the optimal univariate k-means clustering.

  2016-09-24
  
  1. Known issue: loglinear option may generate optimal clustering different
     from linear and quadratic.
  2. The default k estimation method is updated. Updated number of cluster k
     estimation. The main difference is when there are duplicates in the data.
     Otherwise, the estimated k would be the same with previous versions. Added
     an argument estimate.k in fiction Ckmeans.1d.dp() to use the BIC method in
     version 3.4.12 or earlier to estimated k for compatibility.

## Version 3.4.12
  
  2016-08-20
  
  1. The weighted univariate k-means now runs in $O(kn)$, down from $O(kn^2)$.
     This is a result of integrating weighted and unweighted k-means
     clustering into a unified dynamic programming function without sacrificing
     performance. This also fixed a bug in the previous loglinear-time weighted
     k-means implementation.

## Version 3.4.9
  
  2016-07-19
  
  1. If the input array is already sorted, sorting is not performed again.
  2. Added an option method to select either the linear or loglinear algorithm.

  2016-07-16
  
  1. Implemented linear recursive algorithm based on the method described in
  (Aggarwal et al., 1987)

## Version 3.4.8
  
  2016-06-01
  
  1. Implemented an iterative O(nlgn+kn) algorithm. This version completely
     eliminates the involved divide-and-conquer strategy reported in the
     literature and further reduced the overhead.

     This implementation was later determined to be incorrect.

## Version 3.4.7
  
  2016-05-30
  
  1. Implemented an O(nlgn+kn) algorithm combining divide-and-conquer and
     dynamic programming. The space is still O(kn). The runtime is now practical
     for very large sample sizes for any number of clusters.

     This implementation was later determined to be incorrect.

## Version 3.4.6
  
  2016-05-25
  
  1. Implemented an O(kn lg n) algorithm, speeding up the program greatly.

## Version 3.4.5
  
  2016-05-22
  
  1. $s[j,i]$ is now computed in constant time based on pre-computed
     sums of input x and its squares from 0 to i.

## Version 3.4.4
  
  2016-05-22
  
  1. Incorporated a numerically stable method for computing sample variance when
     selecting the number of clusters.
  2. Improved documentation.
  3. Removed a typo in describing time complexity.

  2016-05-18
  
  1. Now Ckmeans.1d.dp() function returns "totss", "tot.withinss", and
     "betweenss" statistics to summarize the optimal clustering obtained.
  2. print.Ckmeans.1d.dp() print out the above statistics.

## Version 3.4.3
  
  2016-05-15
  
  1. Upgraded to support c++11
  2. Introduced optimal k-means clustering for weighted data

## Version 3.4.2
  
  2016-05-14
  
  1. Implemented backward filling of the dynamic programming matrix to
     utilize lower bounds for the optimal cluster boundary. This step
     substantially reduced the runtime by half (two or more times faster).

  2016-05-07
  
  1. Implemented mathematically proven tighter ranges when searching for
     cluster boundaries. The runtime of the function is greatly reduced.
     Most notably, the runtime is roughly constant when number of clusters
     increases after k=2.
  2. Integrated all test cases into one single file.

## Version 3.4.0
  
  2016-05-07
  
  1. Substantial runtime reduction. Added code to check for an upper bound
     for the sum of within cluster square distances. This reduced the runtime
     by half when clustering 100000 points (from standard normal distribution)
     into 10 clusters.
  2. Eliminated the unnecessary calculation of (n-1) elements in the dynamic
     programming matrix that are not needed for the final result. This
     resulted in enormous reduction in run time when the number of cluster
     is 2: assigning one million points into two clusters took half a
     a second on iMac with 2.93 GHz Intel Core i7 processor.
  3. Included a reference to the first description of the dynamic programming
     solution by Richard Bellman (1973).

## Version 3.3.3
  
  2016-05-03
  
  1. Fixed a bug on cluster assignment when there is only one cluster. This
     was a bug introduced in version 3.3.2.

## Version 3.3.2
  
  2016-05-03
  
  1. Added automatic test cases.
  2. Removed an incorrect warning message when the number of clusters is equal
     to the number of unique elements in the input vector.
  3. Changed from 1-based to 0-based C implementation.
  4. Optimized the code by reducing overhead. See 22% reduction in runtime to
     repeatedly cluster seven points into two clusters one million times.

## Version 3.3.1 
  
  2015-02-10
  
  1. Fixed a problem that prevented Windows compilation (now forced the size_t
     type to unsigned long in max() function.

## Version 3.3.0 
  
  2015-02-09

  1. Added automated test cases into the package.
  2. Changed the code to not issue a warning message when the number of clusters
     is estimated to be 1.
  3. When lower bound of the number of clusters is greater than the unique
     number of elements in the input vector, both the min and max numbers of
     clusters are set to the number of unique number of input values.
  4. When the upper bound of the number of clusters is greater than the unique
     number of elements in the input vector, the max number of clusters is set
     to the number of unique elements in the input vector.
  5. Use warning() instead of cat() to display warning messages.
  6. Incorporate changes suggested by a user to speed up the code.
  7. Revised the examples and documentation to improve usability of the package
     in general.
  8. Started the NEWS file.

## Version 3.02 

   2014-03-24 and earlier

  1. The program now automatically determines the number of clusters from a
     given range.
  2. The code is optimized for further speedup.
