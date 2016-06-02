/*
 Ckmeans_1d_dp.cpp -- Performs 1-D k-means by a dynamic programming
 approach that is guaranteed to be optimal.

 Joe Song
 Computer Science Department
 New Mexico State University
 joemsong@cs.nmsu.edu

 Haizhou Wang
 Computer Science Department
 New Mexico State University
 hwang@cs.nmsu.edu

 Created: May 19, 2007
 Updated: September 3, 2009
 Updated: September 6, 2009.  Handle special cases when K is too big or array
 contains identical elements.  Added number of clusters selection by the
 MCLUST package.
 Updated: Oct 6, 2009 by Haizhou Wang Convert this program from R code into
 C++ code.
 Updated: Oct 20, 2009 by Haizhou Wang, add interface to this function, so that
 it could be called directly in R

 Updated: March 20, 2014. Joe Song.
 1.  When the number of clusters is not uniquely specified, the code now
 automatically select an optimal number of clusters by Bayesian
 information criterion.
 2.  Reduced unnecessary sorting performed on the input string in
 kmeans_1d_dp().

 Updated: February 8, 2015.
 1.  Cleaned up the code by removing commented sections (Haizhou Wang)
 2.  Speed up the code slightly as suggested by a user (Joe Song)
 3.  Throw exceptions when for fatal errors (Joe Song)

 Updated: May 3, 2016
 1. Changed from 1-based to 0-based C implementation (MS)
 2. Optimized the code by reducing overhead. See 22% reduction in runtime to
 cluster one million points into two clusters. (MS)
 3. Removed the result class ClusterResult

 Updated: May 7, 2016
 1. Substantial runtime reduction. Added code to check for an upper bound
 for the sum of within cluster square distances. This reduced the runtime
 by half when clustering 100000 points (from standard normal distribution)
 into 10 clusters.
 2. Eliminated the unnecessary calculation of (n-1) elements in the dynamic
 programming matrix that are not needed for the final result. This
 resulted in enormous reduction in run time when the number of cluster
 is 2: assigning one million points into two clusters took half a
 a second on iMac with 2.93 GHz Intel Core i7 processor.

 Updated: May 9, 2016
 1. Added an alternative way to fill in the dynamic programming
 matix by i (columns in the matrix) to achieve speedup

 Updated: May 21, 2016
 1. Moved the weighted solutions to new files. They will be updated separately
 */

#include "Ckmeans.1d.dp.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// #define DEBUG

inline double ssq(const size_t j, const size_t i,
                  const std::vector<double> & sum_x,
                  const std::vector<double> & sum_x_sq)
{
  double sji;
  if(j > 0) {
    double muji = (sum_x[i] - sum_x[j-1]) / (i - j + 1);
    sji = sum_x_sq[i] - sum_x_sq[j-1] - (i - j + 1) * muji * muji;
  } else {
    sji = sum_x_sq[i] - sum_x[i] * sum_x[i] / (i+1);
  }
  sji = (sji < 0) ? 0 : sji;
  return sji;
}

void fill_row_k(int imin, int imax, int k,
            std::vector<std::vector<double>> & S,
            std::vector<std::vector<size_t>> & J,
            const std::vector<double> & sum_x,
            const std::vector<double> & sum_x_sq)
{
  if(imin > imax) {
    return;
  }

  const int N = S[0].size();

  int i = (imin + imax) / 2;

#ifdef DEBUG
  std::cout << std::endl << "  i=" << i << ":" << std::endl;
#endif

  // Initialization of S[k][i]:
  S[k][i] = S[k - 1][i - 1];
  J[k][i] = i;

  int jlow=k; // the lower end for j

  if(imin > k) {
    jlow = std::max(jlow, (int)J[k][imin-1]);
  }
  jlow = std::max(jlow, (int)J[k-1][i]);

  int jhigh = i - 1; // the upper end for j
  if(imax < N-1) {
    jhigh = std::min(jhigh, (int)J[k][imax+1]);
  }

#ifdef DEBUG
  std::cout << "    j-=" << jlow << ", j+=" << jhigh << ":" << std::endl;
#endif

  for(int j=jhigh; j>=jlow; --j) {

    // compute s(j,i)
    double sji = ssq(j, i, sum_x, sum_x_sq);

    // MS May 11, 2016 Added:
    if(sji + S[k-1][jlow-1] >= S[k][i]) break;

    // Examine the lower bound of the cluster border
    // compute s(jlow, i)
    double sjlowi = ssq(jlow, i, sum_x, sum_x_sq);

    double SSQ_jlow = sjlowi + S[k-1][jlow-1];

    if(SSQ_jlow < S[k][i]) {
      // shrink the lower bound
      S[k][i] = SSQ_jlow;
      J[k][i] = jlow;
    }
    jlow ++;

    double SSQ_j = sji + S[k - 1][j - 1];
    if(SSQ_j < S[k][i]) {
      S[k][i] = SSQ_j;
      J[k][i] = j;
    }
  }

#ifdef DEBUG
  std::cout << std::endl << // " k=" << k << ": " <<
    "\t" << S[k][i] << "\t" << J[k][i];
#endif
#ifdef DEBUG
  std::cout << std::endl;
#endif

  fill_row_k(imin, i-1, k, S, J, sum_x, sum_x_sq);
  fill_row_k(i+1, imax, k, S, J, sum_x, sum_x_sq);

}

void fill_dp_matrix(const std::vector<double> & x,
                    std::vector< std::vector< double > > & S,
                    std::vector< std::vector< size_t > > & J)
  /*
   x: One dimension vector to be clustered, must be sorted (in any order).
   S: K x N matrix. S[k][i] is the sum of squares of the distance from
   each x[i] to its cluster mean when there are exactly x[i] is the
   last point in cluster k
   J: K x N backtrack matrix

   NOTE: All vector indices in this program start at position 0
   */
{
  const int K = S.size();
  const int N = S[0].size();

  std::vector<double> sum_x(N), sum_x_sq(N);

  for(int i = 0; i < N; ++i) {
    if(i == 0) {
      sum_x[0] = x[0];
      sum_x_sq[0] = x[0] * x[0];
    } else {
      sum_x[i] = sum_x[i-1] + x[i];
      sum_x_sq[i] = sum_x_sq[i-1] + x[i] * x[i];
    }

    // Initialize for k = 0
    S[0][i] = ssq(0, i, sum_x, sum_x_sq);
    J[0][i] = 0;
  }

  for(int k = 1; k < K; ++k) {
    int imin;
    if(k < K - 1) {
      imin = std::max((size_t)1, (size_t)k);
    } else {
      // No need to compute S[K-1][0] ... S[K-1][N-2]
      imin = N-1;
    }

#ifdef DEBUG
    std::cout << std::endl << "k=" << k << ":";
#endif
    fill_row_k(imin, N-1, k, S, J, sum_x, sum_x_sq);
  }
}

void backtrack(const std::vector<double> & x,
               const std::vector< std::vector< size_t > > & J,
               int* cluster, double* centers, double* withinss, int* count)
{
  const size_t K = J.size();
  const size_t N = J[0].size();
  size_t cluster_right = N-1;
  size_t cluster_left;

  // Backtrack the clusters from the dynamic programming matrix
  for(int k = K-1; k >= 0; --k) {
    cluster_left = J[k][cluster_right];

    for(size_t i = cluster_left; i <= cluster_right; ++i)
      cluster[i] = k;

    double sum = 0.0;

    for(size_t i = cluster_left; i <= cluster_right; ++i)
      sum += x[i];

    centers[k] = sum / (cluster_right-cluster_left+1);

    for(size_t i = cluster_left; i <= cluster_right; ++i)
      withinss[k] += (x[i] - centers[k]) * (x[i] - centers[k]);

    count[k] = cluster_right - cluster_left + 1;

    if(k > 0) {
      cluster_right = cluster_left - 1;
    }
  }
}

void backtrack(const std::vector<double> & x,
               const std::vector< std::vector< size_t > > & J,
               std::vector<size_t> & count)
{
  const size_t K = J.size();
  const size_t N = J[0].size();
  size_t cluster_right = N-1;
  size_t cluster_left;

  // Backtrack the clusters from the dynamic programming matrix
  for(int k = K-1; k >= 0; --k) {
    cluster_left = J[k][cluster_right];
    count[k] = cluster_right - cluster_left + 1;
    if(k > 0) {
      cluster_right = cluster_left - 1;
    }
  }
}

template <class ForwardIterator>
size_t numberOfUnique(ForwardIterator first, ForwardIterator last)
{
  size_t nUnique;

  if (first == last) {
    nUnique = 0;
  } else {
    nUnique = 1;
    for (ForwardIterator itr=first+1; itr!=last; ++itr) {
      if (*itr != *(itr -1)) {
        nUnique ++;
      }
    }
  }
  return nUnique;
}

void kmeans_1d_dp(const double *x, const size_t N, const double *y,
                  size_t Kmin, size_t Kmax,
                  int* cluster, double* centers,
                  double* withinss, int* size)
{
  // Input:
  //  x -- an array of double precision numbers, not necessarily sorted
  //  Kmin -- the minimum number of clusters expected
  //  Kmax -- the maximum number of clusters expected
  // NOTE: All vectors in this program is considered starting at position 0.

  std::vector<double> x_sorted(N);

  std::vector<double> y_sorted;
  auto is_equally_weighted(true);

  std::vector<size_t> order(N);
  std::size_t n(0);
  std::generate(order.begin(), order.end(), [&]{ return n++; });

  //  sort dimension k of X in increasing order
  std::sort(order.begin(), order.end(),
            [&](size_t i1, size_t i2) { return x[i1] < x[i2]; } );

  for(auto i=0; i<order.size(); ++i) {
    x_sorted[i] = x[order[i]];
  }

  // check to see if unequal weight is provided
  if(y != NULL) {
    is_equally_weighted = true;
    for(auto i=1; i<N; ++i) {
      if(y[i] != y[i-1]) {
        is_equally_weighted = false;
        break;
      }
    }
  }

  if(! is_equally_weighted) {
    y_sorted.resize(N);
    for(auto i=0; i<order.size(); ++i) {
      y_sorted[i] = y[order[i]];
    }
  }

  const size_t nUnique = numberOfUnique(x_sorted.begin(), x_sorted.end());

  Kmax = nUnique < Kmax ? nUnique : Kmax;

  if(nUnique > 1) { // The case when not all elements are equal.

    std::vector< std::vector< double > > S( Kmax, std::vector<double>(N) );
    std::vector< std::vector< size_t > > J( Kmax, std::vector<size_t>(N) );

    size_t Kopt;

    // Fill in dynamic programming matrix
    if(is_equally_weighted) {

      fill_dp_matrix(x_sorted, S, J);

      // Choose an optimal number of levels between Kmin and Kmax
      Kopt = select_levels(x_sorted, J, Kmin, Kmax);

    } else {
      fill_weighted_dp_matrix(x_sorted, y_sorted, S, J);

      // Choose an optimal number of levels between Kmin and Kmax
      Kopt = select_levels_weighted(x_sorted, y_sorted, J, Kmin, Kmax);
    }


    if (Kopt < Kmax) { // Reform the dynamic programming matrix S and J
      J.erase(J.begin() + Kopt, J.end());
    }

    std::vector<int> cluster_sorted(N);

    // Backtrack to find the clusters beginning and ending indices
    if(is_equally_weighted) {
      backtrack(x_sorted, J, &cluster_sorted[0], centers, withinss, size);
    } else {
      backtrack_weighted(x_sorted, y_sorted, J, &cluster_sorted[0], centers, withinss, size);
    }

    for(size_t i = 0; i < N; ++i) {
      // Obtain clustering on data in the original order
      cluster[order[i]] = cluster_sorted[i];
    }

  } else {  // A single cluster that contains all elements

    for(size_t i=0; i<N; ++i) {
      cluster[i] = 0;
    }

    centers[0] = x[0];
    withinss[0] = 0.0;
    size[0] = N * (is_equally_weighted ? 1 : y[0]);
  }
}  //end of kmeans_1d_dp()
