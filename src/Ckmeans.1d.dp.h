/*
Ckmeans_1d_dp.h --- Head file for Ckmeans.1d.dp
                    Declare a class "data" and
					wrap function "kmeans_1d_dp()"

  Haizhou Wang
  Computer Science Department
  New Mexico State University
  hwang@cs.nmsu.edu

Created: Oct 10, 2010
*/

#include <cstddef> // For size_t
#include <vector>

/* One-dimensional cluster algorithm implemented in C++ */
/*x is input one-dimensional vector and Kmin and Kmax stand for the range
    for the number of clusters*/
void kmeans_1d_dp(const double *x, const size_t N, size_t Kmin, size_t Kmax,
                  int* cluster, double* centers, double* withinss, int* size);
