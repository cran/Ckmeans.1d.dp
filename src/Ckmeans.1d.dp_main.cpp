/*
 Ckmeans_1d_dp_main.cpp --- wrapper function for "kmeans_1d_dp()"

 Created: Oct 10, 2010

 Haizhou Wang
 Computer Science Department
 New Mexico State University
 hwang@cs.nmsu.edu

 Modified:
 March 20, 2014. Joe Song. Removed parameter int *k from the function.
 Added "int *Ks" and "int *nK" to provide a range of the number of clusters
 to search for. Made other changes.
 March 29, 2014. Haizhou Wang. Replaced "int *Ks" and "int *nK" by
 "int *minK" and "int *maxK".
 */

#include "Ckmeans.1d.dp.h"

/*Wrapper function to call kmeans_1d_dp()*/
extern "C" {
  /*
   x: An array containing input data to be clustered.
   length: length of the one dimensional array.
   minK: Minimum number of clusters.
   maxK: Maximum number of clusters.
   cluster:  An array of cluster IDs for each point in x.
   centers:  An array of centers for each cluster.
   withinss: An array of within-cluster sum of squares for each cluster.
   size:     An array of sizes of each cluster.
   */

  void Ckmeans_1d_dp(double *x, int* length, double *y, int * ylength,
                     int* minK, int *maxK, int* cluster,
                     double* centers, double* withinss, int* size)
  {
    // Call C++ version one-dimensional clustering algorithm*/
    if(*ylength != *length) { y = 0; }

    kmeans_1d_dp(x, (size_t)*length, y, (size_t)(*minK), (size_t)(*maxK),
                 cluster, centers, withinss, size);

    // Change the cluster numbering from 0-based to 1-based
    for(size_t i=0; i< *length; ++i) {
      cluster[i] ++;
    }
  }
} // End of extern "C"
