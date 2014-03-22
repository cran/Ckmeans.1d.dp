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
#include <vector>   
#include "R.h" 		// R memory io
#include "Rmath.h" 	// R math functions

/*Wrapper function to call kmeans_1d_dp()*/
extern "C" {
    /*
     x: One dimensional array containing input data to be clustered.
     s: Size of the one dimensional array.
     minK: Minimum cluster level.
     maxK: Maximum cluster level.
     cluster:  A list of cluster IDs for each point in x.
     centers:  A list of centers for each cluster.
     withinss: A list of within-cluster sum of squares for each cluster.
     size:     A list of sizes of each cluster.
     */
    
	void Ckmeans_1d_dp(double *x, int* s, int* minK, int *maxK, int* cluster,
                       double* centers, double* withinss, int* size)
	{
		size_t length = (size_t) *s;
        
		std::vector<double> input(length+1);
		
		for(size_t i=1; i<input.size(); ++i) {
			input[i] = x[i-1];
        }

		ClusterResult result;  // Clustering result
		
		// Call C++ version one-dimensional clustering algorithm*/
        result = kmeans_1d_dp( input, (size_t)(*minK), (size_t)(*maxK) );
		
        // Since R doesn't allow return value from C/C++ function,
        // we use input pointers to pass the clustering result back.
		for(size_t i=1; i <= length; ++i) {
			cluster[i-1] = (int) result.cluster[i];
        }
        
		for(size_t k=1; k <= result.nClusters; ++k) {
			centers[k-1] = result.centers[k];
			withinss[k-1] = result.withinss[k];
			size[k-1] = (int) result.size[k];
        }
	}
} // End of extern "C"
