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
*/

#include "Ckmeans.1d.dp.h"
#include <vector>   
#include "R.h" 		// R memory io
#include "Rmath.h" 	// R math functions

using namespace std;

/*Wrapper function to call kmeans_1d_dp()*/
extern "C" {
	void Ckmeans_1d_dp(double *x, int* s, int* Ks, int *nK, int* cluster,
                       double* centers, double* withinss, int* size)
	{
		size_t length = (size_t) s[0];
        
		vector<double> input(length+1);
		
		for(size_t i=1; i<input.size(); i++) {
			input[i] = x[i-1];
        }

        vector<size_t> levels(*nK);
        for (size_t k=0; k<*nK; k++) {
            levels[k] = Ks[k];
        }

		ClusterResult result;  // Clustering result
		
		// Call C++ version one-dimensional clustering algorithm*/
        result = kmeans_1d_dp( input, levels );
		
        // Since R doesn't allow return value from C/C++ function,
        // we use input pointers to pass the clustering result back.
		for(size_t i=1; i <= length; i++)
			cluster[i-1] = (int) result.cluster[i];
		for(size_t k=1; k <= result.nClusters; k++)
			centers[k-1] = result.centers[k];
		for(size_t k=1; k <= result.nClusters; k++)
			withinss[k-1] = result.withinss[k];
		for(size_t k=1; k <= result.nClusters; k++)
			size[k-1] = (int) result.size[k];
	}
}
