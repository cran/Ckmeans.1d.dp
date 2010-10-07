/*
Ckmeans_1d_dp_main.cpp --- wrap function for 
						   "kmeans_1d_dp()"
                  
  Haizhou Wang
  Computer Science Department
  New Mexico State University
  hwang@cs.nmsu.edu
				  
Created: Oct 10, 2010
*/

#include "Ckmeans.1d.dp.h"
#include <vector>   
#include "R.h" 		// R memory io
#include "Rmath.h" 	// R math functions

using namespace std;

/*Wrap function to call kmeans_1d_dp()*/
extern "C" {
	void Ckmeans_1d_dp( double *x, int* s, int* k, int* cluster, double* centers, double* withinss, int* size)
	{
		int length = s[0];
		int level = k[0];
		vector<double> input(length+1);
		
		for(int i=1;i<(length+1);i++)
			input[i] = x[i-1];
			
		data result;  /*object of class data, used to convey the final result*/
		
		/*Call C++ version one-dimensional clustering algorithm*/
		result = kmeans_1d_dp( input, level);

		/*Since R doesn't allow return value from C/C++ function, using pointers to give back result*/
		for(int i=1;i<(length+1);i++)
			cluster[i-1] = result.cluster[i];
		for(int i=1;i<(level+1);i++)
			centers[i-1] = result.centers[i];
		for(int i=1;i<(level+1);i++)
			withinss[i-1] = result.withinss[i];
		for(int i=1;i<(level+1);i++)
			size[i-1] = result.size[i];
	}
}
