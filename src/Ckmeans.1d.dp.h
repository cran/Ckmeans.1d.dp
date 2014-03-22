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

/* data that return by kmeans.1d.dp()*/
class ClusterResult {
public:
    size_t nClusters;
    std::vector<size_t> cluster;  	/*record which cluster each point belongs to*/
    std::vector<double> centers;	/*record the center of each cluster*/
    std::vector<double> withinss;   /*within sum of distance square of each cluster*/
    std::vector<size_t> size;		/*size of each cluster*/
};

/* One-dimensional cluster algorithm implemented in C++ */
/*x is input one-dimensional vector and Kmin and Kmax stand for the range 
    for the number of clusters*/
ClusterResult
kmeans_1d_dp(const std::vector<double> & x, size_t Kmin, size_t Kmax);
