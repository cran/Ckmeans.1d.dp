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

#include <vector>
using namespace std;

/* data that return by kmeans.1d.dp()*/
class ClusterResult {
public:
    size_t nClusters;
    vector<size_t> cluster;  	/*record which cluster each point belongs to*/
    vector<double> centers;	/*record the center of each cluster*/
    vector<double> withinss;/*within sum of distance square of each cluster*/
    vector<size_t> size;		/*size of each cluster*/
};

/*one-dimensional cluster algorithm implemented in C*/
/*x is input one-dimensional vector and K stands for the cluster level*/
ClusterResult kmeans_1d_dp(const vector<double> & x, const vector<size_t> & Ks);
