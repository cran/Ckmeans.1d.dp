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
 */

#include "Ckmeans.1d.dp.h"

#include <algorithm>
#include <cmath>
#include <ctime>
#include <iostream>
#include <string>
#include <vector>

#ifndef M_PI
const double M_PI = 3.14159265359;
#endif


void fill_dp_matrix(const std::vector<double> & x,
                    std::vector< std::vector< double > > & D,
                    std::vector< std::vector< size_t > > & B)
/*
 x: One dimension vector to be clustered, must be sorted (in any order).
 D: Distance matrix
 B: Backtrack matrix
 
 NOTE: All vector indices in this program start at position 1,
       position 0 is not used.
 */
{
    const size_t K = D.size() - 1;
    const size_t N = D[0].size() - 1;
    
    for(size_t i = 1; i <= K; ++i) {
        D[i][1] = 0.0;
        B[i][1] = 1;
    }
    
    double mean_x1, mean_xj;
	double d; // d: is the sum of squared distances from x_j ,. . ., x_i to their mean
    
    for(size_t k = 1; k <= K; ++k) {
        mean_x1 = x[1];
        
        for(size_t i = std::max(2ul,(unsigned long)k); i <= N; ++i) { // for(size_t i = 2; i <= N; ++i) {
            if(k == 1) {
                D[1][i] = D[1][i-1] + (i-1) / (double) i *
                          (x[i] - mean_x1) * (x[i] - mean_x1);
                mean_x1 = ((i - 1) * mean_x1 + x[i]) / (double)i;
                
                B[1][i] = 1;
            } else {
                D[k][i] = -1.0;
                d = 0.0;
                mean_xj = 0.0;
                
                for(size_t j=i; j>=k; --j) { // for(size_t j = i; j >= 1; --j) {
                    d = d + (i - j) / (double) (i - j + 1) * 
						(x[j] - mean_xj) * (x[j] - mean_xj);
                    mean_xj = (x[j] + (i - j) * mean_xj) / (double)(i - j + 1);
                    
                    if(D[k][i] == -1.0) { //initialization of D[k,i]
                        if(j == 1) {
                            D[k][i] = d;
                            B[k][i] = j;
                        } else {
                            D[k][i] = d + D[k - 1][j - 1];
                            B[k][i] = j;
                        }
                    } else {
                        if(j == 1) {
                            if(d <= D[k][i]) {
                                D[k][i] = d;
                                B[k][i] = j;
                            }
                        } else {
                            if(d + D[k - 1][j - 1] < D[k][i]) {
                                D[k][i] = d + D[k - 1][j - 1];
                                B[k][i] = j;
                            }
                        }
                    }
                }
            }
        }
    }
}

void backtrack(const std::vector<double> & x,
               const std::vector< std::vector< size_t > > & B,
               ClusterResult & result)
{
    const size_t K = B.size() - 1;
    const size_t N = B[0].size() - 1;
    size_t cluster_right = N;
    size_t cluster_left;
    
    result.nClusters = K;
    
    result.cluster.resize(N + 1);
    result.centers.resize(K + 1);
    result.withinss.resize(K + 1);
    result.size.resize(K + 1);
    
    // Backtrack the clusters from the dynamic programming matrix
    for(size_t k = K; k >= 1; --k) {
        cluster_left = B[k][cluster_right];
        
        for(size_t i = cluster_left; i <= cluster_right; ++i)
            result.cluster[i] = k;
        
        double sum = 0.0;
        
        for(size_t i = cluster_left; i <= cluster_right; ++i)
            sum += x[i];
        
        result.centers[k] = sum/(cluster_right-cluster_left+1);
        
        for(size_t i = cluster_left; i <= cluster_right; ++i)
            result.withinss[k] += (x[i] - result.centers[k]) *
								  (x[i] - result.centers[k]);
        
        result.size[k] = cluster_right - cluster_left + 1;
        
        if(k > 1) {
            cluster_right = cluster_left - 1;
        }
    }
}

// Choose an optimal number of levels between Kmin and Kmax
size_t select_levels(const std::vector<double> & x,
                     const std::vector< std::vector< size_t > > & B,
                     size_t Kmin, size_t Kmax)
{
    if (Kmin == Kmax) {
        return Kmin;
    }
    
    const std::string method = "normal"; // "uniform" or "normal"
    
    size_t Kopt = Kmin;
    
    const size_t base = 1;  // The position of first element in x: 1 or 0.
    const size_t N = x.size() - base;
    
    double maxBIC;
    
    for(size_t K = Kmin; K <= Kmax; ++K) {
        
        std::vector< std::vector< size_t > > BK(B.begin(), B.begin()+K+1);

        ClusterResult result;
        // Backtrack the matrix to determine boundaries between the bins.
        backtrack(x, BK, result);

        size_t indexLeft = base;
        size_t indexRight;
        
        double loglikelihood = 0;
        double binLeft, binRight;
        
        for (size_t k = 0; k < K; ++k) { // Compute the likelihood

            size_t numPointsInBin = result.size[k+base];
            
            indexRight = indexLeft + numPointsInBin - 1;
            
            /* Use mid point inbetween two clusters as boundary
            binLeft = ( indexLeft == base ) ? 
                x[base] : (x[indexLeft-1] + x[indexLeft]) / 2;
            
            binRight = ( indexRight < N-1+base ) ? 
                (x[indexRight] + x[indexRight+1]) / 2 : x[N-1+base];
            */
            
            if(x[indexLeft] < x[indexRight]) {
                binLeft = x[indexLeft];
                binRight = x[indexRight];
            } else if(x[indexLeft] == x[indexRight]) {
                binLeft = ( indexLeft == base ) ?
                    x[base] : (x[indexLeft-1] + x[indexLeft]) / 2;
                binRight = ( indexRight < N-1+base ) ?
                    (x[indexRight] + x[indexRight+1]) / 2 : x[N-1+base];
            } else {
                throw "ERROR: binLeft > binRight";
                // cout << "ERROR: binLeft > binRight" << endl;
            }
            
            double binWidth = binRight - binLeft;
            
            if(method == "uniform") {
                loglikelihood += numPointsInBin * std::log(numPointsInBin / binWidth / N);
            } else if(method == "normal") {

                double mean = 0.0;
                double variance = 0.0;

                for (size_t i = indexLeft; i <= indexRight; ++i) {
                    mean += x[i];
                    variance += x[i] * x[i];
                }
                mean /= numPointsInBin;

                if (numPointsInBin > 1) {
                    variance = (variance - numPointsInBin * mean * mean)
                        / (numPointsInBin - 1);
                } else {
                    variance = 0;
                }
                
                if (variance > 0) {
                    for (size_t i = indexLeft; i <= indexRight; ++i) {
                        loglikelihood += - (x[i] - mean) * (x[i] - mean)
                            / (2.0 * variance);
                    }
                    loglikelihood += numPointsInBin
                        * (std::log(numPointsInBin / (double) N)
                           - 0.5 * std::log ( 2 * M_PI * variance));
                } else {
                    loglikelihood += numPointsInBin * std::log(1.0 / binWidth / N);
                }
                
            } else {
            	throw "ERROR: Wrong likelihood method!";
                // cout << "ERROR: Wrong likelihood method" << endl;
            }
            
            indexLeft = indexRight + 1;
        }
        
        double BIC = 0.0;
        
        // Compute the Bayesian information criterion
        if (method == "uniform") {
            BIC = 2 * loglikelihood - (3 * K - 1) * std::log((double)N);  // K-1
        } else if(method == "normal") {
            BIC = 2 * loglikelihood - (3 * K - 1) * std::log((double)N);  //(K*3-1)
        }
        
        // cout << ", Loglh=" << loglikelihood << ", BIC=" << BIC << endl;
        
        if (K == Kmin) {
            maxBIC = BIC;
            Kopt = Kmin;
        } else {
            if (BIC > maxBIC) {
                maxBIC = BIC;
                Kopt = K;
            }
        }
    }
    
    return Kopt;
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

ClusterResult
kmeans_1d_dp(const std::vector<double> & x, size_t Kmin, size_t Kmax)
{
    // Input:
    //  x -- a vector of numbers, not necessarily sorted
    //  Kmin -- the minimum number of clusters expected
    //  Kmax -- the maximum number of clusters expected
    // NOTE: All vectors in this program is considered starting at position 1,
    //       position 0 is not used.
 
    ClusterResult result;
    const size_t N = x.size() - 1;  // N: is the size of input vector
    
    std::vector<double> x_sorted(x);
    std::sort(x_sorted.begin()+1, x_sorted.end());
    const size_t nUnique = numberOfUnique(x_sorted.begin()+1, x_sorted.end());
    
    Kmax = nUnique < Kmax ? nUnique : Kmax;
    
    if(nUnique > 1) { // The case when not all elements are equal.
        
        std::vector< std::vector< double > > D( (Kmax + 1), std::vector<double>(N + 1) );
        std::vector< std::vector< size_t > > B( (Kmax + 1), std::vector<size_t>(N + 1) );
        
        // Fill in dynamic programming matrix
        fill_dp_matrix(x_sorted, D, B);
        
        // Choose an optimal number of levels between Kmin and Kmax
        size_t Kopt = select_levels(x_sorted, B, Kmin, Kmax);

        if (Kopt < Kmax) { // Reform the dynamic programming matrix D and B
            B.erase(B.begin()+ Kopt + 1, B.end());
        }
        
        // Backtrack to find the clusters beginning and ending indices
        backtrack(x_sorted, B, result);
        
        // Perform clustering on the original data
        for(size_t i = 1; i < x.size(); ++i) {
            size_t indexLeft = 1;
            size_t indexRight;
            
            for (size_t k = 1; k < result.size.size(); ++k) {
                indexRight = indexLeft + result.size[k] - 1;
                if ( x[i] <= x_sorted[indexRight] ) {
                    result.cluster[i] = k;
                    break;
                }
                indexLeft = indexRight + 1;
            }
        }
        
    } else {  // A single cluster that contains all elements
        
        result.nClusters = 1;
        
        result.cluster = std::vector<size_t>(N + 1, 1);
        
        result.centers.resize(2);
        result.withinss.resize(2);
        result.size.resize(2);
        
        result.centers[1] = x[1];
        result.withinss[1] = 0.0;
        result.size[1] = N;
    }
    return result;
}  //end of kmeans_1d_dp()
