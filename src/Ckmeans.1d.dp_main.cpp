/* Ckmeans_1d_dp_main.cpp --- wrapper function for "kmeans_1d_dp()"
 *
 * Copyright (C) 2010-2016 Mingzhou Song and Haizhou Wang
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

/*
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
 May 29, 2017. Joe Song. Change size from integer to double.
 September 21, 2018. Joe Song. Changed from C to Rcpp interface.
 */

#include <Rcpp.h>
using namespace Rcpp;

#include <string>
#include <iostream>

#include "Ckmeans.1d.dp.h"

/*Wrapper function to call kmeans_1d_dp()*/
/*
 x: An array containing input data to be clustered.
 length: length of the one dimensional array.
 minK: Minimum number of clusters.
 maxK: Maximum number of clusters.
 cluster:  An array of cluster IDs for each point in x.
 centers:  An array of centers for each cluster.
 withinss: An array of within-cluster sum of squares for each cluster.
 size:     An array of (weighted) sizes of each cluster.
 */

// [[Rcpp::export]]
List Ckmeans_1d_dp(const NumericVector & x, const size_t length,
                   const NumericVector & y, const size_t ylength,
                   const size_t minK, const size_t maxK,
                   IntegerVector & cluster,
                   NumericVector & centers,
                   NumericVector & withinss,
                   NumericVector & size, // int* size,
                   NumericVector & BICs,
                   const std::string & estimate_k, const std::string & method)
{
  // std::cout << method[0] << std::endl;

  const double * xp(x.begin()), * yp (y.begin());
  double * center_p (centers.begin()), * sp (size.begin());
  double * bp (BICs.begin()), * wp (withinss.begin());
  int * cluster_p (cluster.begin());

  // Call C++ version one-dimensional clustering algorithm*/
  if(ylength != length) { yp = 0; }

  kmeans_1d_dp(xp, length, yp, minK, maxK,
               cluster_p, center_p, wp, sp, bp,
               estimate_k, method, L2);

  // Change the cluster numbering from 0-based to 1-based
  for(size_t i=0; i < length; ++i) {
    cluster[i] ++;
  }

  List result;
  result["centers"] = centers;
  result["cluster"] = cluster;
  result["BIC"] = BICs;
  result["withinss"] = withinss;
  result["size"] = size;

  return(result);
}


// [[Rcpp::export]]
List Ckmedian_1d_dp(const NumericVector & x, const size_t length,
                   const NumericVector & y, const size_t ylength,
                   const size_t minK, const size_t maxK,
                   IntegerVector & cluster,
                   NumericVector & centers,
                   NumericVector & withinss,
                   NumericVector & size, // int* size,
                   NumericVector & BICs,
                   const std::string & estimate_k, const std::string & method)
{
  // std::cout << method[0] << std::endl;

  const double * xp(x.begin()), * yp (y.begin());
  double * center_p (centers.begin()), * sp (size.begin());
  double * bp (BICs.begin()), * wp (withinss.begin());
  int * cluster_p (cluster.begin());

  // Call C++ version one-dimensional clustering algorithm*/
  if(ylength != length) { yp = 0; }

  kmeans_1d_dp(xp, length, yp, minK, maxK,
               cluster_p, center_p, wp, sp, bp,
               estimate_k, method, L1);

  // Change the cluster numbering from 0-based to 1-based
  for(size_t i=0; i < length; ++i) {
    cluster[i] ++;
  }
  List result;
  result["centers"] = centers;
  result["cluster"] = cluster;
  result["BIC"] = BICs;
  result["withinss"] = withinss;
  result["size"] = size;

  return(result);
}

// [[Rcpp::export]]
List Cksegs_1d_dp(const NumericVector & x, const size_t length,
                   const NumericVector & y, const size_t ylength,
                   const size_t minK, const size_t maxK,
                   IntegerVector & cluster,
                   NumericVector & centers,
                   NumericVector & withinss,
                   NumericVector & size, // int* size,
                   NumericVector & BICs,
                   const std::string & estimate_k, const std::string & method)
{
  // std::cout << method[0] << std::endl;

  const double * xp(x.begin()), * yp (y.begin());
  double * center_p (centers.begin()), * sp (size.begin());
  double * bp (BICs.begin()), * wp (withinss.begin());
  int * cluster_p (cluster.begin());

  // Call C++ version one-dimensional clustering algorithm*/
  if(ylength != length) { yp = 0; }

  kmeans_1d_dp(xp, length, yp, minK, maxK,
               cluster_p, center_p, wp, sp, bp,
               estimate_k, method, L2Y);

  // Change the cluster numbering from 0-based to 1-based
  for(size_t i=0; i < length; ++i) {
    cluster[i] ++;
  }

  List result;
  result["centers"] = centers;
  result["cluster"] = cluster;
  result["BIC"] = BICs;
  result["withinss"] = withinss;
  result["size"] = size;

  return(result);

}
