//  Ckmeans.1d.dp
//  MCW_main.cpp
//
//  Created by Hua Zhong on 1/21/18.
//  Copyright © 2018 Hua Zhong. All rights reserved.
//

//#include <numeric>
//#include <algorithm>
#include "MCW_functions.h"
#include <iostream>
#include <Rcpp.h>
using namespace Rcpp;
//using namespace std;

// [[Rcpp::export]]
List MCW_main(const NumericVector & x,
              const NumericMatrix & y,
              size_t Kmin,
              size_t Kmax,
              const std::string estimate_k,
              const std::string method
){

  std::vector< double > x_new(x.length(), 0.0);
  for(size_t i = 0; i < (size_t) x.length(); ++i) {
    x_new[i] = x(i);
  }

  std::vector< std::vector< double > > y_new(y.ncol(), std::vector<double>(y.nrow()));
  for(size_t i = 0; i < (size_t) y.ncol(); ++i) {
    for(size_t j = 0; j < (size_t) y.nrow(); ++j) {
      y_new[i][j] = y(j, i);
    }
  }

  // If negtive weights exist in a sample of weights, scale the sample to positive
  if(false){
    for (size_t i=0; i<y_new.size(); i++) {
      double y_one_min = *std::min_element(y_new[i].begin(), y_new[i].end());
      if(y_one_min < 0){
        for(size_t j = 0; j < y_new[i].size(); j++) {
          y_new[i][j] = y_new[i][j] - y_one_min;
        }
      }
    }
  }

  // Scale the sum of a sample weightsto the length of the sequence
  if(false){
    for (size_t i=0; i<y_new.size(); i++) {
      double factor = (double) x_new.size() / std::accumulate(y_new[i].begin(), y_new[i].end(), 0.0);
      for(size_t j = 0; j < y_new[i].size(); j++) {
        y_new[i][j] = y_new[i][j] * factor;
      }
    }
  }

  std::vector< int > cluster(x.length(), 0);

  std::vector< double > centers(Kmax, 0);

  std::vector< double > withinss(Kmax, 0);

  std::vector< double > size(Kmax, 0);

  std::vector< double > BIC(Kmax-Kmin+1, 0);

  MCW_optimal_zoning_main(x_new, y_new, Kmin, Kmax, cluster, centers,
                          withinss, size, BIC, estimate_k, method);

  return List::create(Named("centers") = centers,
                      Named("cluster") = cluster,
                      Named("BIC") = BIC,
                      Named("withinss") = withinss,
                      Named("size") = size);
}
