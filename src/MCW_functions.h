//  Ckmeans.1d.dp
//  MCW_functions.h
//
//  Created by Hua Zhong on 1/21/18.
//  Copyright © 2018 Hua Zhong. All rights reserved.
//


#include <cstddef> // For size_t
#include <vector>
#include <string>
#include <iostream>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <cstring>
#include <stdlib.h>
#include <limits>
#include <Rcpp.h>

// typedef long double ldouble;
typedef double ldouble;
static const ldouble myinf = std::numeric_limits<ldouble>::max();

ldouble MCW_dissimilarity(const size_t j, const size_t i,
                          const std::vector< std::vector<ldouble> > & sum_x, // running sum of xi
                          const std::vector< ldouble > & sum_x_sq, // running sum of xi^2
                          const std::vector< std::vector<ldouble> > & sum_w  // running sum of weights
);

// Linear algorithm; Link to fill_SMAWK.cpp
void MCW_fill_row_q_SMAWK(int imin, int imax, int q,
                          std::vector< std::vector<ldouble> > & S,
                          std::vector< std::vector<size_t> > & J,
                          const std::vector< std::vector<ldouble> > & sum_x,
                          const std::vector< ldouble > & sum_x_sq,
                          const std::vector< std::vector<ldouble> > & sum_w);

// Quadratic algorithm
void MCW_fill_row_q(int imin, int imax, int q,
                    std::vector< std::vector<ldouble> > & S,
                    std::vector< std::vector<size_t> > & J,
                    const std::vector< std::vector<ldouble> > & sum_x,
                    const std::vector< ldouble > & sum_x_sq,
                    const std::vector< std::vector<ldouble> > & sum_w);

void MCW_fill_dp_matrix(const std::vector<double> & x, // data
                        const std::vector< std::vector<double> > & w, // weight
                        std::vector< std::vector< ldouble > > & S,
                        std::vector< std::vector< size_t > > & J,
                        const std::string & method);

// Used by select_levels() to select a best k via BIC
void MCW_backtrack_BIC_weighted(const std::vector<double> & x, const std::vector<double> & y,
                                const std::vector< std::vector< size_t > > & J,
                                std::vector<size_t> & counts, std::vector<double> & weights,
                                const int K);

void MCW_backtrack_weighted(const std::vector<double> & x,
                            const std::vector< std::vector<double> > & y,
                            const std::vector< std::vector< size_t > > & J,
                            std::vector<int> & cluster,
                            std::vector<double> & centers,
                            std::vector<double> & withinss,
                            std::vector<double> & weights);

// Used by select_levels() to select a best k via BIC
void MCW_shifted_data_variance_weighted(const std::vector<double> & x,
                                        const std::vector<double> & y,
                                        const double total_weight,
                                        const size_t left,
                                        const size_t right,
                                        double & mean, double & variance);

// Hua added, Apr 5, 2018
// Take multiple dimensional of weights
size_t MCW_select_levels_BIC(const std::vector<double> & x,
                             const std::vector< std::vector<double> > & y,
                             const std::vector< std::vector< size_t > > & J,
                             const size_t Kmin, const size_t Kmax,
                             std::vector<double> & BIC);

void MCW_optimal_zoning(const std::vector<double> & x,
                        const std::vector<std::vector<double> > & y,
                        size_t Kmin,
                        size_t Kmax,
                        std::vector<int> & cluster,
                        std::vector<double> & centers,
                        std::vector<double> & withinss,
                        std::vector<double> & size,
                        std::vector<double> & BIC,
                        const std::string estimate_k,
                        const std::string method);

void MCW_optimal_zoning_main(const std::vector<double> & x,
                             const std::vector<std::vector<double> > & y,
                             size_t Kmin,
                             size_t Kmax,
                             std::vector<int> & cluster,
                             std::vector<double> & centers,
                             std::vector<double> & withinss,
                             std::vector<double> & size,
                             std::vector<double> & BIC,
                             const std::string estimate_k,
                             const std::string method);

