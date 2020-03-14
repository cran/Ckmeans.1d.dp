/* EWL2_dynamic_prog.cpp  --- Equally weighted L2 univariate k-means algorithms
 *
 * Copyright (C) 2016-2020 Mingzhou Song
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
//
// Joe Song
// Created: January 3, 2020

#include <string>
#include <iostream>
#include <cassert>
#include <cmath>

#include "EWL2.h"

namespace EWL2 {

void fill_dp_matrix(const std::vector<double> & x, // data
                    const std::vector<double> & w, // weight
                    std::vector< std::vector< ldouble > > & S,
                    std::vector< std::vector< size_t > > & J,
                    const std::string & method)
                    /*
                     x: One dimension vector to be clustered, must be sorted (in any order).
                     S: K x N matrix. S[q][i] is the sum of squares of the distance from
                     each x[i] to its cluster mean when there are exactly x[i] is the
                     last point in cluster q
                     J: K x N backtrack matrix

                     NOTE: All vector indices in this program start at position 0
                     */
{
        const int K = (int) S.size();
        const int N = (int) S[0].size();

        std::vector<ldouble> sum_x(N), sum_x_sq(N);

        std::vector<int> jseq;

        ldouble shift = x[N/2]; // median. used to shift the values of x to
        //  improve numerical stability

        sum_x[0] = x[0] - shift;
        sum_x_sq[0] = (x[0] - shift) * (x[0] - shift);

        S[0][0] = 0;
        J[0][0] = 0;

        for(int i = 1; i < N; ++i) {

                sum_x[i] = sum_x[i-1] + x[i] - shift;
                sum_x_sq[i] = sum_x_sq[i-1] + (x[i] - shift) * (x[i] - shift);

                // Initialize for q = 0
                S[0][i] = dissimilarity(0, i, sum_x, sum_x_sq); // ssq(0, i, sum_x, sum_x_sq, sum_w);
                J[0][i] = 0;
        }

#ifdef DEBUG
        for(size_t i=0; i<x.size(); ++i) {
                std::cout << x[i] << ',';
        }
        std::cout << std::endl;
        std::vector<std::vector<ldouble>> SS(S);
        std::vector<std::vector<size_t>> JJ(J);

#endif

        for(int q = 1; q < K; ++q) {
                int imin;
                if(q < K - 1) {
                        imin = std::max(1, q);
                } else {
                        // No need to compute S[K-1][0] ... S[K-1][N-2]
                        imin = N-1;
                }

#ifdef DEBUG
                // std::cout << std::endl << "q=" << q << ":" << std::endl;
#endif
                // fill_row_k_linear_recursive(imin, N-1, 1, q, jseq, S, J, sum_x, sum_x_sq);
                // fill_row_k_linear(imin, N-1, q, S, J, sum_x, sum_x_sq);
                if(method == "linear") {
                        fill_row_q_SMAWK(imin, N-1, q, S, J, sum_x, sum_x_sq);
                } else if(method == "loglinear") {
                        fill_row_q_log_linear(imin, N-1, q, q, N-1, S, J, sum_x, sum_x_sq);
                } else if(method == "quadratic") {
                        fill_row_q(imin, N-1, q, S, J, sum_x, sum_x_sq);
                } else {
                        throw std::string("ERROR: unknown method") + method + "!";
                }

#ifdef DEBUG

                fill_row_q_log_linear(imin, N-1, q, q, N-1, SS, JJ, sum_x, sum_x_sq);

                for(int i=imin; i<N; ++i) {
                        if(S[q][i] != SS[q][i] || J[q][i] != JJ[q][i]) {
                                std::cout << "ERROR: q=" << q << ", i=" << i << std::endl;
                                std::cout << "\tS=" << S[q][i] << "\tJ=" << J[q][i] << std::endl;
                                std::cout << "Truth\tSS=" << SS[q][i] << "\tJJ=" << JJ[q][i];
                                std::cout << std::endl;
                                assert(false);

                        } else {
                                /*
                                 std::cout << "OK: q=" << q << ", i=" << i << std::endl;
                                 std::cout << "\tS=" << S[q][i] << "\tJ=" << J[q][i] << std::endl;
                                 std::cout << "Truth\tSS=" << SS[q][i] << "\tJJ=" << JJ[q][i];
                                 std::cout << std::endl;
                                 */
                        }

                }
#endif
        }

#ifdef DEBUG
        std::cout << "Linear & log-linear code returned identical dp index matrix."
                  << std::endl;
#endif

}
}
