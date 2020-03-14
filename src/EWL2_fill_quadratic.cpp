/* EWL2_fill_quadratic.cpp --- an iterative algorithm to compute a
 *   row in the dynamic programming matrix in O(n^2) time.
 *
 * Copyright (C) 2016--2020 Mingzhou Song
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
// Revised from fill_quadratic.cpp

#include <vector>
#include <algorithm>
#include <iostream>
#include <cassert>
#include <stack>
#include <cmath>

#include "EWL2.h"

namespace EWL2 {
void fill_row_q(int imin, int imax, int q,
                std::vector< std::vector<ldouble> > & S,
                std::vector< std::vector<size_t> > & J,
                const std::vector<ldouble> & sum_x,
                const std::vector<ldouble> & sum_x_sq)
{
  // Assumption: each cluster must have at least one point.
  for(int i=imin; i<=imax; ++i) {
    S[q][i] = S[q-1][i-1];
    J[q][i] = i;
    int jmin = std::max(q, (int)J[q-1][i]);
    for(int j=i-1; j>=jmin; --j) {
      ldouble Sj(S[q-1][j-1] +
        dissimilarity(j, i, sum_x, sum_x_sq)
                   // ssq(j, i, sum_x, sum_x_sq, sum_w)
      );

      if(Sj < S[q][i]) {
        S[q][i] = Sj;
        J[q][i] = j;
      }
    }
  }
}
}
