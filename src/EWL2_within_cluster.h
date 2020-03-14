/* EWL2_within_cluster.h---dissimilarity within a cluster
 *
 * Copyright (C) 2017--2020 Mingzhou Song
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
// Created: January 3, 2020. simplified from within_cluster.h

#include <vector>
#include <string>

#include "precision.h"

namespace EWL2 {

inline ldouble ssq(
  const size_t j, const size_t i,
  const std::vector<ldouble> & sum_x, // running sum of xi
  const std::vector<ldouble> & sum_x_sq // running sum of xi^2
)
{
  ldouble sji(0.0);

  if(j >= i) {
    sji = 0.0;
  } else if(j > 0) {
    ldouble muji = (sum_x[i] - sum_x[j-1]) / (i - j + 1);
    sji = sum_x_sq[i] - sum_x_sq[j-1] - (i - j + 1) * muji * muji;
  } else {
    sji = sum_x_sq[i] - sum_x[i] * sum_x[i] / (i+1);
  }

  sji = (sji < 0) ? 0 : sji;
  return sji;
}

inline ldouble dissimilarity(
  const size_t j, const size_t i,
  const std::vector<ldouble> & sum_x, // running sum of xi
  const std::vector<ldouble> & sum_x_sq // running sum of xi^2
)
{
  ldouble d=0;

  d = ssq(j, i, sum_x, sum_x_sq);

  return d;
}
}
