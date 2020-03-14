/* EWL2.h --- Head file for equally weighted L2 univariate kmeans clustering
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
 Mingzhou Song
 Computer Science Department
 New Mexico State University
 joemsong@cs.nmsu.edu

 Created: January 3, 2020
 */

#include <cstddef> // For size_t
#include <vector>
#include <string>

#include "EWL2_within_cluster.h"

namespace EWL2 {

void fill_dp_matrix(
    const std::vector<double> & x,
    const std::vector<double> & w,
    std::vector< std::vector< ldouble > > & S,
    std::vector< std::vector< size_t > > & J,
    const std::string & method);

void fill_row_q_SMAWK(
    int imin, int imax, int q,
    std::vector< std::vector<ldouble> > & S,
    std::vector< std::vector<size_t> > & J,
    const std::vector<ldouble> & sum_x,
    const std::vector<ldouble> & sum_x_sq);

void fill_row_q(
    int imin, int imax, int q,
    std::vector< std::vector<ldouble> > & S,
    std::vector< std::vector<size_t> > & J,
    const std::vector<ldouble> & sum_x,
    const std::vector<ldouble> & sum_x_sq);

void fill_row_q_log_linear(
    int imin, int imax, int q, int jmin, int jmax,
    std::vector< std::vector<ldouble> > & S,
    std::vector< std::vector<size_t> > & J,
    const std::vector<ldouble> & sum_x,
    const std::vector<ldouble> & sum_x_sq);

}
