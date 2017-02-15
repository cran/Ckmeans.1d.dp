/* Ckmeans_1d_dp.h --- Head file for Ckmeans.1d.dp
 *  Declare wrap function "kmeans_1d_dp()"
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
 Haizhou Wang
 Computer Science Department
 New Mexico State University
 hwang@cs.nmsu.edu

 Created: Oct 10, 2010
 */

#include <cstddef> // For size_t
#include <vector>
#include <string>

// typedef long double ldouble;
typedef double ldouble;

inline ldouble ssq(
    const size_t j, const size_t i,
    const std::vector<ldouble> & sum_x, // running sum of xi
    const std::vector<ldouble> & sum_x_sq, // running sum of xi^2
    const std::vector<ldouble> & sum_w  // running sum of weights
)
{
  ldouble sji(0.0);

  if(sum_w.empty()) { // equally weighted version
    if(j >= i) {
      sji = 0.0;
    } else if(j > 0) {
      ldouble muji = (sum_x[i] - sum_x[j-1]) / (i - j + 1);
      sji = sum_x_sq[i] - sum_x_sq[j-1] - (i - j + 1) * muji * muji;
    } else {
      sji = sum_x_sq[i] - sum_x[i] * sum_x[i] / (i+1);
    }
  } else { // unequally weighted version
    if(sum_w[j] >= sum_w[i]) {
      sji = 0.0;
    } else if(j > 0) {
      ldouble muji = (sum_x[i] - sum_x[j-1]) / (sum_w[i] - sum_w[j-1]);
      sji = sum_x_sq[i] - sum_x_sq[j-1] - (sum_w[i] - sum_w[j-1]) * muji * muji;
    } else {
      sji = sum_x_sq[i] - sum_x[i] * sum_x[i] / sum_w[i];
    }
  }

  sji = (sji < 0) ? 0 : sji;
  return sji;
}

void fill_dp_matrix(const std::vector<double> & x,
                    const std::vector<double> & w,
                    std::vector< std::vector< ldouble > > & S,
                    std::vector< std::vector< size_t > > & J,
                    const std::string & method);

void backtrack(const std::vector<double> & x,
               const std::vector< std::vector< size_t > > & J,
               int* cluster, double* centers, double* withinss, int* count);

void backtrack(const std::vector<double> & x,
               const std::vector< std::vector< size_t > > & J,
               std::vector<size_t> & count);

void fill_row_q_SMAWK(int imin, int imax, int q,
                      std::vector< std::vector<ldouble> > & S,
                      std::vector< std::vector<size_t> > & J,
                      const std::vector<ldouble> & sum_x,
                      const std::vector<ldouble> & sum_x_sq,
                      const std::vector<ldouble> & sum_w);

void fill_row_q(int imin, int imax, int q,
                std::vector< std::vector<ldouble> > & S,
                std::vector< std::vector<size_t> > & J,
                const std::vector<ldouble> & sum_x,
                const std::vector<ldouble> & sum_x_sq,
                const std::vector<ldouble> & sum_w);


void fill_row_q_log_linear(
    int imin, int imax, int q, int jmin, int jmax,
    std::vector< std::vector<ldouble> > & S,
    std::vector< std::vector<size_t> > & J,
    const std::vector<ldouble> & sum_x,
    const std::vector<ldouble> & sum_x_sq,
    const std::vector<ldouble> & sum_w);

/* One-dimensional cluster algorithm implemented in C++ */
/* x is input one-dimensional vector and
 Kmin and Kmax stand for the range for the number of clusters*/
void kmeans_1d_dp(
    const double *x, const size_t N,
    const double * y,
    size_t Kmin, size_t Kmax,
    int* cluster, double* centers,
    double* withinss, int* size, double* BIC,
    const std::string & estimate_k,
    const std::string & method);


void backtrack(
    const std::vector<double> & x,
    const std::vector< std::vector< size_t > > & J,
    std::vector<size_t> & counts, const int K);

size_t select_levels(
    const std::vector<double> & x,
    const std::vector< std::vector< size_t > > & J,
    size_t Kmin, size_t Kmax, double *BIC);

size_t select_levels_3_4_12(
    const std::vector<double> & x,
    const std::vector< std::vector< size_t > > & J,
    size_t Kmin, size_t Kmax, double *BIC);

void fill_weighted_dp_matrix(
    const std::vector<double> & x,
    const std::vector<double> & y,
    std::vector< std::vector< ldouble > > & S,
    std::vector< std::vector< size_t > > & J);

void backtrack_weighted(
    const std::vector<double> & x, const std::vector<double> & y,
    const std::vector< std::vector< size_t > > & J,
    std::vector<size_t> & counts, std::vector<double> & weights,
    const int K);

void backtrack_weighted(
    const std::vector<double> & x, const std::vector<double> & y,
    const std::vector< std::vector< size_t > > & J,
    int* cluster, double* centers, double* withinss, int* weights);

size_t select_levels_weighted(
    const std::vector<double> & x, const std::vector<double> & y,
    const std::vector< std::vector< size_t > > & J,
    size_t Kmin, size_t Kmax, double *BIC);

size_t select_levels_weighted_3_4_12(
    const std::vector<double> & x, const std::vector<double> & y,
    const std::vector< std::vector< size_t > > & J,
    size_t Kmin, size_t Kmax, double *BIC);

void range_of_variance(
    const std::vector<double> & x,
    double & variance_min, double & variance_max);

