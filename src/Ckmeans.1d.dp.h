/*
 Ckmeans_1d_dp.h --- Head file for Ckmeans.1d.dp
 Declare wrap function "kmeans_1d_dp()"

 Haizhou Wang
 Computer Science Department
 New Mexico State University
 hwang@cs.nmsu.edu

 Created: Oct 10, 2010
 */

#include <cstddef> // For size_t
#include <vector>

/* One-dimensional cluster algorithm implemented in C++ */
/* x is input one-dimensional vector and
 Kmin and Kmax stand for the range for the number of clusters*/
void kmeans_1d_dp(
    const double *x, const size_t N, const double * y,
    size_t Kmin, size_t Kmax,
    int* cluster, double* centers, double* withinss, int* size);


void backtrack(
    const std::vector<double> & x,
    const std::vector< std::vector< size_t > > & J,
    std::vector<size_t> & counts);

size_t select_levels(
    const std::vector<double> & x,
    const std::vector< std::vector< size_t > > & J,
    size_t Kmin, size_t Kmax);

void fill_weighted_dp_matrix(
    const std::vector<double> & x,
    const std::vector<double> & y,
    std::vector< std::vector< double > > & S,
    std::vector< std::vector< size_t > > & J);

void backtrack_weighted(
    const std::vector<double> & x, const std::vector<double> & y,
    const std::vector< std::vector< size_t > > & J,
    std::vector<size_t> & counts, std::vector<double> & weights);

void backtrack_weighted(
    const std::vector<double> & x, const std::vector<double> & y,
    const std::vector< std::vector< size_t > > & J,
    int* cluster, double* centers, double* withinss, int* weights);

size_t select_levels_weighted(
    const std::vector<double> & x, const std::vector<double> & y,
    const std::vector< std::vector< size_t > > & J,
    size_t Kmin, size_t Kmax);
