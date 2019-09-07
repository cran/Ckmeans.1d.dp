//  Ckmeans.1d.dp
//  MCW_functions.cpp
//
//  Created by Hua Zhong on 1/21/18.
//  Copyright © 2018 Hua Zhong. All rights reserved.
//


#include "MCW_functions.h"

// Hua added, Apr 5, 2018
// Take multiple dimensional of weights
ldouble MCW_dissimilarity(const size_t j, const size_t i,
                          const std::vector< std::vector<ldouble> > & sum_x, // running sum of xi
                          const std::vector< ldouble > & sum_x_sq, // running sum of xi^2
                          const std::vector< std::vector<ldouble> > & sum_w   // running sum of weights
)
{
    ldouble sji(0.0);

    if(j >= i){
        sji = 0.0;
    }else if (j == 0) {
        sji = sum_x_sq[i]; // if all sun_w[][i] == 0 --> sum_x_sq[i] == 0
        for (size_t l = 0; l < sum_x.size(); l++) {
            if(sum_w[l][i] > 0) sji -= sum_x[l][i] * sum_x[l][i] / sum_w[l][i];
        }
    }else{
        sji = sum_x_sq[i] - sum_x_sq[j-1];
        for (size_t l = 0; l < sum_x.size(); l++) {
            if(sum_w[l][i] > sum_w[l][j-1]){
                sji -= (sum_x[l][i] - sum_x[l][j-1]) *
                        (sum_x[l][i] - sum_x[l][j-1]) /
                        (sum_w[l][i] - sum_w[l][j-1]);
            }
        }
    }

    sji = (sji < 0) ? 0 : sji;
    return sji;
}

// End of change.


// Hua added, Apr 5, 2018
// Take multiple dimensional of weights
void MCW_fill_row_q(int imin, int imax, int q,
                    std::vector< std::vector<ldouble> > & S,
                    std::vector< std::vector<size_t> > & J,
                    const std::vector< std::vector<ldouble> > & sum_x,
                    const std::vector< ldouble > & sum_x_sq,
                    const std::vector< std::vector<ldouble> > & sum_w)
{
    // Assumption: each cluster must have at least one point.
    for(int i=imin; i<=imax; ++i) {
        S[q][i] = S[q-1][i-1];
        J[q][i] = i;
        int jmin = std::max(q, (int)J[q-1][i]);
        for(int j=i-1; j>=jmin; --j) {
            ldouble Sj(S[q-1][j-1] + MCW_dissimilarity(j, i, sum_x, sum_x_sq, sum_w));
            if(Sj < S[q][i]) {
                S[q][i] = Sj;
                J[q][i] = j;
            }
        }
    }
}
// End of change.


// Hua added, Apr 5, 2018
// Take multiple dimensional of weights
void MCW_fill_dp_matrix(const std::vector<double> & x, // data
                        const std::vector< std::vector<double> > & w, // weight
                        std::vector< std::vector< ldouble > > & S,
                        std::vector< std::vector< size_t > > & J,
                        const std::string & method)
/*
 x: One dimension vector to be clustered, must be sorted (in any order).
 w: Multiple dimensional of weights.
 S: K x N matrix. S[q][i] is the sum of squares of the distance from
 each x[i] to its cluster mean when there are exactly x[i] is the
 last point in cluster q
 J: K x N backtrack matrix

 NOTE: All vector indices in this program start at position 0
 */
{
    const int K = (int) S.size();
    const int N = (int) S[0].size();
    const int W = (int) w.size();

    std::vector< std::vector<ldouble> > sum_x(W, std::vector<ldouble>(N, 0));
    std::vector< ldouble > sum_x_sq(N, 0);
    std::vector< std::vector<ldouble> > sum_w(W, std::vector<ldouble>(N, 0));

    ldouble shift = x[N/2]; // median. used to shift the values of x to
    //  improve numerical stability

    for (size_t l = 0; l < W; l++) {
        sum_x[l][0] = w[l][0] * (x[0] - shift);
        sum_w[l][0] = w[l][0];
        sum_x_sq[0] += w[l][0] * (x[0] - shift) * (x[0] - shift);
    }

    S[0][0] = 0;
    J[0][0] = 0;

    for(int i = 1; i < N; ++i) {
        sum_x_sq[i] = sum_x_sq[i-1];

        for (size_t l = 0; l < W; l++) {
            sum_x[l][i] = sum_x[l][i-1] + w[l][i] * (x[i] - shift);
            sum_w[l][i] = sum_w[l][i-1] + w[l][i];
            sum_x_sq[i] += w[l][i] * (x[i] - shift) * (x[i] - shift);
        }

        // Initialize for q = 0
        S[0][i] = MCW_dissimilarity(0, i, sum_x, sum_x_sq, sum_w);
        J[0][i] = 0;
    }

    for(int q = 1; q < K; ++q) {
        int imin;
        if(q < K - 1) {
            imin = std::max(1, q);
        } else {
            // No need to compute S[K-1][0] ... S[K-1][N-2]
            imin = N-1;
        }

        if(method == "linear") {
            MCW_fill_row_q_SMAWK(imin, N-1, q, S, J, sum_x, sum_x_sq, sum_w);
        } else if(method == "loglinear") {
            //
        } else if(method == "quadratic") {
            MCW_fill_row_q(imin, N-1, q, S, J, sum_x, sum_x_sq, sum_w);
        } else {
            throw std::string("ERROR: unknown method") + method + "!";
        }
    }
}
// End of change.

// Hua added, Apr 5, 2018
// Take multiple dimensional of weights
void MCW_backtrack_weighted(const std::vector<double> & x,
                            const std::vector< std::vector<double> > & y,
                            const std::vector< std::vector< size_t > > & J,
                            std::vector<int> & cluster,
                            std::vector<double> & centers,
                            std::vector<double> & withinss,
                            std::vector<double> & weights)
{
    const int K = (int) J.size();
    const size_t N = J[0].size();
    size_t cluster_right = N-1;
    size_t cluster_left;

    // Backtrack the clusters from the dynamic programming matrix
    for(int k = K-1; k >= 0; --k) {
        weights[k] = 0.0;
        withinss[k] = 0.0;

        cluster_left = J[k][cluster_right];

        for(size_t i = cluster_left; i <= cluster_right; ++i) cluster[i] = k;

        std::vector<double> center_1d(y.size());
        double center_1d_sum = 0.0;

        for(size_t l = 0; l < y.size(); l++){
            double sum = 0.0;
            double weight = 0.0;

            for(size_t i = cluster_left; i <= cluster_right; ++i) {
                sum += x[i] * y[l][i];
                weight += y[l][i];
            }

            center_1d[l] = sum / weight;
            center_1d_sum += center_1d[l];
            weights[k] += weight;

            for(size_t i = cluster_left; i <= cluster_right; ++i){
                withinss[k] += y[l][i] * (x[i] - center_1d[l]) * (x[i] - center_1d[l]);
            }
        }

        centers[k] = center_1d_sum / y.size();

        if(k > 0) {
            cluster_right = cluster_left - 1;
        }
    }
}

// End of change.

void MCW_backtrack_BIC_weighted(const std::vector<double> & x, const std::vector<double> & y,
                                const std::vector< std::vector< size_t > > & J,
                                std::vector<size_t> & counts, std::vector<double> & weights,
                                const int K)
{
    // const int K = (int) J.size();
    const size_t N = J[0].size();
    size_t cluster_right = N-1;
    size_t cluster_left;

    // Backtrack the clusters from the dynamic programming matrix
    for(int k = K-1; k >= 0; --k) {
        cluster_left = J[k][cluster_right];
        counts[k] = cluster_right - cluster_left + 1;

        weights[k] = 0;
        for(size_t i = cluster_left; i <= cluster_right; ++i) {
            weights[k] += y[i];
        }

        if(k > 0) {
            cluster_right = cluster_left - 1;
        }
    }
}

void MCW_shifted_data_variance_weighted(const std::vector<double> & x,
                                        const std::vector<double> & y,
                                        const double total_weight,
                                        const size_t left,
                                        const size_t right,
                                        double & mean, double & variance)
{
    double sum = 0.0;
    double sumsq = 0.0;

    mean = 0.0;
    variance = 0.0;

    size_t n = right - left + 1;

    if(right >= left) {

        double median = x[(left + right) / 2];

        for (size_t i = left; i <= right; ++i) {
            sum += (x[i] - median) * y[i];
            sumsq += (x[i] - median) * (x[i] - median) * y[i];
        }
        mean = sum / total_weight + median;

        if (n > 1) {
            variance = (sumsq - sum * sum / total_weight) / (total_weight - 1);
        }
    }
}

// Hua added, Apr 5, 2018
// Take multiple dimensional of weights
size_t MCW_select_levels_BIC(const std::vector<double> & x,
                             const std::vector< std::vector<double> > & y,
                             const std::vector< std::vector< size_t > > & J,
                             const size_t Kmin, const size_t Kmax,
                             std::vector<double> & BIC)
{
    const size_t N = x.size();

    if (Kmin > Kmax || N < 2) {
        return std::min(Kmin, Kmax);
    }

    size_t Kopt = Kmin;

    long double maxBIC(0.0);


    for(size_t K = Kmin; K <= Kmax; ++K) {
        long double bic = 0.0;

        for (size_t l = 0; l < y.size(); l++) {
            std::vector<double> lambda(Kmax);
            std::vector<double> mu(Kmax);
            std::vector<double> sigma2(Kmax);
            std::vector<double> coeff(Kmax);
            std::vector<size_t> counts(Kmax);
            std::vector<double> weights(Kmax);

            std::vector<double> y_1d = y[l];

            MCW_backtrack_BIC_weighted(x, y_1d, J, counts, weights, (int)K);

            double totalweight = 0.0;

            for(size_t k=0; k<K; k++) {
                totalweight += weights[k];
            }


            size_t indexLeft = 0;
            size_t indexRight;

            for (size_t k = 0; k < K; ++k) { // Estimate GMM parameters first

                lambda[k] = weights[k] / totalweight;

                indexRight = indexLeft + counts[k] - 1;

                MCW_shifted_data_variance_weighted(x, y_1d, weights[k], indexLeft, indexRight, mu[k], sigma2[k]);

                if(sigma2[k] == 0 || counts[k] == 1) {

                    double dmin;

                    if(indexLeft > 0 && indexRight < N-1) {
                        dmin = std::min(x[indexLeft] - x[indexLeft-1], x[indexRight+1] - x[indexRight]);
                    } else if(indexLeft > 0) {
                        dmin = x[indexLeft] - x[indexLeft-1];
                    } else {
                        dmin = x[indexRight+1] - x[indexRight];
                    }

                    // std::cout << "sigma2[k]=" << sigma2[k] << "==>";
                    if(sigma2[k] == 0) sigma2[k] = dmin * dmin / 4.0 / 9.0 ;
                    if(counts[k] == 1) sigma2[k] = dmin * dmin;
                    // std::cout << sigma2[k] << std::endl;
                }

                /*
                 if(sigma2[k] == 0) sigma2[k] = variance_min;
                 if(size[k] == 1) sigma2[k] = variance_max;
                 */

                coeff[k] = lambda[k] / std::sqrt(2.0 * M_PI * sigma2[k]);

                indexLeft = indexRight + 1;
            }

            long double loglikelihood = 0;

            for (size_t i=0; i<N; ++i) {
                long double L=0;
                for (size_t k = 0; k < K; ++k) {
                    L += coeff[k] * std::exp(- (x[i] - mu[k]) * (x[i] - mu[k]) / (2.0 * sigma2[k]));
                }
                loglikelihood += y_1d[i] * std::log(L);

                //Rcpp::Rcout << "L:" << L << std::endl;
            }

            // Compute the Bayesian information criterion
            bic += 2 * loglikelihood - (3 * K - 1) * std::log(totalweight);  //(K*3-1)

            // std::cout << "k=" << K << ": Loglh=" << loglikelihood << ", BIC=" << BIC << std::endl;
            //Rcpp::Rcout << "k=" << K << ": Loglh=" << loglikelihood << ", BIC=" << bic << std::endl;
        }

        if (K == Kmin) {
            maxBIC = bic;
            Kopt = Kmin;
        } else {
            if (bic > maxBIC) {
                maxBIC = bic;
                Kopt = K;
            }
        }

        BIC[K-Kmin] = (double)bic;
    }
    return Kopt;
}
// End of change.


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

// Hua added, Apr 5, 2018
// Take multiple dimensional of weights
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
                        const std::string method)
{
    // Input:
    //  x -- a vector of double precision numbers, not necessarily sorted
    //  y -- an 2D vector, multiple dimensional of weights
    //  Kmin -- the minimum number of clusters expected
    //  Kmax -- the maximum number of clusters expected
    // NOTE: All vectors in this program is considered starting at position 0.

    size_t N = x.size();

    const size_t nUnique = numberOfUnique(x.begin(), x.end());

    Kmax = nUnique < Kmax ? nUnique : Kmax;

    if(nUnique > 1) { // The case when not all elements are equal.

        std::vector< std::vector< ldouble > > S( Kmax, std::vector<ldouble>(N) );
        std::vector< std::vector< size_t > > J( Kmax, std::vector<size_t>(N) );

        size_t Kopt;

        MCW_fill_dp_matrix(x, y, S, J, method);

        if(estimate_k=="BIC") {
            // Choose an optimal number of levels between Kmin and Kmax
            Kopt = MCW_select_levels_BIC(x, y, J, Kmin, Kmax, BIC);
        } else {
            // std::cerr << "ERROR: No such method estimating k!" << std::endl;
            // exit(0);

            Rcpp::stop("ERROR: No such method estimating k!");
        }


        if (Kopt < Kmax) { // Reform the dynamic programming matrix S and J
            J.erase(J.begin() + Kopt, J.end());
        }


        // Backtrack to find the clusters beginning and ending indices
        MCW_backtrack_weighted(x, y, J, cluster, centers, withinss, size);

#ifdef DEBUG
        std::cout << "backtrack done." << std::endl;
#endif

    } else {  // A single cluster that contains all elements

        for(size_t i=0; i<N; ++i) {
            cluster[i] = 0;
        }

        centers[0] = x[0];
        withinss[0] = 0.0;
        size[0] = N * y[0][0];
    }
}
// End of change.

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
                             const std::string method)
{
    // Call C++ version multi-dimensional weighted clustering algorithm*/
    // check to see if unequal weight is provided
    if(y.empty()) {
        // Must not be empty
        // std::cerr << "ERROR: Weight matrix must not be empty!" << std::endl;
        // exit(0);

        Rcpp::stop("ERROR: Weight matrix must not be empty!");
    }else if(x.size() != y[0].size()) {
        // std::cerr << "ERROR: Weight matrix y must have the same rowsize as the length of x!" << std::endl;
        // exit(0);

        Rcpp::stop("ERROR: Weight matrix y must have the same rowsize as the length of x!");
    }

    MCW_optimal_zoning(x, y, Kmin, Kmax, cluster, centers, withinss, size, BIC, estimate_k, method);

    // Change the cluster numbering from 0-based to 1-based
    for(size_t i=0; i< x.size(); ++i) {
        cluster[i] ++;
    }

    //Resize the vectors, remove 0 at the back.
    int cluster_num = cluster.back();
    centers.resize(cluster_num);
    centers.resize(cluster_num);
    size.resize(cluster_num);
}
