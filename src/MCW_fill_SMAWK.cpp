//  Ckmeans.1d.dp
//  MCW_fill_SMAWK.cpp
//
//  Created by Hua Zhong on 1/21/18.
//  Modified from fill_SMAWK.cpp
//
//  Copyright © 2018 Hua Zhong. All rights reserved.
//

#include "MCW_functions.h"

void MCW_reduce_in_place(int imin, int imax, int istep, int q,
                     const std::vector<size_t> & js,
                     std::vector<size_t> & js_red,
                     const std::vector< std::vector<ldouble> > & S,
                     const std::vector< std::vector<size_t> > & J,
                     const std::vector< std::vector<ldouble> > & sum_x,
                     const std::vector< ldouble > & sum_x_sq,
                     const std::vector< std::vector<ldouble> > & sum_w)
{
  int N = (imax - imin) / istep + 1;

  js_red = js;

  if(N >= js.size()) {
    return;
  }

  // Two positions to move candidate j's back and forth
  int left = -1; // points to last favorable position / column
  int right = 0; // points to current position / column

  size_t m = js_red.size();

  while(m > N) { // js_reduced has more than N positions / columns

    int p = (left + 1);

    int i = (imin + p * istep);
    size_t j = (js_red[right]);
    ldouble Sl = (S[q-1][j-1] +
      MCW_dissimilarity(j, i, sum_x, sum_x_sq, sum_w));
    // ssq(j, i, sum_x, sum_x_sq, sum_w));

    size_t jplus1 = (js_red[right+1]);
    ldouble Slplus1 = (S[q-1][jplus1-1] +
      MCW_dissimilarity(jplus1, i, sum_x, sum_x_sq, sum_w));
    // ssq(jplus1, i, sum_x, sum_x_sq, sum_w));

    if(Sl < Slplus1 && p < N-1) {
      js_red[ ++ left ] = j; // i += istep;
      right ++; // move on to next position / column p+1
    } else if(Sl < Slplus1 && p == N-1) {
      js_red[ ++ right ] = j; // delete position / column p+1
      m --;
    } else { // (Sl >= Slplus1)
      if(p > 0) { // i > imin
        // delete position / column p and
        //   move back to previous position / column p-1:
        js_red[right] = js_red[left --];
        // p --; // i -= istep;
      } else {
        right ++; // delete position / column 0
      }
      m --;
    }
  }

  for(int r=(left+1); r < m; ++r) {
    js_red[r] = js_red[right++];
  }

  js_red.resize(m);
  return;
}

inline void MCW_fill_even_positions
  (int imin, int imax, int istep, int q,
   const std::vector<size_t> & js,
   std::vector< std::vector<ldouble> > & S,
   std::vector< std::vector<size_t> > & J,
   const std::vector< std::vector<ldouble> > & sum_x,
   const std::vector< ldouble > & sum_x_sq,
   const std::vector< std::vector<ldouble> > & sum_w)
{
  // Derive j for even rows (0-based)
  size_t n = (js.size());
  int istepx2 = (istep << 1);
  size_t jl = (js[0]);
  for(int i=(imin), r(0); i<=imax; i+=istepx2) {

    // auto jmin = (i == imin) ? js[0] : J[q][i - istep];

    while(js[r] < jl) {
      // Increase r until it points to a value of at least jmin
      r ++;
    }

    // Initialize S[q][i] and J[q][i]
    S[q][i] = S[q-1][js[r]-1] +
      MCW_dissimilarity(js[r], i, sum_x, sum_x_sq, sum_w);
    // ssq(js[r], i, sum_x, sum_x_sq, sum_w);
    J[q][i] = js[r]; // rmin

    // Look for minimum S upto jmax within js
    int jh = (i + istep <= imax)
      ? J[q][i + istep] : js[n-1];

    int jmax = std::min((int)jh, (int)i);

    ldouble sjimin(
        MCW_dissimilarity(jmax, i, sum_x, sum_x_sq, sum_w)
      // ssq(jmax, i, sum_x, sum_x_sq, sum_w)
    );

    for(++ r; r < n && js[r]<=jmax; r++) {

      const size_t & jabs = js[r];

      if(jabs > i) break;

      if(jabs < J[q-1][i]) continue;

      ldouble s =
        MCW_dissimilarity(jabs, i, sum_x, sum_x_sq, sum_w);
      // (ssq(jabs, i, sum_x, sum_x_sq, sum_w));
      ldouble Sj = (S[q-1][jabs-1] + s);

      if(Sj <= S[q][i]) {
        S[q][i] = Sj;
        J[q][i] = js[r];
      } else if(S[q-1][jabs-1] + sjimin > S[q][i]) {
        break;
      } /*else if(S[q-1][js[rmin]-1] + s > S[q][i]) {
 break;
    } */
    }
    r --;
    jl = jh;
  }
}

inline void MCW_find_min_from_candidates
  (int imin, int imax, int istep, int q,
   const std::vector<size_t> & js,
   std::vector< std::vector<ldouble> > & S,
   std::vector< std::vector<size_t> > & J,
   const std::vector< std::vector<ldouble> > & sum_x,
   const std::vector< ldouble > & sum_x_sq,
   const std::vector< std::vector<ldouble> > & sum_w)
{
  size_t rmin_prev = (0);

  for(int i=(imin); i<=imax; i+=istep) {

    size_t rmin = (rmin_prev);

    // Initialization of S[q][i] and J[q][i]
    S[q][i] = S[q-1][js[rmin] - 1] +
      MCW_dissimilarity(js[rmin], i, sum_x, sum_x_sq, sum_w);
    // ssq(js[rmin], i, sum_x, sum_x_sq, sum_w);
    J[q][i] = js[rmin];

    for(size_t r = (rmin+1); r<js.size(); ++r) {

      const size_t & j_abs = (js[r]);

      if(j_abs < J[q-1][i]) continue;
      if(j_abs > i) break;

      ldouble Sj = (S[q-1][j_abs - 1] +
        MCW_dissimilarity(j_abs, i, sum_x, sum_x_sq, sum_w));
      // ssq(j_abs, i, sum_x, sum_x_sq, sum_w));
      if(Sj <= S[q][i]) {
        S[q][i] = Sj;
        J[q][i] = js[r];
        rmin_prev = r;
      }
    }
  }
}

void MCW_SMAWK
  (int imin, int imax, int istep, int q,
   const std::vector<size_t> & js,
   std::vector< std::vector<ldouble> > & S,
   std::vector< std::vector<size_t> > & J,
   const std::vector< std::vector<ldouble> > & sum_x,
   const std::vector< ldouble > & sum_x_sq,
   const std::vector< std::vector<ldouble> > & sum_w)
{
#ifdef DEBUG_REDUCE
  std::cout << "i:" << '[' << imin << ',' << imax << ']' << '+' << istep
            << std::endl;
#endif

  if(imax - imin <= 0 * istep) { // base case only one element left

    MCW_find_min_from_candidates(
      imin, imax, istep, q, js, S, J, sum_x, sum_x_sq, sum_w);

  } else {

    // REDUCE

#ifdef DEBUG //_REDUCE
    std::cout << "js:";
    for (size_t l=0; l < js.size(); ++l) {
      std::cout << js[l] << ",";
    }
    std::cout << std::endl;
    std::cout << std::endl;
#endif

    std::vector<size_t> js_odd;

    MCW_reduce_in_place(imin, imax, istep, q, js, js_odd,
                    S, J, sum_x, sum_x_sq, sum_w);

    int istepx2 = (istep << 1);
    int imin_odd = (imin + istep);
    int imax_odd = (imin_odd + (imax - imin_odd) / istepx2 * istepx2);

    // Recursion on odd rows (0-based):
    MCW_SMAWK(imin_odd, imax_odd, istepx2,q, js_odd,
          S, J, sum_x, sum_x_sq, sum_w);

#ifdef DEBUG // _REDUCE
    std::cout << "js_odd (reduced):";
    for (size_t l=0; l<js_odd.size(); ++l) {
      std::cout << js_odd[l] << ",";
    }
    std::cout << std::endl << std::endl;

    std::cout << "even pos:";
    for (int i=imin; i<imax; i+=istepx2) {
      std::cout << i << ",";
    }
    std::cout << std::endl << std::endl;

#endif

    MCW_fill_even_positions(imin, imax, istep, q, js,
                        S, J, sum_x, sum_x_sq, sum_w);
  }
}

void MCW_fill_row_q_SMAWK(int imin, int imax, int q,
                      std::vector< std::vector<ldouble> > & S,
                      std::vector< std::vector<size_t> > & J,
                      const std::vector< std::vector<ldouble> > & sum_x,
                      const std::vector< ldouble > & sum_x_sq,
                      const std::vector< std::vector<ldouble> > & sum_w)
{
  // Assumption: each cluster must have at least one point.

  std::vector<size_t> js(imax-q+1);
  int abs = (q);
  std::generate(js.begin(), js.end(), [&] { return abs++; } );

  MCW_SMAWK(imin, imax, 1, q, js, S, J, sum_x, sum_x_sq, sum_w);
}
