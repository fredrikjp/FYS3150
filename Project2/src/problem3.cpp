#include "problem3.hpp"

double max_offdiag_symmetric(const arma::mat &A, int &k, int &l) {
  int n = A.n_cols;
  double a = 0.0;

  for (int j = 1; j < n ; j++) {
    for (int i = 0; i < j; i++) {
      if (fabs(A(i, j)) > a) {
        a = fabs(A(i, j));
        k = i;
        l = j;
      }
    }
  }
  return a;
}

