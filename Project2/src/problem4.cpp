#include "problem4.hpp"
#include "problem3.hpp"
#include <iostream>
#include <ostream>

// Performs a single Jacobi rotation, to "rotate away"
// the off-diagonal element at A(k,l).
// - Assumes symmetric matrix, so we only consider k < l
// - Modifies the input matrices A and R
void jacobi_rotate(arma::mat &A, arma::mat &R, int k, int l) {
  // Compute s (sin(thata)) and c (cos(theta))
  double tau = (A(l, l) - A(k, k)) / (2 * A(k, l));
  double t;
  if (tau > 0) {
    t = 1 / (tau + std::sqrt(1 + std::pow(tau, 2)));
  } else {
    t = -1 / (-tau + std::sqrt(1 + std::pow(tau, 2)));
  }
  double c = 1 / (std::sqrt(1 + std::pow(t, 2)));
  double s = c * t;

  // Update matrix A
  double a = A(k, k);
  int n = A.n_cols;

  A(k, k) =
      A(k, k) * std::pow(c, 2) - 2 * A(k, l) * c * s + A(l, l) * std::pow(s, 2);
  A(l, l) = A(l, l) * std::pow(c, 2) + 2 * A(k, l) * c * s + a * std::pow(s, 2);
  A(k, l) = 0.;
  A(l, k) = 0.;

  for (int i = 0; i < n; i++) {
    if (i == k || i == l) {
    } else {
      double b = A(i, k);

      A(i, k) = A(i, k) * c - A(i, l) * s;
      A(k, i) = A(i, k);
      A(i, l) = A(i, l) * c + b * s;
      A(l, i) = A(i, l);
    }
  }
  // Update matrix R
  for (int i = 0; i < n; i++) {
    double r = R(i, k);

    R(i, k) = R(i, k) * c - R(i, l) * s;
    R(i, l) = R(i, l) * c + r * s;
  }
}

// Jacobi method eigensolver:
// - Runs jacobo_rotate until max off-diagonal element < eps
// - Writes the eigenvalues as entries in the vector "eigenvalues"
// - Writes the eigenvectors as columns in the matrix "eigenvectors"
//   (The returned eigenvalues and eigenvectors are sorted using
//   arma::sort_index)
// - Stops if it the number of iterations reaches "maxiter"
// - Writes the number of iterations to the integer "iterations"
// - Sets the bool reference "converged" to true if convergence was reached
// before hitting maxiter
void jacobi_eigensolver(const arma::mat &A, double eps, arma::vec &eigenvalues,
                        arma::mat &eigenvectors, const int maxiter,
                        int &iterations, bool &converged) {
  int n = A.n_cols;
  int k;
  int l;

  arma::mat A_new = A;
  arma::mat R = arma::mat(n, n, arma::fill::eye);
  eigenvalues = arma::vec(n);
  eigenvectors = arma::mat(n, n);
  converged = false;

  max_offdiag_symmetric(A_new, k, l);
  int i = 0;
  while (i < maxiter) {
    jacobi_rotate(A_new, R, k, l);
    max_offdiag_symmetric(A_new, k, l);
    i++;
    if (fabs(A_new(k, l)) <= eps) {
      converged = true;
      iterations = i;
      break;
    }
  }
  if (converged == false) {
    iterations = maxiter;
  }

  for (int i = 0; i < n; i++) {
    eigenvalues(i) = A_new(i, i);
    for (int j = 0; j < n; j++) {
      eigenvectors(i, j) = R(i, j);
    }
  }
  arma::uvec val_index = arma::sort_index(eigenvalues); 

  eigenvalues = eigenvalues.elem(val_index); 
  eigenvectors = eigenvectors.cols(val_index); 
}
