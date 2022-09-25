#ifndef __problem4_hpp__  
#define __problem4_hpp__

#include <armadillo>

void jacobi_rotate(arma::mat& A, arma::mat& R, int k, int l);

void jacobi_eigensolver(const arma::mat& A, double eps, arma::vec& eigenvalues, arma::mat& eigenvectors, const int maxiter, int& iterations, bool& converged);

#endif
