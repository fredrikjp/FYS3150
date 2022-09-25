#ifndef __problem2_hpp__  
#define __problem2_hpp__


#include <armadillo>
#include <cmath>
#include <iostream>
#include <ostream>
#include <vector>


arma::mat create_tridiagonal(const arma::vec& a, const arma::vec& d, const arma::vec& e);

arma::mat create_tridiagonal(int n, double a, double d, double e);

arma::mat create_symmetric_tridiagonal(int n, double a, double d);

void symmetric_tridiagonal_analytic_eigen_solver(arma::vec& eigenval, arma::mat& eigenvec, arma::mat A);
	
void test_symmetric_matrix_eigenequation(arma::vec eigenval1, arma::mat eigenvec1, arma::mat X, double prec);
	
#endif
