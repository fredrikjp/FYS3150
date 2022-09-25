#include "problem2.hpp"

int main()
{
arma::mat X = create_symmetric_tridiagonal(6,-1,2);
std::cout << X << std::endl;

arma::vec eigenval1;
arma::mat eigenvec1;
symmetric_tridiagonal_analytic_eigen_solver(eigenval1, eigenvec1, X);
test_symmetric_matrix_eigenequation(eigenval1, eigenvec1, X, pow(10, -15));
return 0;
} 
