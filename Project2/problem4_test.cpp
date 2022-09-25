#include "problem4.hpp"
#include "problem2.hpp"


int main()
{
int n = 6;
double h = 1./n;
arma::mat A = create_symmetric_tridiagonal(n, -1./std::pow(h, 2), 2./std::pow(h, 2));

double eps = std::pow(10, -100);
arma::vec eigenvalues;
arma::mat eigenvectors;
int maxiter = 100;
int iterations;
bool converged;

jacobi_eigensolver(A, eps, eigenvalues, eigenvectors, maxiter, iterations, converged);

arma::vec eigenval;
arma::mat eigenvec;

symmetric_tridiagonal_analytic_eigen_solver(eigenval, eigenvec, A);

double prec = pow(10, -10);
test_symmetric_matrix_eigenequation(eigenvalues, eigenvectors, A, prec);

return 0;
}
