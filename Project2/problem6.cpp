#include "problem2.hpp"
#include "problem4.hpp"
#include <iomanip>

int main() {
  for (int a = 1; a <= 2; a++) {

    int n = pow(10, a) - 1;
    std::string filename = "problem6_n" + std::to_string(n) + ".txt";
    std::ofstream ofile;
    ofile.open(filename);

    int width = 16;
    int prec = 8;

    double eps = 0.1;
    arma::vec eigenvalues;
    arma::mat eigenvectors;
    int maxiter = 100000;
    int iterations;
    bool converged;
    double h = 1. / (n + 1);
    arma::mat A = create_symmetric_tridiagonal(n, -1. / std::pow(h, 2),
                                               2. / std::pow(h, 2));

    jacobi_eigensolver(A, eps, eigenvalues, eigenvectors, maxiter, iterations,
                       converged);

    arma::vec eigenval;
    arma::mat eigenvec;
    symmetric_tridiagonal_analytic_eigen_solver(eigenval, eigenvec, A);

    arma::vec x = arma::linspace(0, n + 1, n + 2) * h;
    int v_start = 0;
    int v_end = 0;

    for (int j = 0; j < 3; j++) {
      ofile << std::setw(width) << std::setprecision(prec) << std::scientific
            << x(0) << std::setw(width) << std::setprecision(prec)
            << std::scientific << v_start << std::setw(width)
            << std::setprecision(prec) << std::scientific << v_start
            << std::endl;
      for (int i = 0; i < n; i++) {
        double x_i = (i + 1) * h;
        ofile << std::setw(width) << std::setprecision(prec) << std::scientific
              << x(i + 1) << std::setw(width) << std::setprecision(prec)
              << std::scientific << eigenvectors.col(j)(i) << std::setw(width)
              << std::setprecision(prec) << std::scientific
              << eigenvec.col(j)(i) << std::endl;
      }
      ofile << std::setw(width) << std::setprecision(prec) << std::scientific
            << x(n + 1) << std::setw(width) << std::setprecision(prec)
            << std::scientific << v_end << std::setw(width)
            << std::setprecision(prec) << std::scientific << v_end << std::endl;
    }
    ofile.close();
  }

  return 0;
}
