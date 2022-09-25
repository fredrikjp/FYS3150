#include "problem2.hpp"
#include "problem4.hpp"
#include <iomanip>

int main() {

  std::string filename = "problem5.txt";

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
  for (int i = 1; i <= 20; i++) {
    int n = 5 * i;
    double h = 1. / n;
    arma::mat A = create_symmetric_tridiagonal(n, -1. / std::pow(h, 2),
                                               2. / std::pow(h, 2));

    // Running once before timing to increase consistency of results
    jacobi_eigensolver(A, eps, eigenvalues, eigenvectors, maxiter, iterations,
                       converged);

    for (int i = 1; i <= 10; i++) {
      // Start measuring time
      auto t1 = std::chrono::high_resolution_clock::now();

      //
      // The code you want to perform timing on

      jacobi_eigensolver(A, eps, eigenvalues, eigenvectors, maxiter, iterations,
                         converged);

      //

      // Stop measuring time
      auto t2 = std::chrono::high_resolution_clock::now();
      if (converged == false) {
        std::cout << "Solution not converged for N=" << n << std::endl;
        break;
      }
      // Calculate the elapsed time
      // We use chrono::duration<double>::count(), which by default returns
      // duration in seconds
      double duration_seconds = std::chrono::duration<double>(t2 - t1).count();

      ofile << std::setw(width) << std::setprecision(prec) << std::scientific
            << n << std::setw(width) << std::setprecision(prec)
            << std::scientific << duration_seconds << std::endl;
    }
  }
  ofile.close();
  return 0;
}
