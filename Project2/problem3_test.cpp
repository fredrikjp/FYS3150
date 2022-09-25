#include "problem3.hpp"

int main() {
  arma::mat A = arma::mat(4, 4, arma::fill::eye);
  A(3, 0) = 0.5;
  A(0, 3) = 0.5;
  A(2, 1) = -0.7;
  A(1, 2) = -0.7;
  arma::cout << A << std::endl;

  int row;
  int col;
  std::cout << max_offdiag_symmetric(A, row, col) << " row " << row
            << " column " << col << std::endl;
  return 0;
}
