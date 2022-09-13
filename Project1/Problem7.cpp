#include "Problem6.cpp"
#include <math.h>
#include <string>
#include <tuple>

int main() {

  for (int i = 1; i <= 4; i++) {
    int n = pow(10, i);
    double h = 1.0 / (n - 1);
    arma::vec a = arma::vec(n - 3).fill(-1.);
    arma::vec b = arma::vec(n - 2).fill(2.);
    arma::vec c = arma::vec(n - 3).fill(-1.);
    arma::vec x = arma::linspace(h, 1 - h, n - 2);
    arma::vec g = 100 * exp(-10 * x) * pow(h, 2);
    arma::vec v = Thomas(a, b, c, g);
    std::string filename = "Problem7_n" + std::to_string(n) + ".txt";

    std::ofstream ofile;
    ofile.open(filename);

    int width = 16;
    int prec = 8;

    ofile << std::setw(width) << std::setprecision(prec) << std::scientific
          << 0 // Initial bondary
          << std::setw(width) << std::setprecision(prec) << std::scientific << 0
          << std::endl;

    for (int i = 0; i < n - 2; i++) {
      ofile << std::setw(width) << std::setprecision(prec) << std::scientific
            << x(i) << std::setw(width) << std::setprecision(prec)
            << std::scientific << v(i) << std::endl;
    }

    ofile << std::setw(width) << std::setprecision(prec) << std::scientific
          << 1 // End bondary
          << std::setw(width) << std::setprecision(prec) << std::scientific << 0
          << std::endl;

    ofile.close();
  }

  return 0;
}
