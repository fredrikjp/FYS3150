#include "Problem6.cpp"
#include "Problem9.cpp"
#include <chrono>
#include <cmath>

int main() {

  std::string filename = "Thomas_duration.txt";

  std::ofstream ofile;
  ofile.open(filename);

  std::string filename_ = "Special_duration.txt";

  std::ofstream ofile_;
  ofile_.open(filename_);

  for (int i = 1; i <= 6; i++) {
    int n = pow(10, i);
    double h = 1.0 / (n - 1);
    arma::vec a = arma::vec(n - 3).fill(-1.);
    arma::vec b = arma::vec(n - 2).fill(2.);
    arma::vec c = arma::vec(n - 3).fill(-1.);
    arma::vec x = arma::linspace(h, 1 - h, n - 2);
    arma::vec g = 100 * exp(-10 * x) * pow(h, 2);

    for (int i = 1; i <= 10; i++) {
      // Start measuring time
      auto t1 = std::chrono::high_resolution_clock::now();

      //
      // The code you want to perform timing on
      arma::vec v = Thomas(a, b, c, g);
      //

      // Stop measuring time
      auto t2 = std::chrono::high_resolution_clock::now();

      // Calculate the elapsed time
      // We use chrono::duration<double>::count(), which by default returns
      // duration in seconds
      double duration_seconds = std::chrono::duration<double>(t2 - t1).count();

      int width = 16;
      int prec = 8;

      ofile << std::setw(width) << std::setprecision(prec) << std::scientific
            << n << std::setw(width) << std::setprecision(prec)
            << std::scientific << duration_seconds << std::endl;

      // Start measuring time
      auto t3 = std::chrono::high_resolution_clock::now();

      //
      // The code you want to perform timing on
      arma::vec v_ = Special(a, b, c, g);
      //

      // Stop measuring time
      auto t4 = std::chrono::high_resolution_clock::now();

      // Calculate the elapsed time
      // We use chrono::duration<double>::count(), which by default returns
      // duration in seconds
      double duration_seconds_ = std::chrono::duration<double>(t4 - t3).count();

      ofile_ << std::setw(width) << std::setprecision(prec) << std::scientific
             << n << std::setw(width) << std::setprecision(prec)
             << std::scientific << duration_seconds_ << std::endl;
    }
  }
  ofile.close();
  ofile_.close();
  return 0;
}
