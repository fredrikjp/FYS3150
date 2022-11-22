#include "../include/MCMC.hpp"
using namespace std;

// Constructor
MCMC::MCMC(int L, double T, bool ordered) {
  L_ = L;
  N_ = pow(L, 2);
  T_ = T;
  beta_ = 1 / (kB_ * T_);

  // Calculate Boltzmann factor
  arma::vec Delta_E = 2 * arma::vec("-4 -3 -2 -1 0 1 2 3 4");
  Boltz_ = exp(-beta_ * Delta_E);
  // unsigned int seed = 1840371;
  // Get a seed from the system clock
  unsigned int seed =
      std::chrono::system_clock::now().time_since_epoch().count();
  generator_.seed(seed);

  // Generate initial random state
  MCMC::random_state(L_, ordered);
}
void MCMC::random_state(int L, bool ordered) {
  uniform_int_distribution<int> dist(0, 1);
  // Surrounding the actual state matrix in a square with zeros to simplify code
  s_ = arma::mat(L + 2, L + 2).fill(0);
  for (int i = 1; i < L + 1; ++i) {
    for (int j = 1; j < L + 1; j++) {
      if (ordered) {
        s_(i, j) = 1;
      } else {
        s_(i, j) = (dist(generator_) == 1) ? 1 : -1;
      }
    }
  }
}

void MCMC::flip_spin() {
  uniform_int_distribution<int> rand_index(1, L_);
  int i = rand_index(generator_);
  int j = rand_index(generator_);
  s_(i, j) *= -1; // Flip random spin element
  int Boltz_index =
      4 -
      (s_(i + 1, j) + s_(i - 1, j) + s_(i, j - 1) + s_(i, j + 1)) * s_(i, j);
  double rel_prob =
      Boltz_(Boltz_index); // Relative probability p(s_after)/p(s_before) =
                           // Boltzmann factor
  uniform_real_distribution<double> rand_r(0, 1);
  double r = rand_r(generator_);
  // Trying to avoid if-test (not sure if using bools is beneficial)
  double x = rel_prob - r;
  s_(i, j) *= ((x > 0) - (x <= 0)); // Reverting to previous s if r >= rel_prob
}

void MCMC::Monte_Carlo(int n, int burn_in) {
  for (int burn = 0; burn < N_ * burn_in; burn++) {
      MCMC::flip_spin();
  }
  for (int a = 0; a < n; a++) {
    for (int b = 0; b < N_; b++) {
      MCMC::flip_spin();
      double M = arma::accu(s_);
      double E = 0;
      for (int i = 1; i < L_ + 1; i++) {
        for (int j = 1; j < L_ + 1; j++) {
          E += -s_(i, j) * s_(i + 1, j);
          E += -s_(i, j) * s_(i + 1, j);
        }
      }
      m_.push_back(M / N_);
      e_.push_back(E / N_);
    }
    arma::vec m = arma::conv_to<arma::vec>::from(m_);
    arma::vec e = arma::conv_to<arma::vec>::from(e_);

    mean_e_.push_back(arma::mean(e));
    mean_m_.push_back(arma::mean(arma::abs(m)));

    CV_.push_back(
        1 / (N_ * kB_ * T_ * T_) *
        (arma::mean(arma::square(e * N_)) - pow(mean_e_.back() * N_, 2)));
    chi_.push_back(
        1 / (N_ * kB_ * T_) *
        (arma::mean(arma::square(m * N_)) - pow(mean_m_.back() * N_, 2)));
  }
}

void MCMC::write_to_file(string filename) {
  // Create and open the output file. Or, technically, create
  // an "output file stream" (type std::ofstream) and connect it to our
  // filename.
  std::ofstream ofile;
  ofile.open(filename);

  // Some width and precision parameters we will use to format the output
  int width = 16;
  int prec = 8;
  int n = mean_e_.size();
  for (int i = 0; i < n; i++) {
    ofile << std::setw(width) << std::setprecision(prec) << std::scientific
          << mean_e_[i] << std::setw(width) << std::setprecision(prec)
          << std::scientific << mean_m_[i] << std::setw(width)
          << std::setprecision(prec) << std::scientific << CV_[i]
          << std::setw(width) << std::setprecision(prec) << std::scientific
          << chi_[i] << endl;
  }
  ofile.close();
}

void MCMC::write_m_e_to_file(string filename) {
  // Create and open the output file. Or, technically, create
  // an "output file stream" (type std::ofstream) and connect it to our
  // filename.
  std::ofstream ofile;
  ofile.open(filename);

  // Some width and precision parameters we will use to format the output
  int width = 16;
  int prec = 8;
  int n = e_.size();
  for (int i = 0; i < n; i++) {
    ofile << std::setw(width) << std::setprecision(prec) << std::scientific
          << e_[i] << std::setw(width) << std::setprecision(prec)
          << std::scientific << m_[i] << endl;
  }
  ofile.close();
}
