
#ifndef __MCMC_hpp__
#define __MCMC_hpp__

#include <armadillo>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <omp.h>
#include <random>
#include <vector>
using namespace std;

class MCMC {
private:
  // Construct a Mersenne Twister 19937 random number generator
  mt19937 generator_;
  double kB_ = 1; //.380649e-23; // [J/K], J = J joules

public:
  int L_, N_;
  double T_, beta_;
  arma::vec Boltz_;

  // State matrix
  arma::mat s_;

  // Magnetization and energy per spin
  std::vector<double> m_, e_;
  std::vector<double> mean_m_, mean_e_; // mean of |m| and e
  std::vector<double> chi_, CV_;        // Susceptibility and heat capacity

  // Constructor
  MCMC(int L, double T, bool ordered);

  // Create random state
  void random_state(int L, bool ordered);

  // Try flip single spin
  void flip_spin();

  // Perform n Monte Carlo cycles and store quantities
  void Monte_Carlo(int n, int burn_in);

  // Parallelized model fro different temperatures
  void T_Parallelized(int n, int burn_in, int n_T_values);

  // Writes means, CV and chi to file
  void write_to_file(string filename);

  // Writes e and m to file
  void write_m_e_to_file(string filename);
};

#endif
