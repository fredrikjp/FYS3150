
#ifndef __double_slit_hpp__
#define __double_slit_hpp__

#include <armadillo>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <omp.h>
#include <random>
#include <vector>
using namespace std;

class double_slit {

public:
  // Number of positions in both x and y direction M_ and number of time steps
  // N_
  int M_, N_;

  // Position step in x and y direction h_, and time step dt_
  double h_, dt_;

  // Complex matricies for Crank-Nicolson
  arma::sp_cx_mat A_;
  arma::sp_cx_mat B_;

  // Complex solution matrix
  arma::cx_mat U_;

  // Potential matrix
  arma::mat V_;

  // Constructor
  double_slit(double h, double dt, double t, double cx, double cy, double dx,
              double dy, double px, double py, double v0, double wall_width,
              double wall_postion, double slit_separation, double slit_opening,
              int n_slits);

  // Transform position indecies
  int index_transform(int i, int j);

  // Generates Crank-Nicolson matricies 
  void generate_matricies();

  // Creates initial wave function state
  void generate_initial_state(double cx, double cy, double dx, double dy,
                              double px, double py);

  // Sets up potential
  void generate_potential(double v0, double wall_width,
                               double wall_postion, double slit_seperation,
                               double slit_opening, int n_slits);

  // Solves the Crank-Nicolson scheme
  void solve_Crank_Nicolson();

  // Writes time vector t and solution matrix U_ to file
  void write_to_file(string filename, int prec);

};

#endif
