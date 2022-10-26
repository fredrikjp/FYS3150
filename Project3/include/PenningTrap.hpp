#include "Particle.hpp"
#include <armadillo>
#include <iostream>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>

class PenningTrap {
private:
  double ke_ =
      1.38935333 * pow(10, 5); // Coulomb constant [u(mu m)^3/((mu s)^2 e)]

public:
  double B0_, V0_, d_, f_, w_;
  std::vector<Particle> particles_;
  bool particle_interaction_ = false;
  bool time_dependent_E_ = false;
  
  // Simulation matricies storing positions R_ and velocities V_
  arma::mat R_;
  arma::mat V_;

  // Simulation time
  std::vector<double> t_;

  // Constructor
  PenningTrap(double B0_in, double V0_in, double d_in);

  // Constructor
  PenningTrap(double B0_in, double V0_in, double d_in, double f_in, double w_in);

  // Add a particle to the trap
  void add_particle(Particle p_in);

  // Add n random particles to the trap
  void add_random_particles(int n, int q, double m);

  // Count of particles in trap
  int particles_in_trap();

  // External electric field at point r=(x,y,z)
  arma::vec external_E_field(arma::vec r);

  // External magnetic field at point r=(x,y,z)
  arma::vec external_B_field(arma::vec r);

  // Force on particle_i from particle_j
  arma::vec force_particle(int i, int j);

  // The total force on particle_i from the external fields
  arma::vec total_force_external(int i);

  // The total force on particle_i from the other particles
  arma::vec total_force_particles(int i);

  // The total force on particle_i from both external fields and other particles
  arma::vec total_force(int i);

  // Evolve the system one time step (dt) using Runge-Kutta 4th order
  void evolve_RK4(double dt);

  // Evolve the system one time step (dt) using Forward Euler
  void evolve_forward_Euler(double dt);

  // Simulate the system for n timesteps dt
  void simulation(double dt, int n, std::string method);

  // Write to file
  void write_to_file(std::string filename);
};
