#include "../include/PenningTrap.hpp"
#include <algorithm>
#include <functional>
#include <iostream>
#include <iterator>
#include <math.h>
#include <ostream>
#include <string>

// Constructor
PenningTrap::PenningTrap(double B0_in, double V0_in, double d_in) {
  B0_ = B0_in;
  V0_ = V0_in;
  d_ = d_in;
}

// Constructor for time dependent electric field
PenningTrap::PenningTrap(double B0_in, double V0_in, double d_in, double f_in,
                         double w_in) {
  B0_ = B0_in;
  V0_ = V0_in;
  d_ = d_in;
  f_ = f_in;
  w_ = w_in;
  time_dependent_E_ = true;
}
// Add a particle to the trap
void PenningTrap::add_particle(Particle p_in) { particles_.push_back(p_in); }

// Add n randomly generated particles of charge q and mass m to the trap
void PenningTrap::add_random_particles(int n, int q, double m) {
  for (int i = 0; i < n; i++) {
    arma::vec r = arma::vec(3).randn() * 0.1 * d_; // random initial position
    arma::vec v = arma::vec(3).randn() * 0.1 * d_; // random initial velocity
    particles_.push_back(Particle(q, m, r, v));
  }
}

// Count of particles in trap
int PenningTrap::particles_in_trap() {
  int count = 0;
  for (int i = 0; i < particles_.size(); i++) {
    arma::vec R = particles_[i].r_;
    if (arma::norm(R) < d_) {
      count += 1;
    }
  }
  return count;
}

// External electric field at point r=(x,y,z)
arma::vec PenningTrap::external_E_field(arma::vec r) {
  arma::vec E_field = arma::vec(r.n_elem).fill(0);
  double V0;
  if (time_dependent_E_) {
    V0 = V0_ * (1 + f_ * cos(w_ * t_.back()));
  } else {
    V0 = V0_;
  }
  if (arma::norm(r) < d_) {
    E_field(0) = V0 / pow(d_, 2) * r(0);
    E_field(1) = V0 / pow(d_, 2) * r(1);
    E_field(2) = -2 * V0 / pow(d_, 2) * r(2);
  }
  return E_field;
}

// External magnetic field at point r=(x,y,z)
arma::vec PenningTrap::external_B_field(arma::vec r) {
  arma::vec B_field = arma::vec(r.n_elem).fill(0.);
  if (arma::norm(r) < d_) {
    B_field(2) = B0_;
  }
  return B_field;
}

// Force on particle_i from particle_j
arma::vec PenningTrap::force_particle(int i, int j) {
  Particle particle_i = particles_[i];
  Particle particle_j = particles_[j];
  double q_i = particle_i.q_;
  double q_j = particle_j.q_;
  arma::vec r_i = particle_i.r_;
  arma::vec r_j = particle_j.r_;
  arma::vec force =
      ke_ * q_i * q_j * (r_i - r_j) / pow(arma::norm(r_i - r_j), 3);
  return force;
}

// The total force on particle_i from the external fields
arma::vec PenningTrap::total_force_external(int i) {
  Particle particle = particles_[i];
  double q = particle.q_;
  arma::vec r = particle.r_;
  arma::vec v = particle.v_;
  arma::vec force = q * PenningTrap::external_E_field(r) +
                    arma::cross(q * v, PenningTrap::external_B_field(r));

  return force;
}

// The total force on particle_i from the other particles
arma::vec PenningTrap::total_force_particles(int i) {
  arma::vec force = arma::vec(3).fill(0);
  for (int j = 0; j < particles_.size(); j++) {
    if (i != j) {
      force += PenningTrap::force_particle(i, j);
    }
  }
  return force;
}

// The total force on particle_i from both external fields and other particles_
arma::vec PenningTrap::total_force(int i) {
  arma::vec force;
  if (particle_interaction_) {
    force = PenningTrap::total_force_particles(i) +
            PenningTrap::total_force_external(i);
  } else {
    force = PenningTrap::total_force_external(i);
  }
  return force;
}

// Evolve the system one time step (dt) using Runge-Kutta 4th order
void PenningTrap::evolve_RK4(double dt) {
  double m;
  arma::vec v;
  arma::vec r;
  arma::vec kr1;
  arma::vec kv1;
  arma::vec kr2;
  arma::vec kv2;
  arma::vec kr3;
  arma::vec kv3;
  arma::vec kr4;
  arma::vec kv4;
  for (int i = 0; i < particles_.size(); i++) {
    m = particles_[i].m_;
    v = particles_[i].v_;
    r = particles_[i].r_;

    // Time index for last element
    int j = t_.size() - 1;

    kr1 = dt * v;
    kv1 = dt * PenningTrap::total_force(i) / m;

    particles_[i].r_ = r + kr1 / 2.;
    particles_[i].v_ = v + kv1 / 2.;

    // Temporary time for total_force
    t_.push_back(t_[j] + dt / 2.);

    kr2 = dt * (v + kv1 / 2.);
    kv2 = dt * PenningTrap::total_force(i) / m;

    particles_[i].r_ = r + kr2 / 2.;
    particles_[i].v_ = v + kv2 / 2.;

    kr3 = dt * (v + kv2 / 2.);
    kv3 = dt * PenningTrap::total_force(i) / m;

    particles_[i].r_ = r + kr3;
    particles_[i].v_ = v + kv3;

    t_[j+1] = t_[j] + dt;

    kr4 = dt * (v + kv3);
    kv4 = dt * PenningTrap::total_force(i) / m;

    particles_[i].r_ = r + 1. / 6 * (kr1 + 2 * kr2 + 2 * kr3 + kr4);
    particles_[i].v_ = v + 1. / 6 * (kv1 + 2 * kv2 + 2 * kv3 + kv4);

    // Delete temporary time element
    t_.pop_back();
  }
}

// Evolve the system one time step (dt) using Forward Euler
void PenningTrap::evolve_forward_Euler(double dt) {
  for (int i = 0; i < particles_.size(); i++) {
    double m = particles_[i].m_;
    arma::vec v = particles_[i].v_;
    arma::vec r = particles_[i].r_;

    arma::vec r_next = r + dt * v;
    arma::vec v_next = v + dt * PenningTrap::total_force(i) / m;

    particles_[i].v_ = v_next;
    particles_[i].r_ = r_next;
  }
}

void PenningTrap::simulation(double dt, int n, std::string method) {
  int m = 3 * particles_.size();

  // Positions R_ and velocities V_ for all particles at all times t_
  R_ = arma::mat(m, n + 1);
  V_ = arma::mat(m, n + 1);
  // Initial conditions
  t_.push_back(0);
  for (int p = 0; p < particles_.size(); p++) {
    for (int j = 0; j < 3; j++) {
      R_(3 * p + j, 0) = particles_[p].r_(j);
      V_(3 * p + j, 0) = particles_[p].v_(j);
    }
  }
  // Rest of timesteps
  for (int i = 1; i <= n; i++) {
    if (method == "RK4") {
      PenningTrap::evolve_RK4(dt);
    } else if (method == "Forward Euler") {
      PenningTrap::evolve_forward_Euler(dt);
    } else {
      throw std::invalid_argument("Must use method 'Forward Euler' or 'RK4'");
    }
    t_.push_back(dt * i);
    for (int p = 0; p < particles_.size(); p++) {
      for (int j = 0; j < 3; j++) {
        R_(3 * p + j, i) = particles_[p].r_(j);
        V_(3 * p + j, i) = particles_[p].v_(j);
      }
    }
  }
}

void PenningTrap::write_to_file(std::string filename) {

  // Create and open the output file. Or, technically, create
  // an "output file stream" (type std::ofstream) and connect it to our
  // filename.
  std::ofstream ofile;
  ofile.open(filename);

  // Some width and precision parameters we will use to format the output
  int width = 28;
  int prec = 20;

  for (int i = 0; i < R_.n_rows; i++) {
    for (int j = 0; j < R_.n_cols; j++) {
      // Write a line with the current x and y values (nicely formatted) to file
      ofile << std::setw(width) << std::setprecision(prec) << std::scientific
            << R_(i, j);
    }
    ofile << std::endl;
  }

  ofile << std::endl;
  for (int i = 0; i < R_.n_rows; i++) {
    for (int j = 0; j < R_.n_cols; j++) {
      // Write a line with the current x and y values (nicely formatted) to file
      ofile << std::setw(width) << std::setprecision(prec) << std::scientific
            << V_(i, j);
    }
    ofile << std::endl;
  }

  ofile << std::endl;
  for (int i = 0; i < R_.n_cols; i++) {
    // Write a line with the current x and y values (nicely formatted) to file
    ofile << std::setw(width) << std::setprecision(prec) << std::scientific
          << t_[i];
  }
  ofile.close();
}
