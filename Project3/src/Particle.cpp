#include "Particle.hpp"
#include <armadillo>

// Constructor
Particle::Particle(double q, double m, arma::vec r, arma::vec v) {
  q_ = q;    // charge
  m_ = m;    // mass
  r_ = r; // position
  v_ = v; // velocity
}
