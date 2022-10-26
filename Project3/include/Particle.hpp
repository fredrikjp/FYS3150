#ifndef __Particle_hpp__
#define __Particle_hpp__

#include <armadillo>

class Particle {
public:
  double q_, m_;    // Charge and mass
  arma::vec r_, v_; // Position and velocity vector

  // Constructor
  Particle(double q, double m, arma::vec r, arma::vec v);

};

#endif
