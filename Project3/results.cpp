#include "include/Particle.hpp"
#include "include/PenningTrap.hpp"

int main() {
    PenningTrap obj1 = PenningTrap(9.65 * pow(10, 1), 2.41 * pow(10, 6), 500.);
    PenningTrap obj2 = PenningTrap(9.65 * pow(10, 1), 2.41 * pow(10, 6), 500.);
    PenningTrap obj3 = PenningTrap(9.65 * pow(10, 1), 2.41 * pow(10, 6), 500.);
    PenningTrap obj4 = PenningTrap(9.65 * pow(10, 1), 2.41 * pow(10, 6), 500.);
    PenningTrap obj5 = PenningTrap(9.65 * pow(10, 1), 2.41 * pow(10, 6), 500.);
    PenningTrap obj6 = PenningTrap(9.65 * pow(10, 1), 2.41 * pow(10, 6), 500.);
    PenningTrap obj7 = PenningTrap(9.65 * pow(10, 1), 2.41 * pow(10, 6), 500.);
    PenningTrap obj8 = PenningTrap(9.65 * pow(10, 1), 2.41 * pow(10, 6), 500.);
    PenningTrap obj9 = PenningTrap(9.65 * pow(10, 1), 2.41 * pow(10, 6), 500.);
    PenningTrap obj10 = PenningTrap(9.65 * pow(10, 1), 2.41 * pow(10, 6), 500.);
    PenningTrap obj11 = PenningTrap(9.65 * pow(10, 1), 2.41 * pow(10, 6), 500.);

    arma::vec r1 = arma::vec("20. 0. 20.");
    arma::vec v1 = arma::vec("0. 25. 0.");
    Particle P1 = Particle(1., 40.078, r1, v1);

    arma::vec r2 = arma::vec("25. 25. 0.");
    arma::vec v2 = arma::vec("0. 40. 5.");
    Particle P2 = Particle(1., 40.078, r2, v2);

    obj1.add_particle(P1);

    obj2.add_particle(P1);
    obj2.add_particle(P2);

    obj3.add_particle(P1);
    obj3.add_particle(P2);

    obj4.add_particle(P1);

    obj5.add_particle(P1);

    obj6.add_particle(P1);

    obj7.add_particle(P1);

    obj8.add_particle(P1);

    obj9.add_particle(P1);

    obj10.add_particle(P1);

    obj11.add_particle(P1);




    obj1.simulation(0.01, 5000, "RK4");
    obj1.write_to_file("P1_t50_.txt");

    obj2.simulation(0.01, 5000, "RK4");
    obj2.write_to_file("P2_t50_PIfalse.txt");

    obj3.particle_interaction_ = true;
    obj3.simulation(0.01, 5000, "RK4");
    obj3.write_to_file("P2_t50_PItrue.txt");

    obj4.simulation(50./4000, 4000, "RK4");
    obj4.write_to_file("P1_RK4_t50_n4000.txt");

    obj5.simulation(50./8000, 8000, "RK4");
    obj5.write_to_file("P1_RK4_t50_n8000.txt");

    obj6.simulation(50./16000, 16000, "RK4");
    obj6.write_to_file("P1_RK4_t50_n16000.txt");

    obj7.simulation(50./32000, 32000, "RK4");
    obj7.write_to_file("P1_RK4_t50_n32000.txt");

    obj8.simulation(50./4000, 4000, "Forward Euler");
    obj8.write_to_file("P1_FE_t50_n4000.txt");

    obj9.simulation(50./8000, 8000, "Forward Euler");
    obj9.write_to_file("P1_FE_t50_n8000.txt");

    obj10.simulation(50./16000, 16000, "Forward Euler");
    obj10.write_to_file("P1_FE_t50_n16000.txt");

    obj11.simulation(50./32000, 32000, "Forward Euler");
    obj11.write_to_file("P1_FE_t50_n32000.txt");

arma::arma_rng::set_seed(0);
int n = 116;
arma::vec w = arma::linspace(0.2, 2.5, n);
arma::vec f = arma::vec("0.1 0.4 0.7");
arma::mat particles_left = arma::mat(3, n);
for (int i = 0; i < 3; i++) {
  for (int j = 0; j < n; j++) {
    PenningTrap obj = PenningTrap(9.25e+1, 2.41e+6, 500., f(i), w(j));
    obj.add_random_particles(100, 1, 40.078);
    obj.simulation(500. / 5000, 5000, "RK4");
    particles_left(i, j) = obj.particles_in_trap();
  }
}

std::ofstream ofile;
ofile.open("Particles_left.txt");
ofile << particles_left << std::endl << f.as_row() << std::endl << w.as_row();
ofile.close();

  int n1 = 20;
  arma::vec w1 = arma::linspace(1.94, 2.22, n1);
  arma::vec f1 = arma::vec("0.1 0.4 0.7");
  arma::mat particles_left1 = arma::mat(3, n1);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < n1; j++) {
      PenningTrap obj = PenningTrap(9.25e+1, 2.41e+6, 500., f1(i), w1(j));
      obj.add_random_particles(100, 1, 40.078);
      obj.simulation(500. / 5000, 5000, "RK4");
      particles_left1(i, j) = obj.particles_in_trap();
    }
  }

  std::ofstream ofile1;
  ofile1.open("Particles_left_Coulomb_off.txt");
  ofile1 << particles_left1 << std::endl
         << f1.as_row() << std::endl
         << w1.as_row();
  ofile1.close();

  int n2 = 20;
  arma::vec w2 = arma::linspace(1.94, 2.22, n2);
  arma::vec f2 = arma::vec("0.1 0.4 0.7");
  arma::mat particles_left2 = arma::mat(3, n2);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < n2; j++) {
      PenningTrap obj = PenningTrap(9.25e+1, 2.41e+6, 500., f2(i), w2(j));
      obj.particle_interaction_ = true;
      obj.add_random_particles(100, 1, 40.078);
      obj.simulation(500. / 5000, 5000, "RK4");
      particles_left2(i, j) = obj.particles_in_trap();
    }
  }

  std::ofstream ofile2;
  ofile2.open("Particles_left_Coulomb_on.txt");
  ofile2 << particles_left2 << std::endl
         << f2.as_row() << std::endl
         << w2.as_row();
  ofile2.close();
  return 0;
}
