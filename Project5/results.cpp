#include "include/double_slit.hpp"

int main() {
  double h = 0.005;
  double dt = 0.000025;
  double t = 0.008;
  double cx = 0.25;
  double cy = 0.5;
  double dx = 0.05;
  double dy = 0.05;
  double px = 200;
  double py = 0;
  double v0 = 0;
  double wall_width = 0.02;
  double wall_postion = 0.5;
  double slit_separation = 0.05;
  double slit_opening = 0.05;
  int n_slits = 2;

  // Problem 7 
  ///*
  double_slit obj1 =
      double_slit(h, dt, t, cx, cy, dx, dy, px, py, v0, wall_width,
                  wall_postion, slit_separation, slit_opening, n_slits);
  obj1.solve_Crank_Nicolson();
  obj1.write_to_file("problem7a.txt", 15);

  v0 = 1e+10;
  dy = 0.10;

  double_slit obj2 =
      double_slit(h, dt, t, cx, cy, dx, dy, px, py, v0, wall_width,
                  wall_postion, slit_separation, slit_opening, n_slits);
  obj2.solve_Crank_Nicolson();
  obj2.write_to_file("problem7b.txt", 15);
  //*/

  // Problem 8
  ///*
  t = 0.002;
  dy = 0.20; 
  v0 = 1e+10;

  double_slit obj3 =
      double_slit(h, dt, t, cx, cy, dx, dy, px, py, v0, wall_width,
                  wall_postion, slit_separation, slit_opening, n_slits);
  obj3.solve_Crank_Nicolson();
  obj3.write_to_file("problem8.txt", 15);
  //*/

  // Problem 9
  ///*
  t = 0.002;
  dy = 0.20; 
  v0 = 1e+10;

  n_slits = 1;
  double_slit obj4 = 
      double_slit(h, dt, t, cx, cy, dx, dy, px, py, v0, wall_width,
                  wall_postion, slit_separation, slit_opening, n_slits);
  obj4.solve_Crank_Nicolson();
  obj4.write_to_file("problem9_single.txt", 15);

  n_slits = 3;
  double_slit obj5 = 
      double_slit(h, dt, t, cx, cy, dx, dy, px, py, v0, wall_width,
                  wall_postion, slit_separation, slit_opening, n_slits);
  obj5.solve_Crank_Nicolson();
  obj5.write_to_file("problem9_triple.txt", 15);
  //*/

  return 0;
}
