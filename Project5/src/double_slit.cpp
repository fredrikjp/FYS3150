#include "../include/double_slit.hpp"
using namespace std;

// Constructor
double_slit::double_slit(double h, double dt, double t, double cx, double cy,
                         double dx, double dy, double px, double py, double v0,
                         double wall_width, double wall_postion,
                         double slit_separation, double slit_opening,
                         int n_slits) {
  h_ = h;
  M_ = 1 / h + 1;
  dt_ = dt;
  N_ = t / dt + 1;

  double_slit::generate_potential(v0, wall_width, wall_postion, slit_separation,
                                  slit_opening, n_slits);
  double_slit::generate_matricies();
  double_slit::generate_initial_state(cx, cy, dx, dy, px, py);
}

int double_slit::index_transform(int i, int j) {
  // i is the x-direction index and j is the y-direction index
  int k = i + (M_ - 2) * j;
  return k;
}

void double_slit::generate_matricies() {
  int n = M_ - 2;
  arma::cx_vec a = arma::cx_vec(n * n);
  arma::cx_vec b = arma::cx_vec(n * n);
  complex<double> r(0, dt_ / (2 * h_ * h_));
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      int k = double_slit::index_transform(i, j);
      complex<double> im(0, dt_ * V_(j, i) / 2);
      a(k) = 1. + 4. * r + im;
      b(k) = 1. - 4. * r - im;
    }
  }

  A_ = arma::sp_cx_mat(n * n, n * n);
  B_ = arma::sp_cx_mat(n * n, n * n);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (j < n - 1) {
        A_(i * n + j, i * n + j) = a(i * n + j);
        A_(i * n + j + 1, i * n + j) = -r;
        A_(i * n + j, i * n + j + 1) = -r;
        B_(i * n + j, i * n + j) = b(i * n + j);
        B_(i * n + j + 1, i * n + j) = r;
        B_(i * n + j, i * n + j + 1) = r;
      } else {
        A_(i * n + j, i * n + j) = a(i * n + j);
        B_(i * n + j, i * n + j) = b(i * n + j);
      }
      if (i * n + j < n * n - n) {
        A_(i * n + n + j, i * n + j) = -r;
        A_(i * n + j, i * n + n + j) = -r;
        B_(i * n + n + j, i * n + j) = r;
        B_(i * n + j, i * n + n + j) = r;
      }
    }
  }
}

void double_slit::generate_potential(double v0, double wall_width,
                                     double wall_postion,
                                     double slit_separation,
                                     double slit_opening, int n_slits) {
  int n = M_ - 2;

  V_ = arma::mat(n, n).fill(0);

  int wall_mid_index = round(wall_postion / h_ - 1); // Wall postion x-direction
  int wall_half_width_index_size =
      round((wall_width / h_ - 1) / 2); // Wall width x-direction
  int slit_separation_index_size =
      round(slit_separation / h_); // Seperation wall between slits y-direction
  int slit_opening_index_size =
      round(slit_opening / h_); // Slit opening y-direction
  int wall_center_index =
      round(0.5 / h_ - 1); // Symmetry-axis for the slits in y-direction

  int first_slit_index;
  if (n_slits % 2 == 0) {
    first_slit_index =
        round(wall_center_index -
              (slit_separation_index_size / 2 + slit_opening_index_size +
               (n_slits / 2 - 1) *
                   (slit_separation_index_size + slit_opening_index_size)));
  } else {
    first_slit_index = round(wall_center_index - slit_opening_index_size / 2 -
                             ((n_slits - 1) / 2) * (slit_separation_index_size +
                                                    slit_opening_index_size));
  }
  for (int i = 0; i <= wall_half_width_index_size; i++) {
    // Create wall
    V_.col(wall_mid_index + i) = arma::vec(n).fill(v0);
    V_.col(wall_mid_index - i) = arma::vec(n).fill(v0);
    for (int j = 0; j < n_slits; j++) {
      for (int k = 0; k < slit_opening_index_size; k++) {
        // Create slit openings
        V_(first_slit_index + k +
               j * (slit_separation_index_size + slit_opening_index_size),
           wall_mid_index + i) = 0;
        V_(first_slit_index + k +
               j * (slit_separation_index_size + slit_opening_index_size),
           wall_mid_index - i) = 0;
      }
    }
  }
}

void double_slit::generate_initial_state(double xc, double yc, double dx,
                                         double dy, double px, double py) {

  U_ = arma::cx_mat((M_ - 2) * (M_ - 2), N_);
  for (int i = 0; i < M_ - 2; i++) {
    for (int j = 0; j < M_ - 2; j++) {
      int k = double_slit::index_transform(i, j);
      complex<double> expo(
          -(h_ * (i + 1) - xc) * (h_ * (i + 1) - xc) / (2 * dx * dx) -
              (h_ * (j + 1) - yc) * (h_ * (j + 1) - yc) / (2 * dy * dy),
          px * (h_ * (i + 1) - xc) + py * (h_ * (j + 1) - yc));
      U_(k, 0) = exp(expo);
    }
  }
  U_.col(0) = arma::normalise(U_.col(0));
}

void double_slit::solve_Crank_Nicolson() {
  for (int i = 0; i < N_ - 1; i++) {
    arma::cx_vec u_current = U_.col(i);
    arma::cx_vec u_next = arma::spsolve(A_, B_ * u_current);
    U_.col(i + 1) = u_next;
  }
}

void double_slit::write_to_file(string filename, int prec) {
  // Create and open the output file. Or, technically, create
  // an "output file stream" (type std::ofstream) and connect it to our
  // filename.
  std::ofstream ofile;
  ofile.open(filename);

  // Time
  arma::vec t = arma::linspace(0, (N_ - 1) * dt_, N_);
  // Positions in common for x and y direction
  arma::vec positions = arma::linspace(h_, 1 - h_, M_ - 2);
  
  ofile << t.as_row();
  ofile << positions.as_row();

  // Some width parameters
  int width = 2 * prec + 25;
  for (int i = 0; i < (M_ - 2) * (M_ - 2); i++) {
    for (int j = 0; j < N_; j++) {
      ofile << std::setw(width) << std::setprecision(prec) << std::scientific
            << U_(i, j);
    }
    ofile << endl;
  }
  ofile.close();
}
