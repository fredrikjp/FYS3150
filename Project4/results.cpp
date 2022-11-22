#include "include/MCMC.hpp"
#include <string>

int main() {

  // 2X2 ising model
   /*
  MCMC obj1 = MCMC(2, 1, false);
  obj1.Monte_Carlo(3000, 0);
  obj1.write_to_file("problem4.txt");
   */

  // 20X20 ising model ordered vs unordered
  /*
 // Ordered starting point
 MCMC obj2 = MCMC(20, 1, true);
 obj2.Monte_Carlo(5000, 0);
 obj2.write_to_file("T1_ordered.txt");

 // Unordered starting point
 MCMC obj3 = MCMC(20, 1, false);
 obj3.Monte_Carlo(5000, 0);
 obj3.write_to_file("T1_unordered.txt");

 // Ordered starting point
 MCMC obj4 = MCMC(20, 2.4, true);
 obj4.Monte_Carlo(5000, 0);
 obj4.write_to_file("T2_ordered.txt");

 // Unordered starting point
 MCMC obj5 = MCMC(20, 2.4, false);
 obj5.Monte_Carlo(5000, 0);
 obj5.write_to_file("T2_unordered.txt");
  */

  // 20X20 ising model for generating probability distribution
  /*
 // Ordered starting point
 MCMC obj6 = MCMC(20, 1, true);
 obj6.Monte_Carlo(3000, 250);
 obj6.write_m_e_to_file("T1_e.txt");


 // Ordered starting point
 MCMC obj7 = MCMC(20, 2.4, true);
 obj7.Monte_Carlo(3000, 250);
 obj7.write_m_e_to_file("T2_e.txt");
  */

  // Parallelelization dividing temperature models to different threads
  // time measurements
  /*
    int nT = 8;
    arma::vec T_ = arma::linspace(2.1, 2.4, nT);

for (int i = 1; i<4; i++)
{
    int l = 1000*i;

    auto t1 = std::chrono::high_resolution_clock::now();
  #pragma omp parallel // Start parallel region
    {
  // Here we start the parallelized loop over T
  #pragma omp for
      for (int i = 0; i < nT; i++) {
        MCMC thread_obj = MCMC(20, T_(i), false);
        thread_obj.Monte_Carlo(l, 0);
      } // End parallelized loop over T

    } // End entire parallel region
    auto t2 = std::chrono::high_resolution_clock::now();
    double parallel_duration_seconds =
        std::chrono::duration<double>(t2 - t1).count();

    // Same without parallelization
    auto t3 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < nT; i++) {
      MCMC thread_obj = MCMC(20, T_(i), false);
      thread_obj.Monte_Carlo(l, 0);
    }

    auto t4 = std::chrono::high_resolution_clock::now();
    double plain_duration_seconds =
        std::chrono::duration<double>(t4 - t3).count();

    cout << "Parallelelization time[s]:" << parallel_duration_seconds
         << " | Plain time[s]" << plain_duration_seconds << endl;
}
   */

// Quantities for different L and T
 /*   
  int n_T = 20;
  arma::vec T = arma::linspace(2.1, 2.4, n_T);

  int n = 1000;
  arma::vec L_vec = arma::vec("40 60 80 100");
  int m = L_vec.n_elem;
  arma::mat quantities_T = arma::mat(n_T, 4);
  for (int j = 0; j < m; j++) {
    int L = L_vec(j);
#pragma omp parallel // Start parallel region
    {
// Here we start the parallelized loop over T
#pragma omp for
      for (int i = 0; i < n_T; i++) {
        MCMC thread_obj = MCMC(L, T(i), true);
        thread_obj.Monte_Carlo(n, 250);
        quantities_T(i, 0) = thread_obj.mean_e_.back();
        quantities_T(i, 1) = thread_obj.mean_m_.back();
        quantities_T(i, 2) = thread_obj.CV_.back();
        quantities_T(i, 3) = thread_obj.chi_.back();
      } // End parallelized loop over T

    } // End entire parallel region

    string filename = "L" + to_string(L) + "_model.txt";

    std::ofstream ofile;
    ofile.open(filename);

    ofile << quantities_T;
    ofile.close();
  }
  */ 

  return 1;
}
