#include "problem2.hpp"

// Create tridiagonal matrix from vectors.
// - lower diagonal: vector a, lenght n-1
// - main diagonal:  vector d, lenght n
// - upper diagonal: vector e, lenght n-1
arma::mat create_tridiagonal(const arma::vec& a, const arma::vec& d, const arma::vec& e)
{
  int n = d.n_elem; 

  // Start from identity matrix
  arma::mat A = arma::mat(n, n, arma::fill::eye);
  
  // Fill first row (row index 0)
  A(0,0) = d(0);
  A(0,1) = e(0);

  // Loop that fills rows 2 to n-1 (row indices 1 to n-2)
  for (int i=1; i <= n-2; i++){
	A(i, i-1) = a(i-1);
	A(i, i) = d(i);
	A(i, i+1) = e(i);
	
  }
  
  // Fill last row (row index n-1)
  A(n-1, n-2) = a(n-2);
  A(n-1, n-1) = d(n-1);

  return A;
}


// Create a tridiagonal matrix tridiag(a,d,e) of size n*n
// from scalar input a, d and e
arma::mat create_tridiagonal(int n, double a, double d, double e)
{
  arma::vec a_vec = arma::vec(n-1, arma::fill::ones) * a;
  arma::vec d_vec = arma::vec(n, arma::fill::ones) * d;
  arma::vec e_vec = arma::vec(n-1, arma::fill::ones) * e;

  // Call the vector version of this function and return the result
  return create_tridiagonal(a_vec, d_vec, e_vec);
}


// Create a symmetric tridiagonal matrix tridiag(a,d,a) of size n*n
// from scalar input a and d.
arma::mat create_symmetric_tridiagonal(int n, double a, double d)
{
  // Call create_tridiagonal and return the result
  return create_tridiagonal(n, a, d, a);
}

void symmetric_tridiagonal_analytic_eigen_solver(arma::vec& eigenval, arma::mat& eigenvec, arma::mat A)
{
  double pi = 2*acos(0.0);	
  double a = A(1,0);
  double d = A(0,0);
  int n = A.n_rows;
  
  eigenval = arma::vec(n);
  eigenvec = arma::mat(n,n);
  

  for (int i = 1; i <= n; i++)
  {
  eigenval(i-1) = d + 2*a*cos(i*pi/(n+1));
  for (int j = 1; j <= n; j++)
  {
  eigenvec(j-1, i-1) = sin(j*i*pi/(n+1));
  } 
  }
  eigenvec = arma::normalise(eigenvec);
}

void test_symmetric_matrix_eigenequation(arma::vec eigenval1, arma::mat eigenvec1, arma::mat X, double prec)
{
arma::vec eigenval2;
arma::mat eigenvec2;
arma::eig_sym(eigenval2, eigenvec2, X);

std::vector<int> eigenval_index;
std::vector<int> eigenvec_index;

for (int i = 0; i <= 5; i++)
{
if (std::abs(eigenval1(i) - eigenval2(i)) <= prec){
}
else{
	eigenval_index.push_back(i);
}
if (std::abs(eigenvec1(0, i) - eigenvec2(0, i)) < prec){
	for (int j = 1; j <= 5; j++)
	{
	if (std::abs(eigenvec1(j, i) - eigenvec2(j, i)) < prec){
	}
	else{
		eigenvec_index.push_back(i);
		break;
	}
	}
}
else if (std::abs(eigenvec1(0, i) + eigenvec2(0, i)) < prec){
	for (int j = 1; j <= 5; j++)
	{
	if (std::abs(eigenvec1(j, i) + eigenvec2(j, i)) < prec){
	}
	else{
		eigenvec_index.push_back(i);
		break;
	}
	}
}
else{
	eigenvec_index.push_back(i);
}
}
if (eigenval_index.size() > 0)
{
std::cout << "Different eigenvalues on index:"; 
for (int i = 0; i < eigenval_index.size(); i++) {
	std::cout << " " << eigenval_index[i];
}
std::cout << std::endl;
}
else
{
std::cout << "The eigenvalues are equal" << std::endl; 
}
if (eigenvec_index.size() > 0)
{
std::cout << "Different eigenvectors on column index:"; 
for (int i = 0; i < eigenvec_index.size(); i++) {
	std::cout << " " << eigenvec_index[i];
}
std::cout << std::endl;
}
else
{
std::cout << "The eigenvectors are equal" << std::endl; 
}
std::cout << "Test ran with equality to precission "<< prec << std::endl;
}








