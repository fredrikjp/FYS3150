#include <armadillo>
#include <cmath>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <unistd.h>

//a, b and c are respectively sub, main and super diagonal of a tridiagonal matrix, while g is the rhs of the matrix eq.
arma::vec Special(arma::vec a, arma::vec b, arma::vec c, arma::vec g) 
{
	int n = b.n_elem;
	arma::vec _b = arma::vec(n);
	arma::vec _g = arma::vec(n);
	arma::vec v = arma::vec(n);
	_b(0) = b(0);	
	_g(0) = g(0);

	for (int i = 1; i <= n-1; i++)
	{
		_b(i) = 2 - 1./_b(i-1);	//Element (i-1) in a corresponds to row i in the matrix
		_g(i) = g(i) + _g(i-1)/_b(i-1);
	}
	
	v(n-1) = _g(n-1)/_b(n-1);
	for (int i = n-2; i >= 0; i--)
	{
		v(i) = (_g(i) + v(i+1))/_b(i);
	}

	return v;
}

 











