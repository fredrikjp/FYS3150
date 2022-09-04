#include <armadillo>
#include <cmath>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>

arma::vec u(arma::vec x){
	arma::vec u = 1-(1-exp(-10))*x-exp(-10*x);
	return u;
}

int main()
{
	int n = 1000;
	arma::vec x = arma::linspace(0, 1, n);
	arma::vec a = u(x);

	std::string filename = "x_u(x).txt";

	std::ofstream ofile;
	ofile.open(filename);

	int width = 10;
	int prec = 2;
	for (int i = 0; i < n; i++)
	{
		ofile << std::setw(width) << std::setprecision(prec) << std::scientific << x(i)
			<< std::setw(width) << std::setprecision(prec) << std::scientific << a(i)
			<< std::endl;

	}
	
	ofile.close();

	return 0;
}


