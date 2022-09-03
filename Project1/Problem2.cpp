#include <armadillo>
#include <cmath>


arma::vec u(arma::vec x){
	arma::vec u = 1-(1-exp(-10))*x-exp(-10*x);
	return u;
}

int main()
{
	arma::vec x = arma::linspace(0, 1, 1000);
	arma::vec a = u(x);
	a.save("u.bin"); 
	x.save("x.bin"); 

	return 0;
}


