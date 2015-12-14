#include<iostream>
#include<math.h>
using namespace std;

int main() {
	double rho,u,mu,tau,delta,Re,x;

	rho = 1.20084;
	u = 133;
	x = 1;
	mu = 1.77637e-5;

	Re = rho*u*x/mu;
	delta = 4.91 * x/sqrt(Re);
	cout << "Laminar flat plate BL thickness: "<<delta<<endl;

	tau = mu * u/delta;
	u = sqrt(tau/mu);

	double y_min;
	y_min = mu/(rho*u);
	cout << "Resolution needed is: "<<y_min<<endl;

}
