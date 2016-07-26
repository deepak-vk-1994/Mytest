#ifndef _GLOBALS_H
#define _GLOBALS_H

#include <fstream>
#include <string>
#include <math.h>
#include <vector>

enum type{INTERIOR,LEFT,RIGHT,TOP,BOTTOM,OTHER};
enum mach{SUBSONIC,SUPERSONIC};
enum BC{INLET,OUTLET,SLIP,NO_SLIP};
enum inviscidflux{AUSM};
enum accuracyorder{FO,HO,RK4};
enum viscous{NO_VISC};
enum turbu{NO_TURB,SST};
enum other{PI};

int debug = 0;
int SIMULATION_TIME;
int in_flux,ord_accuracy,visc_flux,turbulence,time_int,other;
std::string gridfile;
std::string restartfile;
int mach_switch = 0; //Default Subsonic
int global_time = 0;
int N_elem = 0;
int N_pts = 0;
double n_facet = 1.0;
double xleft,xright,ybottom,ytop,D,chord;
double p_inf;
double rho_inf;
double p0,rho0,M;
double u_inf = 0.0; 
double v_inf = 0.0; 
double w_inf = 0.0;
double k_inf = 0.0;
double omega_inf = 0.0;
double CFL;
double Re = 0.0;
double visc_inf = 0.0;
double nu_inf = 0.0;
double T_inf = 0.0;
double gam = 1.4;
double R = 287.0;
double mu_ref;
double T_ref;
double S1 = 110.4; 
double Pr = 0.72;

int BCL,BCR,BCT,BCB,BCO;
int ghostoutput = 0;
int restart = 0;
int ic = 0;
std::vector<int> bcs, bcns;
std::vector<int>::iterator bcs_it, bcns_it;

void setInputVariables() {
	std::ifstream file("./Input.txt");
	std::string line;
	getline (file,line);
	getline (file,line);
	if (stod(line) == 1) {
		debug = 1;
	}
	getline (file,line);
	getline (file,line);
	gridfile = line;
	getline (file,line);
	getline (file,line);
	n_facet = stoi(line);
	getline (file,line);
	getline (file,line);
	xleft = stod(line);
	getline (file,line);
	getline (file,line);
	xright = stod(line);
	getline (file,line);
	getline (file,line);
	ybottom = stod(line);
	getline (file,line);
	getline (file,line);
	ytop = stod(line);
	getline (file,line);
	getline (file,line);
	D = stod(line);
	getline (file,line);
	getline (file,line);
	p0 = stod(line);
	getline (file,line);
	getline (file,line);
	rho0 = stod(line);
	getline (file,line);
	getline (file,line);
	M = stod(line);
	getline (file,line);
	getline (file,line);
	T_ref = stod(line);
	getline (file,line);
	getline (file,line);
	mu_ref = stod(line);
	getline (file,line);
	getline (file,line);
	SIMULATION_TIME = stoi(line);
	getline (file,line);
	getline (file,line);
	mach_switch = stoi(line);
	getline (file,line);
	getline (file,line);
	if (line == "AUSM") 
		in_flux = AUSM;
	getline (file,line);
	getline (file,line);
	if (line == "FO") {
		ord_accuracy = FO;
	}
	else {
		ord_accuracy = HO;
	}
	getline (file,line);
	getline (file,line);
	if (line == "HO")
		visc_flux = HO;
	else
		visc_flux = NO_VISC;
	getline (file,line);
	getline (file,line);
	if (line == "RK4")
		time_int = RK4;
	getline (file,line);
	getline (file,line);
	CFL = stod(line);
	getline (file,line);
	getline (file,line);
	if (line == "SST")
		turbulence = SST;
	else
		turbulence = NO_TURB;
	getline (file,line);
	getline (file,line);
	if (line == "PI")
		other = PI;
	getline (file,line);
	getline (file,line);
	BCL = stoi(line);
	getline (file,line);
	getline (file,line);
	BCR = stoi(line);
	getline (file,line);
	getline (file,line);
	BCT = stoi(line);
	getline (file,line);
	getline (file,line);
	BCB = stoi(line);
	getline (file,line);
	getline (file,line);
	BCO = stoi(line);
	getline (file,line);
	getline (file,line);
	ghostoutput = stoi(line);
	getline (file,line);
	getline (file,line);
	restart = stoi(line);
	if (restart == 1) {
		getline (file,line);
		getline (file,line);
		getline (file,line);
		restartfile = line;
		getline (file,line);
		getline (file,line);
		N_elem = stoi(line);
		getline (file,line);
		getline (file,line);
		N_pts = stoi(line);
		getline (file,line);
		getline (file,line);
		ic = stoi(line);
	}
	chord = D;
	p_inf = p0/(pow((1 + (gam-1)*0.5*M*M),(gam/(gam-1))));
	rho_inf = rho0/(pow((1 + (gam-1)*0.5*M*M),(1/(gam-1))));
	u_inf = M*sqrt(gam*p_inf/rho_inf);
	T_inf = p_inf/(rho_inf*R);
	visc_inf = ( mu_ref * pow((T_inf/T_ref),1.5) * (T_ref + S1)/(T_inf + S1) );
	nu_inf = visc_inf/rho_inf;
	k_inf = 1.5*(0.05*u_inf)*(0.05*u_inf);
	omega_inf = k_inf/(1e-3*nu_inf);
	Re = rho_inf*u_inf*D/visc_inf;
	std::cout << "P: " << p_inf <<" RHO: " << rho_inf <<" U: " << u_inf <<" Re: "<< Re <<" K: "<< k_inf <<" OMEGA: " << omega_inf << std::endl;
}


double turbphi_k1 = 0.85;
double turbphi_k2 = 1.0;
double turbphi_omega1 = 0.5;
double turbphi_omega2 = 0.856;
double betastar = 0.09;
double beta1 = 0.075;
double beta2 = 0.0828;
double kappa = 0.41;
double a1 = 0.31;
double Pr_t = 0.9;


double max(double a, double b) {
	return (a<b)?b:a;
}
double min(double a, double b) {
	return (a<b)?a:b;
}

#endif

