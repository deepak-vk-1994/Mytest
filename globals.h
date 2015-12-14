#ifndef _GLOBALS_H
#define _GLOBALS_H

//#define DEBUG

//Flux splitting scheme
#define LOUISTEFFAN
#define FIRSTORDER
//#define SECONDORDER
//#define VISCOUS


enum type{INTERIOR,LEFT,RIGHT,TOP,BOTTOM,OTHER,GHOST};
enum mach{SUBSONIC,SUPERSONIC};


int SIMULATION_TIME;
int mach_switch = 0; //Default Subsonic

double nx_facet = 0;
double ny_facet = 0;
double nz_facet = 1;

int N_elem;
int N_pts;
/*
double D=1;
double xleft = -3.5;
double xright = 3.5;
double ybottom = 0.0;
double ytop = 2.0;
*/

//Ramp
double D=1;
double xleft = 0.0;
double xright = 1.0;
double ytop = 0.18;
double ybottom = 0.0;

/*
//CYLINDER
double D = 1e-1;
double xleft = -1;
double xright = 3;
double ybottom = -1;
double ytop = 1;
*/
//Poiseuelle
// double D = 1e-5;
//double xleft = -0.0003;
//double xright = 0.0003;
//double ybottom = -0.0001;
//double ytop = 0.0001;

double p_inf;
double rho_inf;
double u_inf = 0.0; 
double v_inf = 0.0; 
double w_inf = 0.0;
double CFL;
double Re = 0.0;

double gam = 1.4;
double R = 287.0;

//Viscous references
double mu_ref;
double T_ref;
double S1 = 110.4; // Sutherland's constant
double Pr;  //Prandlt Number


/*  --DEBUG grid-- */

//int N_tri = 4;
//int N_pts = 5;
//double xleft = 0;
//double xright = 3.4;
//double ybottom = 0;
//double ytop = 2.2;


double max(double a, double b) {
	return (a<b)?b:a;
}
double min(double a, double b) {
	return (a<b)?a:b;
}

#endif
