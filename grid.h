#ifndef _GRID_H
#define _GRID_H
#include <math.h>
#include <float.h>

struct Point {
	double x,y,z;
};

class Element {
	private:
			
	public:
		//Geometry
		int vertex[3]; //Vertex index
		int face[3];
		double volume,xc,yc,dist;
		double delta_t;

		//Primitive variables
		double p,u,v,w,rho;  //Laminar variables
		double k,omega; //Turbulent variables
		double p_temp,u_temp,v_temp,w_temp,rho_temp;
		double k_temp,omega_temp;
		double F1,mu_t;
		//Other variables
		double source1,source2,beta; //magnitude of vorticity

		//Gradients
		double gradu[2];
		double gradv[2];
		double gradp[2];
		double gradrho[2];
		double gradtemp[2];
		double gradentropy[2];

		//turbulent gradients
		double gradk[2],gradomega[2];

		double phi[6]; // Limiter
		double E[4][6],residue[6]; 

		double check;
		double pk,dk,pw,dw,lam;
		//Default constructor
		Element() {
			p = p_inf;
			rho = rho_inf;
			u = u_inf;
			v = v_inf;
			w = w_inf;
			k = k_inf;
			omega = omega_inf;

			volume = 0.0;
			xc = yc = 0.0;
			dist = 0.0;

			delta_t = 0.0;
			for (int i = 0; i < 3; i++) {
				vertex[i] = -1;
				face[i] = -1;
			}
			for (int i = 0; i < 4; i++) {
				phi[i] = 0.0;
				
				for (int j = 0; j < 4; j++) E[i][j] = 0.0;
			}
			for (int i = 0; i < 6; i++) {
				residue[i] = 0.0;
			}
			gradu[0] = gradu[1] = 0.0;
			gradv[0] = gradv[1] = 0.0;
			gradtemp[0] = gradtemp[1] = 0.0;
			gradp[0] = gradp[1] = 0.0;
			gradrho[0] = gradrho[1] = 0.0;
			gradk[0] = gradk[1] = 0.0;
			gradomega[0] = gradomega[1] = 0.0;
			gradentropy[0] = gradentropy[1] = 0.0;

			F1 = mu_t = 0.0;
			source1 = source2 = 0.0;
			beta = 0.0;

			check = 0.0;
			
		}

		//Member Fns
		double calEnthalpy() {
			return ( gam/(gam-1)*(p/rho) + 0.5*(u*u + v*v + w*w) );
		}
		double calTemp() {
			return ( p/(rho*R) );
		}
		double calViscosity() {
			double T = this->calTemp();
			return ( mu_ref * pow((T/T_ref),1.5) * (T_ref + S1)/(T + S1) );  //Sutherland's Law
		}

		double calVorticity() {
			return fabs((gradv[0] - gradu[1]));
		}

		double calStrainrate() {
			return sqrt(  2*(gradu[0]*gradu[0]) + 2*(gradv[1]*gradv[1]) + (gradu[1] + gradv[0])*(gradu[1] + gradv[0]) );
		}

		double calCrossDiffusion() {
			return (2.0*rho*turbphi_omega2*(1.0/omega)*(gradk[0]*gradomega[0] + gradk[1]*gradomega[1]));
		}

		double calEntropy() {
			return ( p/pow(rho,gam));
		}
		
		void calTurbulentVisc() {
			double arg2,F2,nu,vort;
			nu = this->calViscosity()/rho;
			vort = this->calVorticity();
			arg2 = std::max(2.0*sqrt(k)/(betastar*omega*dist),500.0*nu/(dist*dist*omega));
			F2 = tanh(arg2*arg2);
			mu_t = rho*a1*k/(std::max(a1*omega,vort*F2));
		}

		void calF1() {
			double arg1,nu,CDkw;
			nu = this->calViscosity()/rho;
			CDkw = std::max(this->calCrossDiffusion(),1e-20);
			arg1 = std::min(std::max(sqrt(k)/(betastar*omega*dist),500.0*nu/(dist*dist*omega)) , 4.0*rho*turbphi_omega2*k/(CDkw*dist*dist));
			F1 = tanh(arg1*arg1*arg1*arg1);
		}
};

struct Face {
	int vertex1, vertex2;
	int elementL, elementR;
	int marker;
	//Flux
	double k_flux, omega_flux;
	double mass_flux, mom_x_flux, mom_y_flux, energy_flux;
	//Turbulent fluxes

	double nx,ny,nz,area,xmp,ymp;
	double rho_face,u_face,v_face,p_face,T_face,en_face;  //For calculating gradients
	//turbulence gradients
	double k_face,omega_face;
	Face() {
		vertex1 = vertex2 = elementL = elementR = -1;
		mass_flux = mom_x_flux = mom_y_flux = energy_flux = 0.0;
		k_flux = omega_flux = 0.0;
		nx = ny = nz = area = 0.0;
		rho_face = u_face = v_face = p_face = T_face = k_face = omega_face = en_face = 0.0;
		xmp = ymp = 0.0;
		marker = INTERIOR;
	}
};


#endif

