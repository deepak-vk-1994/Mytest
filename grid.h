#ifndef _GRID_H
#define _GRID_H
#include<math.h>

struct Point {
	double x,y,z;
};

class Element {
	private:
			
	public:
		int vertex[3]; //Vertex index
		int face[3];
		double p,u,v,w,rho;
		double volume,xc,yc;
		double gradvel[4];
		double gradp[2];
		double gradrho[2];
		double gradtemp[2];
		double phi[4]; // Limiter

		Element() {
			p = p_inf;
			rho = rho_inf;
			u = u_inf;
			v = v_inf;
			w = w_inf;

			volume = 0.0;
			xc = yc = 0.0;
			for (int i = 0; i < 3; i++) {
				vertex[i] = -1;
				face[i] = -1;
			}
			for (int i = 0; i < 4; i++) {
				gradvel[i] = 0.0;
				phi[i] = 0.0;
			}
			gradtemp[0] = gradtemp[1] = 0.0;
			gradp[0] = gradp[1] = 0.0;
			gradrho[0] = gradrho[1] = 0.0;
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
//			return mu_ref;
		}
};

struct Face {
	int vertex1, vertex2;
	int elementL, elementR;
	int marker;
	double mass_flux, mom_x_flux, mom_y_flux, energy_flux;
	double nx,ny,nz,area,xmp,ymp;
	double rho_face,u_face,v_face,p_face,T_face;
	Face() {
		vertex1 = vertex2 = elementL = elementR = -1;
		mass_flux = mom_x_flux = mom_y_flux = energy_flux = 0.0;
		nx = ny = nz = area = 0.0;
		rho_face = u_face = v_face = p_face = T_face = 0.0;
		xmp = ymp = 0.0;
		marker = INTERIOR;
	}
};


#endif

