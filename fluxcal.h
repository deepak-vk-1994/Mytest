#ifndef _FLUXCAL_H
#define _FLUXCAL_H

#include<math.h>
#include<assert.h>
#include "globals.h"
#include "grid.h"

double maximum(double a, double b, double c, double d) {
	return std::max(std::max(a,b),std::max(c,d));
}

double minimum(double a, double b, double c, double d) {
	return std::min(std::min(a,b),std::min(c,d));
}

void calLimiters(int index, std::vector<Element > &elements, std::vector<Face > &faces, double phi[4]) {
	double q[4], qmax[4], qmin[4], dot[4];
	double xmc, ymc, delmax, delmin,lim;
	// std::vector<double> q(4),qmax(4), qmin(4), dot(4);
	int n = 0;
	double epsilon,K,h; //Venkatakrishnan parameters
	h = cbrt(elements[index].volume);
	K = 5; //Change according to solution
	epsilon = (K*h)*(K*h)*(K*h);

	for (int i = 0; i < 4; i++) phi[i] = 0.0;
	if (index >= N_elem) {//Do nothing
	}
	else {
		std::vector<int> nbour(3);
		for (int i = 0; i < 3; i++) {
			n = faces[elements[index].face[i]].elementL;
			if (n == index)
				n = faces[elements[index].face[i]].elementR;
			nbour[i] = n;
			assert(nbour[i] != index);
		}

		q[0] = elements[index].rho; q[1] = elements[index].u; q[2] = elements[index].v; q[3] = elements[index].p;
		qmax[0] = maximum(elements[nbour[0]].rho,elements[nbour[1]].rho,elements[nbour[2]].rho,elements[index].rho);
		qmax[1] = maximum(elements[nbour[0]].u,elements[nbour[1]].u,elements[nbour[2]].u,elements[index].u);
		qmax[2] = maximum(elements[nbour[0]].v,elements[nbour[1]].v,elements[nbour[2]].v,elements[index].v);
		qmax[3] = maximum(elements[nbour[0]].p,elements[nbour[1]].p,elements[nbour[2]].p,elements[index].p);

		qmin[0] = minimum(elements[nbour[0]].rho,elements[nbour[1]].rho,elements[nbour[2]].rho,elements[index].rho);
		qmin[1] = minimum(elements[nbour[0]].u,elements[nbour[1]].u,elements[nbour[2]].u,elements[index].u);
		qmin[2] = minimum(elements[nbour[0]].v,elements[nbour[1]].v,elements[nbour[2]].v,elements[index].v);
		qmin[3] = minimum(elements[nbour[0]].p,elements[nbour[1]].p,elements[nbour[2]].p,elements[index].p);

		for (int i = 0; i < 3; i++) {
			xmc = faces[elements[index].face[i]].xmp - elements[index].xc;
			ymc = faces[elements[index].face[i]].ymp - elements[index].yc;
			dot[0] = elements[index].gradrho[0]*xmc + elements[index].gradrho[1]*ymc;
			dot[1] = elements[index].gradvel[0]*xmc + elements[index].gradvel[1]*ymc;
			dot[2] = elements[index].gradvel[2]*xmc + elements[index].gradvel[3]*ymc;
			dot[3] = elements[index].gradp[0]*xmc + elements[index].gradp[1]*ymc;

			for (int j = 0; j < 4; j++) {
				delmax = qmax[j] - q[j];
				delmin = qmin[j] - q[j];


				if (dot[j] > 0)
				 	lim = (1/dot[j]) * ( ((delmax*delmax + epsilon)*dot[j] + 2*dot[j]*dot[j]*delmax)/(delmax*delmax + 2*dot[j]*dot[j] + delmax*dot[j] + epsilon) );

				else if (dot[j] < 0)
					lim = (1/dot[j]) * ( ((delmin*delmin + epsilon)*dot[j] + 2*dot[j]*dot[j]*delmin)/(delmin*delmin + 2*dot[j]*dot[j] + delmin*dot[j] + epsilon) );

				else
					lim = 1.0;

				if (i == 0) phi[j] = lim;
				else phi[j] = std::min(phi[j],lim);

			}

		}
	}
	
}

void calFlux(std::vector<Point > &points, std::vector<Element > &elements, std::vector<Face > &faces) {
	int EL, ER;
	int neighbour;

	//Liou steffan variables
	double ML, MR, a_face;
	double alpha_plus, alpha_minus, beta_plus, beta_minus, C_VL_plus, C_VL_minus, C_LS_plus, C_LS_minus;
	double D_plus, D_minus;
	double E_t_plus, E_t_minus;

	//Viscous flux variables
	double txx,tyy,txy,mu,K,qx,qy,dudx,dudy,dvdx,dvdy,dTdx,dTdy,dudx_f,dudy_f,dvdx_f,dvdy_f,dTdx_f,dTdy_f;

	//Others
	double u_f,v_f;
	double uL,uR,vL,vR,rhoL,rhoR,pL,pR,TL,TR;
	double xL,xR,yL,yR,r,xmcL,ymcL,xmcR,ymcR;
	double phiL[4],phiR[4]; //Limiters for each primitive variable
	int num = 0;
	

	for (int i = 0; i < N_elem; i++) {


		for (int k = 0; k < 4; k++) {
			elements[i].gradvel[k] = 0.0;
 		}


		elements[i].gradtemp[0] = elements[i].gradtemp[1] = 0.0;
		elements[i].gradrho[0] = elements[i].gradrho[1] = 0.0;
		elements[i].gradp[0] = elements[i].gradp[1] = 0.0;



		for (int j = 0; j < 3; j++) {
			num = elements[i].face[j];
			EL = faces[num].elementL;
			ER = faces[num].elementR;
			
			faces[num].rho_face = 0.5 * (elements[EL].rho + elements[ER].rho);
			faces[num].u_face = 0.5 * (elements[EL].u + elements[ER].u);
			faces[num].v_face = 0.5 * (elements[EL].v + elements[ER].v);
			faces[num].p_face= 0.5 * (elements[EL].p + elements[ER].p);
			faces[num].T_face = 0.5 * (elements[EL].calTemp() + elements[ER].calTemp());


			if (i == ER) {
			elements[i].gradvel[0] += faces[num].area/elements[i].volume * faces[num].u_face*faces[num].nx;
			elements[i].gradvel[1] += faces[num].area/elements[i].volume * faces[num].u_face*faces[num].ny;
			elements[i].gradvel[2] += faces[num].area/elements[i].volume * faces[num].v_face*faces[num].nx;
			elements[i].gradvel[3] += faces[num].area/elements[i].volume * faces[num].v_face*faces[num].ny;
        	elements[i].gradtemp[0] += faces[num].area/elements[i].volume * faces[num].T_face*faces[num].nx;
			elements[i].gradtemp[1] += faces[num].area/elements[i].volume * faces[num].T_face*faces[num].ny;
			elements[i].gradrho[0] += faces[num].area/elements[i].volume * faces[num].rho_face*faces[num].nx;
			elements[i].gradrho[1] += faces[num].area/elements[i].volume * faces[num].rho_face*faces[num].ny;
			elements[i].gradp[0] += faces[num].area/elements[i].volume * faces[num].p_face*faces[num].nx;
			elements[i].gradp[1] += faces[num].area/elements[i].volume * faces[num].p_face*faces[num].ny;

            }

            else if (i == EL) {
			elements[i].gradvel[0] += -1.0*faces[num].area/elements[i].volume * faces[num].u_face*faces[num].nx;
			elements[i].gradvel[1] += -1.0*faces[num].area/elements[i].volume * faces[num].u_face*faces[num].ny;
			elements[i].gradvel[2] += -1.0*faces[num].area/elements[i].volume * faces[num].v_face*faces[num].nx;
			elements[i].gradvel[3] += -1.0*faces[num].area/elements[i].volume * faces[num].v_face*faces[num].ny;
        	elements[i].gradtemp[0] += -1.0*faces[num].area/elements[i].volume * faces[num].T_face*faces[num].nx;
			elements[i].gradtemp[1] += -1.0*faces[num].area/elements[i].volume * faces[num].T_face*faces[num].ny;
			elements[i].gradrho[0] += -1.0*faces[num].area/elements[i].volume * faces[num].rho_face*faces[num].nx;
			elements[i].gradrho[1] += -1.0*faces[num].area/elements[i].volume * faces[num].rho_face*faces[num].ny;
			elements[i].gradp[0] += -1.0*faces[num].area/elements[i].volume * faces[num].p_face*faces[num].nx;
			elements[i].gradp[1] += -1.0*faces[num].area/elements[i].volume * faces[num].p_face*faces[num].ny;
            }
			
		}

		calLimiters(i,elements,faces,elements[i].phi);
		
	}


/*------------------------------------------------------------------------------------------------------------------------------------*/

	//FLUX CALCULATION
	for (std::vector<Face >::iterator it = faces.begin(); it != faces.end(); it++) {
		EL = (*it).elementL;
		ER = (*it).elementR;		
		
		xmcL = (*it).xmp - elements[EL].xc; ymcL = (*it).ymp - elements[EL].yc;
		xmcR = (*it).xmp - elements[ER].xc; ymcR = (*it).ymp - elements[ER].yc;

		for (int i = 0; i < 4; i++){
			phiL[i] = elements[EL].phi[i];
			phiR[i] = elements[ER].phi[i];
		}
		//Reconstruction of primitive variables
		//First Order
		#ifdef FIRSTORDER
			rhoL = elements[EL].rho;
			rhoR = elements[ER].rho;
			uL = elements[EL].u;
			uR = elements[ER].u;
			vL = elements[EL].v;
			vR = elements[ER].v;
			pL = elements[EL].p;
			pR = elements[ER].p;
		#endif
		//Second Order
		#ifdef SECONDORDER
			rhoL = elements[EL].rho + phiL[0]*(elements[EL].gradrho[0]*xmcL + elements[EL].gradrho[1]*ymcL);
			rhoR = elements[ER].rho + phiR[0]*(elements[ER].gradrho[0]*xmcR + elements[ER].gradrho[1]*ymcR);
			uL = elements[EL].u + phiL[1]*(elements[EL].gradvel[0]*xmcL + elements[EL].gradvel[1]*ymcL);
			uR = elements[ER].u + phiR[1]*(elements[ER].gradvel[0]*xmcR + elements[ER].gradvel[1]*ymcR);
			vL = elements[EL].v + phiL[2]*(elements[EL].gradvel[2]*xmcL + elements[EL].gradvel[3]*ymcL);
			vR = elements[ER].v + phiR[2]*(elements[ER].gradvel[2]*xmcR + elements[ER].gradvel[3]*ymcR);
			pL = elements[EL].p + phiL[3]*(elements[EL].gradp[0]*xmcL + elements[EL].gradp[1]*ymcL);
			pR = elements[ER].p + phiR[3]*(elements[ER].gradp[0]*xmcR + elements[ER].gradp[1]*ymcR);
		#endif
		TL = pL/(rhoL*R);
		TR = pR/(rhoR*R);
		E_t_plus = gam/(gam-1)*(pR/rhoR) + 0.5*(uR*uR + vR*vR);
		E_t_minus = gam/(gam-1)*(pL/rhoL) + 0.5*(uL*uL + vL*vL);	
		
		a_face = 0.5 * ( sqrt(gam*(pL/rhoL)) + sqrt(gam*(pR/rhoR)) );
		ML = (uL*((*it).nx) + vL*((*it).ny))/a_face;
		MR = (uR*((*it).nx) + vR*((*it).ny))/a_face;


		//INVISCID FLUXES
		#ifdef LOUISTEFFAN
			alpha_plus = 0.5 * (1.0 + ((MR > 0) - (MR < 0)) );
			alpha_minus = 0.5 * (1.0 - ((ML > 0) - (ML < 0)) );

			beta_plus = -1.0 * max(0.0, 1.0-floor(fabs(MR)));
			beta_minus = -1.0 * max(0.0, 1.0-floor(fabs(ML)));
		
			C_VL_plus = alpha_plus*(1 + beta_plus)*MR - 0.25*beta_plus*(MR + 1)*(MR + 1);
			C_VL_minus = alpha_minus*(1 + beta_minus)*ML + 0.25*beta_minus*(ML - 1)*(ML - 1);

			C_LS_plus = max((C_VL_plus + C_VL_minus),0);
			C_LS_minus = min((C_VL_plus + C_VL_minus),0);

			D_plus = alpha_plus*(1 + beta_plus) - beta_plus*(0.25*(1 + MR)*(1 + MR)*(2 - MR));
			D_minus = alpha_minus*(1 + beta_minus) - beta_minus*(0.25*(1 - ML)*(1 - ML)*(2 + ML));		

			
	
			(*it).mass_flux = ((*it).area) * (rhoR*a_face*C_LS_plus*1.0 + rhoL*a_face*C_LS_minus*1.0 + D_plus*pR*0.0 + D_minus*pL*0.0 );
			(*it).mom_x_flux = ((*it).area) * (rhoR*a_face*C_LS_plus*uR + rhoL*a_face*C_LS_minus*uL + D_plus*pR*((*it).nx) + D_minus*pL*((*it).nx) );
			(*it).mom_y_flux = ((*it).area) * (rhoR*a_face*C_LS_plus*vR + rhoL*a_face*C_LS_minus*vL + D_plus*pR*((*it).ny) + D_minus*pL*((*it).ny) );
			(*it).energy_flux = ((*it).area) * (rhoR*a_face*C_LS_plus*E_t_plus + rhoL*a_face*C_LS_minus*E_t_minus + D_plus*pR*0.0 + D_minus*pL*0.0 );

			if ((*it).marker == TOP || (*it).marker == BOTTOM) {
			 	(*it).mass_flux = 0.0;
			 	(*it).mom_x_flux = ((*it).area) * (D_plus*pR*((*it).nx) + D_minus*pL*((*it).nx) );
			 	(*it).mom_y_flux = ((*it).area) * (D_plus*pR*((*it).ny) + D_minus*pL*((*it).ny) );
			 	(*it).energy_flux = 0.0;
			 }

		#endif	



		//VISCOUS FLUXES
		#ifdef VISCOUS
			xL = elements[EL].xc;
			yL = elements[EL].yc;		
			xR = elements[ER].xc;
			yR = elements[ER].yc;
			r = sqrt((xR-xL)*(xR-xL) + (yR-yL)*(yR-yL));  //Vector magnitude from centroidL to centroidR

			mu = 0.5 * (elements[EL].calViscosity() + elements[ER].calViscosity());
			dudx = 0.5 * (elements[EL].gradvel[0] + elements[ER].gradvel[0]);
			dudy = 0.5 * (elements[EL].gradvel[1] + elements[ER].gradvel[1]);
			dvdx = 0.5 * (elements[EL].gradvel[2] + elements[ER].gradvel[2]);
			dvdy = 0.5 * (elements[EL].gradvel[3] + elements[ER].gradvel[3]);

			dTdx = 0.5 * (elements[EL].gradtemp[0] + elements[ER].gradtemp[0]);
			dTdy = 0.5 * (elements[EL].gradtemp[1] + elements[ER].gradtemp[1]);

			if ((*it).marker == INTERIOR) {
				dudx_f = dudx + (uL - uR - dudx*(xL-xR) - dudy*(yL-yR))*(xL-xR)/(r*r);
				dudy_f = dudy + (uL - uR - dudx*(xL-xR) - dudy*(yL-yR))*(yL-yR)/(r*r);
				dvdx_f = dvdx + (vL - vR - dvdx*(xL-xR) - dvdy*(yL-yR))*(xL-xR)/(r*r);
				dvdy_f = dvdy + (vL - vR - dvdx*(xL-xR) - dvdy*(yL-yR))*(yL-yR)/(r*r);
				dTdx_f = dTdx + (TL - TR - dTdx*(xL-xR) - dTdy*(yL-yR))*(xL-xR)/(r*r);
				dTdy_f = dTdy + (TL - TR - dTdx*(xL-xR) - dTdy*(yL-yR))*(yL-yR)/(r*r);
			}

			else {
				dudx_f = (uL - uR)*(xL-xR)/(r*r);
				dudy_f = (uL - uR)*(yL-yR)/(r*r);
				dvdx_f = (vL - vR)*(xL-xR)/(r*r);
				dvdy_f = (vL - vR)*(yL-yR)/(r*r);
				dTdx_f = (TL - TR)*(xL-xR)/(r*r);
				dTdy_f = (TL - TR)*(yL-yR)/(r*r);
			}

			txx = 2.0*mu*dudx_f - 2.0/3.0*mu*(dudx_f + dvdy_f);
			tyy = 2.0*mu*dvdy_f - 2.0/3.0*mu*(dudx_f + dvdy_f);
			txy = mu*(dudy_f + dvdx_f);
			K = mu*gam*R/((gam-1)*Pr);
			qx = -K*dTdx_f;
			qy = -K*dTdy_f;

			(*it).u_face = 0.5*(uL + uR);
			(*it).v_face = 0.5*(vL + vR);

	//		Subtract viscous fluxes from inviscid fluxes
			(*it).mom_x_flux -= ((*it).area) * (txx*((*it).nx) + txy*((*it).ny));
			(*it).mom_y_flux -= ((*it).area) * (txy*((*it).nx) + tyy*((*it).ny));
			(*it).energy_flux -= ((*it).area) * ((((*it).u_face)*txx + ((*it).v_face)*txy - qx)*((*it).nx) + (((*it).u_face)*txy + ((*it).v_face)*tyy - qy)*((*it).ny) );	
			
			//No Penetration
			if ((*it).marker == TOP || (*it).marker == BOTTOM || (*it).marker == OTHER) {
				(*it).mass_flux = 0.0;
				(*it).mom_x_flux = ((*it).area) * (D_plus*pR*((*it).nx) + D_minus*pL*((*it).nx) - (txx*((*it).nx) + txy*((*it).ny)) );
				(*it).mom_y_flux = ((*it).area) * (D_plus*pR*((*it).ny) + D_minus*pL*((*it).ny) - (txy*((*it).nx) + tyy*((*it).ny)) );
				(*it).energy_flux = -1.0*((*it).area) * ((-1.0*qx)*((*it).nx) + (-1.0*qy)*((*it).ny) );
			}

		#endif
		
		
	}	

	
	
}


				
#endif

