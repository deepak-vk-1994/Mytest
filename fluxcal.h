#ifndef _FLUXCAL_H
#define _FLUXCAL_H

#include <math.h>
#include <assert.h>
#include <float.h>
#include "grid.h"

inline double maximum(double a, double b, double c, double d) {
	return std::max(std::max(a,b),std::max(c,d));
}

inline double minimum(double a, double b, double c, double d) {
	return std::min(std::min(a,b),std::min(c,d));
}

inline double calGamma(double b, double sigma) {
	return (b/betastar - sigma*kappa*kappa/sqrt(betastar));
}

inline void calLimiters(int index, std::vector<Point > &points, std::vector<Element > &elements, std::vector<Face > &faces, double phi[6]) {
	double q[6], qmax[6], qmin[6], dot[6];
	double xmc, ymc, delmax, delmin,lim;
	// std::vector<double> q(4),qmax(4), qmin(4), dot(4);
	int n = 0;
	double epsilon,K,h; //Venkatakrishnan parameters
	h = sqrt(elements[index].volume);
	K = 5.0; //Change according to solution
	epsilon = (K*h)*(K*h)*(K*h);

	for (int i = 0; i < 6; i++) phi[i] = 0.0;
	if (index >= N_elem) {//Do nothing
	}
	else {
		std::vector<int> nbour(3);
		for (int i = 0; i < 3; i++) {
			n = faces[elements[index].face[i]].elementL;
			if (n == index)
				n = faces[elements[index].face[i]].elementR;
			nbour[i] = n;
			if (nbour[i] == index) cout << index;

			assert(nbour[i] != index);
		}

		q[0] = elements[index].rho; q[1] = elements[index].u; q[2] = elements[index].v; q[3] = elements[index].p;
		q[4] = elements[index].k; q[5] = elements[index].omega;

		qmax[0] = maximum(elements[nbour[0]].rho,elements[nbour[1]].rho,elements[nbour[2]].rho,elements[index].rho);
		qmax[1] = maximum(elements[nbour[0]].u,elements[nbour[1]].u,elements[nbour[2]].u,elements[index].u);
		qmax[2] = maximum(elements[nbour[0]].v,elements[nbour[1]].v,elements[nbour[2]].v,elements[index].v);
		qmax[3] = maximum(elements[nbour[0]].p,elements[nbour[1]].p,elements[nbour[2]].p,elements[index].p);
		qmax[4] = maximum(elements[nbour[0]].k,elements[nbour[1]].k,elements[nbour[2]].k,elements[index].k);
		qmax[5] = maximum(elements[nbour[0]].omega,elements[nbour[1]].omega,elements[nbour[2]].omega,elements[index].omega);

		qmin[0] = minimum(elements[nbour[0]].rho,elements[nbour[1]].rho,elements[nbour[2]].rho,elements[index].rho);
		qmin[1] = minimum(elements[nbour[0]].u,elements[nbour[1]].u,elements[nbour[2]].u,elements[index].u);
		qmin[2] = minimum(elements[nbour[0]].v,elements[nbour[1]].v,elements[nbour[2]].v,elements[index].v);
		qmin[3] = minimum(elements[nbour[0]].p,elements[nbour[1]].p,elements[nbour[2]].p,elements[index].p);
		qmin[4] = minimum(elements[nbour[0]].k,elements[nbour[1]].k,elements[nbour[2]].k,elements[index].k);
		qmin[5] = minimum(elements[nbour[0]].omega,elements[nbour[1]].omega,elements[nbour[2]].omega,elements[index].omega);

		for (int i = 0; i < 3; i++) {
			xmc = faces[elements[index].face[i]].xmp - elements[index].xc;
			ymc = faces[elements[index].face[i]].ymp - elements[index].yc;
			dot[0] = elements[index].gradrho[0]*xmc + elements[index].gradrho[1]*ymc;
			dot[1] = elements[index].gradu[0]*xmc + elements[index].gradu[1]*ymc;
			dot[2] = elements[index].gradv[0]*xmc + elements[index].gradv[1]*ymc;
			dot[3] = elements[index].gradp[0]*xmc + elements[index].gradp[1]*ymc;
			dot[4] = elements[index].gradk[0]*xmc + elements[index].gradk[1]*ymc;
			dot[5] = elements[index].gradomega[0]*xmc + elements[index].gradomega[1]*ymc;

			// xmc = points[elements[index].vertex[i]].x - elements[index].xc;
			// ymc = points[elements[index].vertex[i]].y - elements[index].yc;
			// dot[0] = elements[index].rho + (elements[index].gradrho[0]*xmc + elements[index].gradrho[1]*ymc) - q[0];
			// dot[1] = elements[index].u + (elements[index].gradu[0]*xmc + elements[index].gradu[1]*ymc) - q[1];
			// dot[2] = elements[index].v + (elements[index].gradv[0]*xmc + elements[index].gradv[1]*ymc) - q[2];
			// dot[3] = elements[index].p + (elements[index].gradp[0]*xmc + elements[index].gradp[1]*ymc) - q[3];
			// dot[4] = elements[index].k + (elements[index].gradk[0]*xmc + elements[index].gradk[1]*ymc) - q[4];
			// dot[5] = elements[index].omega + (elements[index].gradomega[0]*xmc + elements[index].gradomega[1]*ymc) - q[5];


			for (int j = 0; j < 6; j++) {
				delmax = qmax[j] - q[j];
				delmin = qmin[j] - q[j];


				if (dot[j] > 0)
				 	lim =  ( ((delmax*delmax + epsilon) + 2*dot[j]*delmax)/(delmax*delmax + 2*dot[j]*dot[j] + delmax*dot[j] + epsilon) );

				else if (dot[j] < 0)
					lim =  ( ((delmin*delmin + epsilon) + 2*dot[j]*delmin)/(delmin*delmin + 2*dot[j]*dot[j] + delmin*dot[j] + epsilon) );

				else
					lim = 1.0;

				if (i == 0) phi[j] = lim;
				else phi[j] = std::min(phi[j],lim);

			}

		}
	}
	
}

inline void calLimitersFace(int elem_index, int face_index, std::vector<Point > &points, std::vector<Element > &elements, std::vector<Face > &faces, double phi[6]) {
	double q[6], qmax[6], qmin[6], dot[6], distance1, distance2;
	double xmc, ymc, delmax, delmin,lim;
	int n = 0;
	double epsilon,K,h; //Venkatakrishnan parameters
	h = sqrt(elements[elem_index].volume);
	K = 5.0; //Change according to solution
	epsilon = (K*h)*(K*h)*(K*h);

	std::vector<int> nbour(3);
	for (int i = 0; i < 3; i++) {
		n = faces[elements[elem_index].face[i]].elementL;
		if (n == elem_index)
			n = faces[elements[elem_index].face[i]].elementR;
		nbour[i] = n;
		assert(nbour[i] != elem_index);
	}

	q[0] = elements[elem_index].rho; q[1] = elements[elem_index].u; q[2] = elements[elem_index].v; q[3] = elements[elem_index].p;
	q[4] = elements[elem_index].k; q[5] = elements[elem_index].omega;

	qmax[0] = maximum(elements[nbour[0]].rho,elements[nbour[1]].rho,elements[nbour[2]].rho,elements[elem_index].rho);
	qmax[1] = maximum(elements[nbour[0]].u,elements[nbour[1]].u,elements[nbour[2]].u,elements[elem_index].u);
	qmax[2] = maximum(elements[nbour[0]].v,elements[nbour[1]].v,elements[nbour[2]].v,elements[elem_index].v);
	qmax[3] = maximum(elements[nbour[0]].p,elements[nbour[1]].p,elements[nbour[2]].p,elements[elem_index].p);
	qmax[4] = maximum(elements[nbour[0]].k,elements[nbour[1]].k,elements[nbour[2]].k,elements[elem_index].k);
	qmax[5] = maximum(elements[nbour[0]].omega,elements[nbour[1]].omega,elements[nbour[2]].omega,elements[elem_index].omega);

	qmin[0] = minimum(elements[nbour[0]].rho,elements[nbour[1]].rho,elements[nbour[2]].rho,elements[elem_index].rho);
	qmin[1] = minimum(elements[nbour[0]].u,elements[nbour[1]].u,elements[nbour[2]].u,elements[elem_index].u);
	qmin[2] = minimum(elements[nbour[0]].v,elements[nbour[1]].v,elements[nbour[2]].v,elements[elem_index].v);
	qmin[3] = minimum(elements[nbour[0]].p,elements[nbour[1]].p,elements[nbour[2]].p,elements[elem_index].p);
	qmin[4] = minimum(elements[nbour[0]].k,elements[nbour[1]].k,elements[nbour[2]].k,elements[elem_index].k);
	qmin[5] = minimum(elements[nbour[0]].omega,elements[nbour[1]].omega,elements[nbour[2]].omega,elements[elem_index].omega);

	xmc = points[faces[face_index].vertex1].x - elements[elem_index].xc;
	ymc = points[faces[face_index].vertex1].y - elements[elem_index].yc;

	dot[0] = elements[elem_index].gradrho[0]*xmc + elements[elem_index].gradrho[1]*ymc;
	dot[1] = elements[elem_index].gradu[0]*xmc + elements[elem_index].gradu[1]*ymc;
	dot[2] = elements[elem_index].gradv[0]*xmc + elements[elem_index].gradv[1]*ymc;
	dot[3] = elements[elem_index].gradp[0]*xmc + elements[elem_index].gradp[1]*ymc;
	dot[4] = elements[elem_index].gradk[0]*xmc + elements[elem_index].gradk[1]*ymc;
	dot[5] = elements[elem_index].gradomega[0]*xmc + elements[elem_index].gradomega[1]*ymc;

	for (int j = 0; j < 6; j++) {
		delmax = qmax[j] - q[j];
		delmin = qmin[j] - q[j];
		if (dot[j] > 0)
		 	lim =  ( ((delmax*delmax + epsilon) + 2*dot[j]*delmax)/(delmax*delmax + 2*dot[j]*dot[j] + delmax*dot[j] + epsilon) );

		else if (dot[j] < 0)
			lim =  ( ((delmin*delmin + epsilon) + 2*dot[j]*delmin)/(delmin*delmin + 2*dot[j]*dot[j] + delmin*dot[j] + epsilon) );

		else
			lim = 1.0;

		phi[j] = lim;
	}

	xmc = points[faces[face_index].vertex2].x - elements[elem_index].xc;
	ymc = points[faces[face_index].vertex2].y - elements[elem_index].yc;

	dot[0] = elements[elem_index].gradrho[0]*xmc + elements[elem_index].gradrho[1]*ymc;
	dot[1] = elements[elem_index].gradu[0]*xmc + elements[elem_index].gradu[1]*ymc;
	dot[2] = elements[elem_index].gradv[0]*xmc + elements[elem_index].gradv[1]*ymc;
	dot[3] = elements[elem_index].gradp[0]*xmc + elements[elem_index].gradp[1]*ymc;
	dot[4] = elements[elem_index].gradk[0]*xmc + elements[elem_index].gradk[1]*ymc;
	dot[5] = elements[elem_index].gradomega[0]*xmc + elements[elem_index].gradomega[1]*ymc;

	for (int j = 0; j < 6; j++) {
		delmax = qmax[j] - q[j];
		delmin = qmin[j] - q[j];
		if (dot[j] > 0)
		 	lim =  ( ((delmax*delmax + epsilon) + 2*dot[j]*delmax)/(delmax*delmax + 2*dot[j]*dot[j] + delmax*dot[j] + epsilon) );

		else if (dot[j] < 0)
			lim =  ( ((delmin*delmin + epsilon) + 2*dot[j]*delmin)/(delmin*delmin + 2*dot[j]*dot[j] + delmin*dot[j] + epsilon) );

		else
			lim = 1.0;

		if (lim < phi[j])
			phi[j] = lim;
	}
}

void calFlux(std::vector<Point > &points, std::vector<Element > &elements, std::vector<Face > &faces) {
	int EL, ER;

	//Liou steffan variables
	double ML, MR, a_face;
	double alpha_plus, alpha_minus, beta_plus, beta_minus, C_VL_plus, C_VL_minus, C_LS_plus, C_LS_minus;
	double D_plus, D_minus;
	double E_t_plus, E_t_minus;

	//Viscous flux variables
	double txx,tyy,txy,mu,K,qx,qy,dudx,dudy,dvdx,dvdy,dTdx,dTdy,dudx_f,dudy_f,dvdx_f,dvdy_f,dTdx_f,dTdy_f;
	double distance1, distance2;
	// double tauxx,tauyy,tauxy;

	//Others
	double uL,uR,vL,vR,rhoL,rhoR,pL,pR,TL,TR;
	double xL,xR,yL,yR,r,xmcL,ymcL,xmcR,ymcR;
	double phiL[6],phiR[6]; //Limiters for each primitive variable
	int num = 0;
	
	//Turbulence
	double kL,kR,omegaL,omegaR,mu_t,g,lamda,F1;
	double dkdx,dkdy,domegadx,domegady,dkdx_f,dkdy_f,domegadx_f,domegady_f;
	double turbphik,turbphiomega,P_k,P_w,D_k,D_w;

	for (int i = 0; i < N_elem; i++) {

		elements[i].gradu[0] = elements[i].gradu[1] = 0.0;
		elements[i].gradv[0] = elements[i].gradv[1] = 0.0;
		elements[i].gradtemp[0] = elements[i].gradtemp[1] = 0.0;
		elements[i].gradrho[0] = elements[i].gradrho[1] = 0.0;
		elements[i].gradp[0] = elements[i].gradp[1] = 0.0;
		elements[i].gradk[0] = elements[i].gradk[1] = 0.0;
		elements[i].gradomega[0] = elements[i].gradomega[1] = 0.0;
		elements[i].gradentropy[0] = elements[i].gradentropy[1] = 0.0;
		

		for (int j = 0; j < 3; j++) {
			num = elements[i].face[j];
			EL = faces[num].elementL;
			ER = faces[num].elementR;
			
			faces[num].rho_face = 0.5 * (elements[EL].rho + elements[ER].rho);
			faces[num].u_face = 0.5 * (elements[EL].u + elements[ER].u);
			faces[num].v_face = 0.5 * (elements[EL].v + elements[ER].v);
			faces[num].p_face = 0.5 * (elements[EL].p + elements[ER].p);
			faces[num].k_face = 0.5 * (elements[EL].k + elements[ER].k);
			faces[num].omega_face = 0.5 * (elements[EL].omega + elements[ER].omega);
			faces[num].T_face = 0.5 * (elements[EL].calTemp() + elements[ER].calTemp());
			faces[num].en_face = 0.5 * (elements[EL].calEntropy() + elements[ER].calEntropy());

			if (i == ER) {
				elements[i].gradu[0] += faces[num].area/elements[i].volume * faces[num].u_face*faces[num].nx;
				elements[i].gradu[1] += faces[num].area/elements[i].volume * faces[num].u_face*faces[num].ny;
				elements[i].gradv[0] += faces[num].area/elements[i].volume * faces[num].v_face*faces[num].nx;
				elements[i].gradv[1] += faces[num].area/elements[i].volume * faces[num].v_face*faces[num].ny;
	        	elements[i].gradtemp[0] += faces[num].area/elements[i].volume * faces[num].T_face*faces[num].nx;
				elements[i].gradtemp[1] += faces[num].area/elements[i].volume * faces[num].T_face*faces[num].ny;
				elements[i].gradrho[0] += faces[num].area/elements[i].volume * faces[num].rho_face*faces[num].nx;
				elements[i].gradrho[1] += faces[num].area/elements[i].volume * faces[num].rho_face*faces[num].ny;
				elements[i].gradp[0] += faces[num].area/elements[i].volume * faces[num].p_face*faces[num].nx;
				elements[i].gradp[1] += faces[num].area/elements[i].volume * faces[num].p_face*faces[num].ny;
				elements[i].gradk[0] += faces[num].area/elements[i].volume * faces[num].k_face*faces[num].nx;
				elements[i].gradk[1] += faces[num].area/elements[i].volume * faces[num].k_face*faces[num].ny;
				elements[i].gradomega[0] += faces[num].area/elements[i].volume * faces[num].omega_face*faces[num].nx;
				elements[i].gradomega[1] += faces[num].area/elements[i].volume * faces[num].omega_face*faces[num].ny;
				elements[i].gradentropy[0] += faces[num].area/elements[i].volume * faces[num].en_face*faces[num].nx;
				elements[i].gradentropy[1] += faces[num].area/elements[i].volume * faces[num].en_face*faces[num].ny;
            }

            else if (i == EL) {
				elements[i].gradu[0] += -1.0*faces[num].area/elements[i].volume * faces[num].u_face*faces[num].nx;
				elements[i].gradu[1] += -1.0*faces[num].area/elements[i].volume * faces[num].u_face*faces[num].ny;
				elements[i].gradv[0] += -1.0*faces[num].area/elements[i].volume * faces[num].v_face*faces[num].nx;
				elements[i].gradv[1] += -1.0*faces[num].area/elements[i].volume * faces[num].v_face*faces[num].ny;
	        	elements[i].gradtemp[0] += -1.0*faces[num].area/elements[i].volume * faces[num].T_face*faces[num].nx;
				elements[i].gradtemp[1] += -1.0*faces[num].area/elements[i].volume * faces[num].T_face*faces[num].ny;
				elements[i].gradrho[0] += -1.0*faces[num].area/elements[i].volume * faces[num].rho_face*faces[num].nx;
				elements[i].gradrho[1] += -1.0*faces[num].area/elements[i].volume * faces[num].rho_face*faces[num].ny;
				elements[i].gradp[0] += -1.0*faces[num].area/elements[i].volume * faces[num].p_face*faces[num].nx;
				elements[i].gradp[1] += -1.0*faces[num].area/elements[i].volume * faces[num].p_face*faces[num].ny;
				elements[i].gradk[0] += -1.0*faces[num].area/elements[i].volume * faces[num].k_face*faces[num].nx;
				elements[i].gradk[1] += -1.0*faces[num].area/elements[i].volume * faces[num].k_face*faces[num].ny;
				elements[i].gradomega[0] += -1.0*faces[num].area/elements[i].volume * faces[num].omega_face*faces[num].nx;
				elements[i].gradomega[1] += -1.0*faces[num].area/elements[i].volume * faces[num].omega_face*faces[num].ny;
				elements[i].gradentropy[0] += -1.0*faces[num].area/elements[i].volume * faces[num].en_face*faces[num].nx;
				elements[i].gradentropy[1] += -1.0*faces[num].area/elements[i].volume * faces[num].en_face*faces[num].ny;
            }
			
		}

		calLimiters(i,points,elements,faces,elements[i].phi);

		if (turbulence == SST) {
			elements[i].calF1();
			elements[i].calTurbulentVisc();
			D_k = betastar*elements[i].rho*elements[i].omega*elements[i].k;
			P_k = elements[i].mu_t*elements[i].calVorticity()*elements[i].calVorticity();
			P_k = std::min(P_k,20.0*D_k);	
			elements[i].source1 = P_k - D_k;
			elements[i].pk = P_k;
			elements[i].dk = D_k;
			g = elements[i].F1*calGamma(beta1,turbphi_omega1) + (1 - elements[i].F1)*calGamma(beta2,turbphi_omega2);  
			elements[i].beta = elements[i].F1*beta1 + (1.0 - elements[i].F1)*beta2;    
			lamda = (1.0 - elements[i].F1)*elements[i].calCrossDiffusion();     
			D_w = elements[i].beta*elements[i].rho*elements[i].omega*elements[i].omega;
			P_w = std::min(g*elements[i].rho*P_k/elements[i].mu_t,20.0*D_w);	
			elements[i].source2 = P_w - D_w + lamda;
			elements[i].pw = P_w;
			elements[i].dw = D_w;
			elements[i].lam = lamda;
		}
	}


/*------------------------------------------------------------------------------------------------------------------------------------*/

	//FLUX CALCULATION
	for (std::vector<Face >::iterator it = faces.begin(); it != faces.end(); it++) {
		EL = (*it).elementL;
		ER = (*it).elementR;		
		
		xmcL = (*it).xmp - elements[EL].xc; ymcL = (*it).ymp - elements[EL].yc;
		xmcR = (*it).xmp - elements[ER].xc; ymcR = (*it).ymp - elements[ER].yc;

		for (int i = 0; i < 6; i++){
			phiL[i] = elements[EL].phi[i];
			phiR[i] = elements[ER].phi[i];
		}
		
		//Reconstruction of primitive variables
		//First Order

		// calLimitersFace(EL,it - faces.begin(),points,elements,faces,phiL);
		// calLimitersFace(ER,it - faces.begin(),points,elements,faces,phiR);

		if (ord_accuracy == FO) {
			rhoL = elements[EL].rho;
			rhoR = elements[ER].rho;
			uL = elements[EL].u;
			uR = elements[ER].u;
			vL = elements[EL].v;
			vR = elements[ER].v;
			pL = elements[EL].p;
			pR = elements[ER].p;
			kL = elements[EL].k;
			kR = elements[ER].k;
			omegaL = elements[EL].omega;
			omegaR = elements[ER].omega;
		}
		//Second Order
		if (ord_accuracy == HO) {
			rhoL = elements[EL].rho + phiL[0]*(elements[EL].gradrho[0]*xmcL + elements[EL].gradrho[1]*ymcL);
			rhoR = elements[ER].rho + phiR[0]*(elements[ER].gradrho[0]*xmcR + elements[ER].gradrho[1]*ymcR);
			uL = elements[EL].u + phiL[1]*(elements[EL].gradu[0]*xmcL + elements[EL].gradu[1]*ymcL);
			uR = elements[ER].u + phiR[1]*(elements[ER].gradu[0]*xmcR + elements[ER].gradu[1]*ymcR);
			vL = elements[EL].v + phiL[2]*(elements[EL].gradv[0]*xmcL + elements[EL].gradv[1]*ymcL);
			vR = elements[ER].v + phiR[2]*(elements[ER].gradv[0]*xmcR + elements[ER].gradv[1]*ymcR);
			pL = elements[EL].p + phiL[3]*(elements[EL].gradp[0]*xmcL + elements[EL].gradp[1]*ymcL);
			pR = elements[ER].p + phiR[3]*(elements[ER].gradp[0]*xmcR + elements[ER].gradp[1]*ymcR);
			kL = elements[EL].k + phiL[4]*(elements[EL].gradk[0]*xmcL + elements[EL].gradk[1]*ymcL);
			kR = elements[ER].k + phiR[4]*(elements[ER].gradk[0]*xmcR + elements[ER].gradk[1]*ymcR);
			omegaL = elements[EL].omega + phiL[5]*(elements[EL].gradomega[0]*xmcL + elements[EL].gradomega[1]*ymcL);
			omegaR = elements[ER].omega + phiR[5]*(elements[ER].gradomega[0]*xmcR + elements[ER].gradomega[1]*ymcR);
		}

		TL = pL/(rhoL*R);
		TR = pR/(rhoR*R);
		E_t_plus = gam/(gam-1)*(pR/rhoR) + 0.5*(uR*uR + vR*vR);
		E_t_minus = gam/(gam-1)*(pL/rhoL) + 0.5*(uL*uL + vL*vL);	
		
		a_face = 0.5 * ( sqrt(gam*(pL/rhoL)) + sqrt(gam*(pR/rhoR)) );
		ML = (uL*((*it).nx) + vL*((*it).ny))/a_face;
		MR = (uR*((*it).nx) + vR*((*it).ny))/a_face;

		//INVISCID FLUXES
		if (in_flux == AUSM) {
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
	
			(*it).mass_flux = ((*it).area) * (rhoR*a_face*C_LS_plus*1.0 + rhoL*a_face*C_LS_minus*1.0);
			(*it).mom_x_flux = ((*it).area) * (rhoR*a_face*C_LS_plus*uR + rhoL*a_face*C_LS_minus*uL + D_plus*pR*((*it).nx) + D_minus*pL*((*it).nx));
			(*it).mom_y_flux = ((*it).area) * (rhoR*a_face*C_LS_plus*vR + rhoL*a_face*C_LS_minus*vL + D_plus*pR*((*it).ny) + D_minus*pL*((*it).ny));
			(*it).energy_flux = ((*it).area) * (rhoR*a_face*C_LS_plus*E_t_plus + rhoL*a_face*C_LS_minus*E_t_minus);

			//Turbulent fluxes
			if (turbulence == SST) {
				(*it).k_flux = ((*it).area) * (rhoR*a_face*C_LS_plus*kR + rhoL*a_face*C_LS_minus*kL);
				(*it).omega_flux = ((*it).area) * (rhoR*a_face*C_LS_plus*omegaR + rhoL*a_face*C_LS_minus*omegaL);
			}
			
        	if (visc_flux == NO_VISC) {
	        	if ((*it).marker != INTERIOR) {
	        		for (bcs_it = bcs.begin(); bcs_it != bcs.end(); bcs_it++) {
						if ((*it).marker == (*bcs_it)) {
						 	(*it).mass_flux = 0.0;
						 	(*it).mom_x_flux = ((*it).area) * (D_plus*pR*((*it).nx) + D_minus*pL*((*it).nx) );
						 	(*it).mom_y_flux = ((*it).area) * (D_plus*pR*((*it).ny) + D_minus*pL*((*it).ny) );
						 	(*it).energy_flux = 0.0;
						 	if (turbulence == SST) {
							 	(*it).k_flux = 0.0;
							 	(*it).omega_flux = 0.0;
							 }
						 }
	        		}
	        	}
            }
		
		}	

		//VISCOUS FLUXES
		if (visc_flux == HO) {
			//Redefine variables as first order
			rhoL = elements[EL].rho;
			rhoR = elements[ER].rho;
			uL = elements[EL].u;
			uR = elements[ER].u;
			vL = elements[EL].v;
			vR = elements[ER].v;
			pL = elements[EL].p;
			pR = elements[ER].p;
			kL = elements[EL].k;
			kR = elements[ER].k;
			omegaL = elements[EL].omega;
			omegaR = elements[ER].omega;

			xL = elements[EL].xc;
			yL = elements[EL].yc;		
			xR = elements[ER].xc;
			yR = elements[ER].yc;
			r = sqrt((xR-xL)*(xR-xL) + (yR-yL)*(yR-yL));  //Vector magnitude from centroidL to centroidR

			mu = 0.5 * (elements[EL].calViscosity() + elements[ER].calViscosity());

			dudx = 0.5 * (elements[EL].gradu[0] + elements[ER].gradu[0]);
			dudy = 0.5 * (elements[EL].gradu[1] + elements[ER].gradu[1]);
			dvdx = 0.5 * (elements[EL].gradv[0] + elements[ER].gradv[0]);
			dvdy = 0.5 * (elements[EL].gradv[1] + elements[ER].gradv[1]);

			dTdx = 0.5 * (elements[EL].gradtemp[0] + elements[ER].gradtemp[0]);
			dTdy = 0.5 * (elements[EL].gradtemp[1] + elements[ER].gradtemp[1]);

			dudx_f = dudx + (uL - uR - dudx*(xL-xR) - dudy*(yL-yR))*(xL-xR)/(r*r);
			dudy_f = dudy + (uL - uR - dudx*(xL-xR) - dudy*(yL-yR))*(yL-yR)/(r*r);
			dvdx_f = dvdx + (vL - vR - dvdx*(xL-xR) - dvdy*(yL-yR))*(xL-xR)/(r*r);
			dvdy_f = dvdy + (vL - vR - dvdx*(xL-xR) - dvdy*(yL-yR))*(yL-yR)/(r*r);
			dTdx_f = dTdx + (TL - TR - dTdx*(xL-xR) - dTdy*(yL-yR))*(xL-xR)/(r*r);
			dTdy_f = dTdy + (TL - TR - dTdx*(xL-xR) - dTdy*(yL-yR))*(yL-yR)/(r*r);

			//Turbulence gradients
			dkdx = 0.5 * (elements[EL].gradk[0] + elements[ER].gradk[0]);
			dkdy = 0.5 * (elements[EL].gradk[1] + elements[ER].gradk[1]);
			domegadx = 0.5 * (elements[EL].gradomega[0] + elements[ER].gradomega[0]);
			domegady = 0.5 * (elements[EL].gradomega[1] + elements[ER].gradomega[1]);

			dkdx_f = dkdx + (kL - kR - dkdx*(xL-xR) - dkdy*(yL-yR))*(xL-xR)/(r*r);
			dkdy_f = dkdy + (kL - kR - dkdx*(xL-xR) - dkdy*(yL-yR))*(yL-yR)/(r*r);
			domegadx_f = domegadx + (omegaL - omegaR - domegadx*(xL-xR) - domegady*(yL-yR))*(xL-xR)/(r*r);
			domegady_f = domegady + (omegaL - omegaR - domegadx*(xL-xR) - domegady*(yL-yR))*(yL-yR)/(r*r);

			if (turbulence == SST) {
				F1 = 0.5*(elements[EL].F1 + elements[ER].F1);
				mu_t = 0.5*(elements[EL].mu_t + elements[ER].mu_t);

				turbphik = F1*turbphi_k1 + (1.0-F1)*turbphi_k2;
				turbphiomega = F1*turbphi_omega1 + (1.0-F1)*turbphi_omega2;
			}
			else {
				mu_t = 0.0;
			}

			txx = 2.0*mu*dudx_f - 2.0/3.0*mu*(dudx_f + dvdy_f) + 2.0*mu_t*dudx_f - 2.0/3.0*mu_t*(dudx_f + dvdy_f);
			tyy = 2.0*mu*dvdy_f - 2.0/3.0*mu*(dudx_f + dvdy_f) + 2.0*mu_t*dvdy_f - 2.0/3.0*mu_t*(dudx_f + dvdy_f);
			txy = mu*(dudy_f + dvdx_f) + mu_t*(dudy_f + dvdx_f);
			
			K = mu*gam*R/((gam-1)*Pr) + mu_t*gam*R/((gam-1)*Pr_t);
			qx = -K*dTdx_f;
			qy = -K*dTdy_f;

			(*it).u_face = 0.5*(uL + uR);
			(*it).v_face = 0.5*(vL + vR);

			(*it).mom_x_flux -= ((*it).area) * (txx*((*it).nx) + txy*((*it).ny));
			(*it).mom_y_flux -= ((*it).area) * (txy*((*it).nx) + tyy*((*it).ny));
			(*it).energy_flux -= ((*it).area) * ((((*it).u_face)*txx + ((*it).v_face)*txy - qx)*((*it).nx) + (((*it).u_face)*txy + ((*it).v_face)*tyy - qy)*((*it).ny) );	
			if (turbulence == SST) {
				(*it).k_flux -= ((*it).area) * ((mu + turbphik*mu_t)*dkdx_f*((*it).nx) + (mu + turbphik*mu_t)*dkdy_f*((*it).ny));
				(*it).omega_flux -= ((*it).area) * ((mu + turbphiomega*mu_t)*domegadx_f*((*it).nx) + (mu + turbphiomega*mu_t)*domegady_f*((*it).ny));
			}
			//No Penetration
			if ((*it).marker != INTERIOR) {
        		for (bcs_it = bcs.begin(); bcs_it != bcs.end(); bcs_it++) {
					if ((*it).marker == (*bcs_it)) {
						(*it).mass_flux = 0.0;
						(*it).mom_x_flux = ((*it).area) * (D_plus*pR*((*it).nx) + D_minus*pL*((*it).nx) - (txx*((*it).nx) + txy*((*it).ny)) );
						(*it).mom_y_flux = ((*it).area) * (D_plus*pR*((*it).ny) + D_minus*pL*((*it).ny) - (txy*((*it).nx) + tyy*((*it).ny)) );
						(*it).energy_flux = -1.0*((*it).area) * ((-1.0*qx)*((*it).nx) + (-1.0*qy)*((*it).ny));
						if (turbulence == SST) {
							(*it).k_flux = -1.0*((*it).area) * ((mu + turbphik*mu_t)*dkdx_f*((*it).nx) + (mu + turbphik*mu_t)*dkdy_f*((*it).ny));
							(*it).omega_flux = -1.0*((*it).area) * ((mu + turbphiomega*mu_t)*domegadx_f*((*it).nx) + (mu + turbphiomega*mu_t)*domegady_f*((*it).ny));
						}
					}
				}
			}
		}		
	}	
}
				
#endif


