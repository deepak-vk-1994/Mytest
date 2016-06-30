#ifndef _BOUNDARY_H
#define _BOUNDARY_H

#include <math.h>
#include <float.h>
#include "grid.h"
#include "fluxcal.h"

inline void applyBCface(std::vector<Element > &elements, int indexL, int indexR, int BCno,int mach_switch, double nx, double ny) {
	if (mach_switch == SUBSONIC) {
		if (BCno == INLET) {
			elements[indexL].p = elements[indexR].p;
			elements[indexL].u = u_inf;
			elements[indexL].v = v_inf;
			elements[indexL].w = w_inf;
			elements[indexL].rho = rho_inf;
			elements[indexL].k = k_inf;
			elements[indexL].omega = omega_inf;
			elements[indexL].mu_t = elements[indexR].mu_t;
			elements[indexL].F1 = elements[indexR].F1;
		}
		else if (BCno ==OUTLET) {
			elements[indexL].p = p_inf;
			elements[indexL].u = elements[indexR].u;
			elements[indexL].v = elements[indexR].v;
			elements[indexL].w = elements[indexR].w;
			elements[indexL].rho = elements[indexR].rho;
			elements[indexL].k = elements[indexR].k;
			elements[indexL].omega = elements[indexR].omega;
			elements[indexL].mu_t = elements[indexR].mu_t;
			elements[indexL].F1 = elements[indexR].F1;
		}
		else if (BCno == SLIP) {
			elements[indexL].rho = elements[indexR].rho;
			elements[indexL].p = elements[indexR].p;
			elements[indexL].u = elements[indexR].u - 2.0 *( elements[indexR].u*(nx)*(nx) + elements[indexR].v*(nx)*(ny) );
			elements[indexL].v = elements[indexR].v - 2.0 *( elements[indexR].u*(nx)*(ny) + elements[indexR].v*(ny)*(ny) );
			elements[indexL].k = elements[indexR].k;
			elements[indexL].omega = elements[indexR].omega;
			elements[indexL].mu_t = elements[indexR].mu_t;
			elements[indexL].F1 = elements[indexR].F1;
		}
		else if (BCno == NO_SLIP) {
			double omega_wall = 60.0*elements[indexR].calViscosity()/(beta1*elements[indexR].rho*4.0*elements[indexR].dist*elements[indexR].dist);
			elements[indexL].rho = elements[indexR].rho;
			elements[indexL].p = elements[indexR].p;
			elements[indexL].u = -elements[indexR].u;
			elements[indexL].v = -elements[indexR].v;
			elements[indexL].k = -elements[indexR].k;
			elements[indexL].mu_t = -elements[indexR].mu_t;
			elements[indexL].F1 = elements[indexR].F1;
			elements[indexL].omega = 2.0*omega_wall - elements[indexR].omega;
		}
	}

	else if (mach_switch == SUPERSONIC) {
		if (BCno == INLET) {
			elements[indexL].p = p_inf;
			elements[indexL].u = u_inf;
			elements[indexL].v = v_inf;
			elements[indexL].w = w_inf;
			elements[indexL].rho = rho_inf;
			elements[indexL].k = k_inf;
			elements[indexL].omega = omega_inf;
			elements[indexL].mu_t = elements[indexR].mu_t;
			elements[indexL].F1 = elements[indexR].F1;
		}
		else if (BCno ==OUTLET) {
			elements[indexL].p = elements[indexR].p;
			elements[indexL].u = elements[indexR].u;
			elements[indexL].v = elements[indexR].v;
			elements[indexL].w = elements[indexR].w;
			elements[indexL].rho = elements[indexR].rho;
			elements[indexL].k = elements[indexR].k;
			elements[indexL].omega = elements[indexR].omega;
			elements[indexL].mu_t = elements[indexR].mu_t;
			elements[indexL].F1 = elements[indexR].F1;
		}
		else if (BCno == SLIP) {
			elements[indexL].rho = elements[indexR].rho;
			elements[indexL].p = elements[indexR].p;
			elements[indexL].u = elements[indexR].u - 2.0 *( elements[indexR].u*(nx)*(nx) + elements[indexR].v*(nx)*(ny) );
			elements[indexL].v = elements[indexR].v - 2.0 *( elements[indexR].u*(nx)*(ny) + elements[indexR].v*(ny)*(ny) );
			elements[indexL].k = elements[indexR].k;
			elements[indexL].omega = elements[indexR].omega;
			elements[indexL].mu_t = elements[indexR].mu_t;
			elements[indexL].F1 = elements[indexR].F1;
		}
		else if (BCno == NO_SLIP) {
			double omega_wall = 60.0*elements[indexR].calViscosity()/(beta1*elements[indexR].rho*4.0*elements[indexR].dist*elements[indexR].dist);
			elements[indexL].rho = elements[indexR].rho;
			elements[indexL].p = elements[indexR].p;
			elements[indexL].u = -elements[indexR].u;
			elements[indexL].v = -elements[indexR].v;
			elements[indexL].k = -elements[indexR].k;
			elements[indexL].mu_t = -elements[indexR].mu_t;
			elements[indexL].F1 = elements[indexR].F1;
			elements[indexL].omega = 2.0*omega_wall - elements[indexR].omega;
		}
	}
}
void applyBC(std::vector<Point > &points, std::vector<Element > &elements, std::vector<Face > &faces) {	
	for (std::vector<Face >::iterator it = faces.begin(); it != faces.end(); it++) {
		if ((*it).marker == INTERIOR) continue;
		//INLET and OUTLET
		for (int i = 0; i < 2; i++) {
			elements[(*it).elementL].gradu[i] = -elements[(*it).elementR].gradu[i];
			elements[(*it).elementL].gradv[i] = -elements[(*it).elementR].gradv[i];
			elements[(*it).elementL].gradp[i] = -elements[(*it).elementR].gradp[i];
			elements[(*it).elementL].gradrho[i] = -elements[(*it).elementR].gradrho[i];
			elements[(*it).elementL].gradtemp[i] = -elements[(*it).elementR].gradtemp[i];
			//Turbulence
			elements[(*it).elementL].gradk[i] = -elements[(*it).elementR].gradk[i];
			elements[(*it).elementL].gradomega[i] = -elements[(*it).elementR].gradomega[i];
		}
		
		if ((*it).marker == LEFT) {
			applyBCface(elements,(*it).elementL,(*it).elementR,BCL,mach_switch,(*it).nx,(*it).ny);
		}
		else if ((*it).marker == RIGHT) {
			applyBCface(elements,(*it).elementL,(*it).elementR,BCR,mach_switch,(*it).nx,(*it).ny);
		}
		else if ((*it).marker == TOP) {
			applyBCface(elements,(*it).elementL,(*it).elementR,BCT,mach_switch,(*it).nx,(*it).ny);
		}
		else if ((*it).marker == BOTTOM) {
			applyBCface(elements,(*it).elementL,(*it).elementR,BCB,mach_switch,(*it).nx,(*it).ny);
		}
		else if ((*it).marker == OTHER) {
			applyBCface(elements,(*it).elementL,(*it).elementR,BCO,mach_switch,(*it).nx,(*it).ny);
		}
	}

	if (debug == 1) {
		ofstream file("./Debug/BC.txt");
		for (std::vector<Element >::iterator it = elements.begin(); it != elements.end(); it++) {
			file << (*it).rho <<" "<<(*it).p<<" "<<(*it).u<<" "<<(*it).v<<std::endl;
		}
	}
}

void calLocalTimeStep(std::vector<Element > &elements,std::vector<Face > &faces) {
	 double lambda,a;
	 Face e; 
	 int EL,ER;

	 for (int i = 0; i < N_elem; i++) {
		 lambda = 0;
		 for (int k = 0; k < 3; k++) {
		 	e = faces[elements[i].face[k]];
		 	EL = (e).elementL;
		 	ER = (e).elementR;
		 	a = 0.5 * ( sqrt(gam*(elements[EL].p/elements[EL].rho)) + sqrt(gam*(elements[ER].p/elements[ER].rho)) );
		 	lambda += (fabs(elements[i].u*(e.nx) + elements[i].v*(e.ny)) + a)*(e.area);
		 }
		 elements[i].delta_t = elements[i].volume*1.0*CFL/lambda;
	}
}

void computeResNorm(double &resnorm_laminar, double &resnorm_turbulent,std::vector<Element > &elements,int i) {
	resnorm_laminar += ( (elements[i].residue[0]/(rho_inf*u_inf))*(elements[i].residue[0]/(rho_inf*u_inf)) + (elements[i].residue[1]/(rho_inf*u_inf*u_inf))*(elements[i].residue[1]/(rho_inf*u_inf*u_inf)) + (elements[i].residue[2]/(rho_inf*u_inf*u_inf))*(elements[i].residue[2]/(rho_inf*u_inf*u_inf)) + (elements[i].residue[3]/(rho_inf*rho_inf*u_inf*u_inf*u_inf))*(elements[i].residue[3]/(rho_inf*rho_inf*u_inf*u_inf*u_inf)) );
	resnorm_turbulent += ( ((elements[i].residue[4]- elements[i].source1*elements[i].volume)/(rho_inf*u_inf*k_inf))*((elements[i].residue[4]- elements[i].source1*elements[i].volume)/(rho_inf*u_inf*k_inf)) + ((elements[i].residue[5]- elements[i].source2*elements[i].volume)/(rho_inf*u_inf*omega_inf))*((elements[i].residue[5]- elements[i].source2*elements[i].volume)/(rho_inf*u_inf*omega_inf)) );
}

void computeResidue(std::vector<Element > &elements,std::vector<Face > &faces) {
	for (int i = 0; i < N_elem; i++) {
		for (int k = 0; k < 6; k++) elements[i].residue[k] = 0.0;
		for (int j = 0; j < 3; j++) {
			if (faces[elements[i].face[j]].elementR == i) {
				elements[i].residue[0] += faces[elements[i].face[j]].mass_flux;
				elements[i].residue[1] += faces[elements[i].face[j]].mom_x_flux;
				elements[i].residue[2] += faces[elements[i].face[j]].mom_y_flux;
				elements[i].residue[3] += faces[elements[i].face[j]].energy_flux;
				elements[i].residue[4] += faces[elements[i].face[j]].k_flux;
				elements[i].residue[5] += faces[elements[i].face[j]].omega_flux;
			}
			else if (faces[elements[i].face[j]].elementL == i) {
				elements[i].residue[0] -= faces[elements[i].face[j]].mass_flux;
				elements[i].residue[1] -= faces[elements[i].face[j]].mom_x_flux;
				elements[i].residue[2] -= faces[elements[i].face[j]].mom_y_flux;
				elements[i].residue[3] -= faces[elements[i].face[j]].energy_flux;
				elements[i].residue[4] -= faces[elements[i].face[j]].k_flux;
				elements[i].residue[5] -= faces[elements[i].face[j]].omega_flux;
			}
			else 
				cout<<"ERROR in face-triangle indexing in UPDATE"<<endl;
		}
	}
}

void getE(std::vector<Element > &elements, std::vector<Face > &faces,int index) {
	double delta_t = 0.0;

	for (int i = 0; i < N_elem; i++) {
		if (index == 0) {
			elements[i].rho_temp = elements[i].rho;
			elements[i].u_temp = elements[i].u;
			elements[i].v_temp = elements[i].v;
			elements[i].p_temp = elements[i].p;
			elements[i].k_temp = elements[i].k;
			elements[i].omega_temp = elements[i].omega;
		}
		delta_t = elements[i].delta_t;
		
		if (index == 0) delta_t *= 1.0/2.0;
		if (index == 1) delta_t *= 1.0/2.0;
		if (index == 2) delta_t *= 1.0;
		if (index == 3) delta_t *= 1.0;	
		
		//Modified RK4
		// if (index == 0) delta_t *= 1.0/4.0;
		// if (index == 1) delta_t *= 1.0/3.0;
		// if (index == 2) delta_t *= 1.0/2.0;
		// if (index == 3) delta_t *= 1.0;

		elements[i].E[index][0] = (-1.0*delta_t/elements[i].volume) * elements[i].residue[0];
		elements[i].E[index][1] = (-1.0*delta_t/elements[i].volume) * (-1.0*elements[i].u_temp*elements[i].residue[0]/elements[i].rho_temp + elements[i].residue[1]/elements[i].rho_temp);
		elements[i].E[index][2] = (-1.0*delta_t/elements[i].volume) * (-1.0*elements[i].v_temp*elements[i].residue[0]/elements[i].rho_temp+ elements[i].residue[2]/elements[i].rho_temp);
		elements[i].E[index][3] = (-1.0*delta_t/elements[i].volume) * (0.5*(gam-1)*(elements[i].u_temp*elements[i].u_temp + elements[i].v_temp*elements[i].v_temp)*elements[i].residue[0] - (gam-1)*elements[i].u_temp*elements[i].residue[1] - (gam-1)*elements[i].v_temp*elements[i].residue[2] + (gam-1)*elements[i].residue[3]);
		elements[i].E[index][4] = (-1.0*delta_t/elements[i].volume) * (-1.0*elements[i].k_temp*elements[i].residue[0]/elements[i].rho_temp + 1.0/(1.0 + betastar*elements[i].omega_temp*delta_t) * (elements[i].residue[4] - elements[i].source1*elements[i].volume)/elements[i].rho_temp);
		elements[i].E[index][5] = (-1.0*delta_t/elements[i].volume) * (-1.0*elements[i].omega_temp*elements[i].residue[0]/elements[i].rho_temp + 1.0/(1.0 + 2.0*elements[i].beta*elements[i].omega_temp*delta_t) * (elements[i].residue[5]- elements[i].source2*elements[i].volume)/elements[i].rho_temp);
	}
}

void update(std::vector<Element > &elements,std::vector<Face > &faces,int stage) {
	for (int i = 0; i < N_elem; i++) {
		elements[i].rho = elements[i].rho_temp + elements[i].E[stage][0];
		elements[i].u = elements[i].u_temp + elements[i].E[stage][1];
		elements[i].v = elements[i].v_temp + elements[i].E[stage][2];
		elements[i].p = elements[i].p_temp + elements[i].E[stage][3];
		if (turbulence != NO_TURB) {
			elements[i].k = elements[i].k_temp + elements[i].E[stage][4];
			elements[i].omega = elements[i].omega_temp + elements[i].E[stage][5];

			if (elements[i].k < 0) elements[i].k = DBL_EPSILON;
			if (elements[i].omega < 0) {
				elements[i].omega = 0.0;
				int n;
				for (int in = 0; in < 3; in++) {
					n = faces[elements[i].face[in]].elementL;
					if (n == i)
						n = faces[elements[i].face[in]].elementR;
					elements[i].omega += elements[n].omega;
				}
				elements[i].omega /= 3.0;
				if (elements[i].omega < 0) elements[i].omega = elements[i].omega_temp;
				if (elements[i].omega < 0) {
					std::cout << global_time << "omega: "<<elements[i].omega << std::endl;
					exit(0);
				}
			}
		}
	}
}

void updateStateRK4(std::vector<Point > &points,std::vector<Element > &elements,std::vector<Face > &faces,std::ofstream &file) {
	double resnorm_laminar,resnorm_turbulent;
	resnorm_laminar = resnorm_turbulent = 0.0;
	//Stage 1
	calLocalTimeStep(elements,faces);	

	computeResidue(elements,faces);
	getE(elements,faces,0);
	update(elements,faces,0);

	//Stage 2
	applyBC(points,elements,faces);
	calFlux(points,elements,faces);
	computeResidue(elements,faces);
	getE(elements,faces,1);
	update(elements,faces,1);

	//Stage 3
	applyBC(points,elements,faces);
	calFlux(points,elements,faces);
	computeResidue(elements,faces);
	getE(elements,faces,2);
	update(elements,faces,2);

	//Final update
	applyBC(points,elements,faces);
	calFlux(points,elements,faces);
	computeResidue(elements,faces);
	getE(elements,faces,3);
	// update(elements,faces,3);   //For modified RK4


	for (int i = 0; i < N_elem; i++) {
		elements[i].rho = elements[i].rho_temp + (1.0/6.0*elements[i].E[0][0] + 1.0/3.0*elements[i].E[1][0] + 1.0/3.0*elements[i].E[2][0] + 1.0/6.0*elements[i].E[3][0]);
		elements[i].u = elements[i].u_temp + (1.0/6.0*elements[i].E[0][1] + 1.0/3.0*elements[i].E[1][1] + 1.0/3.0*elements[i].E[2][1] + 1.0/6.0*elements[i].E[3][1]);
		elements[i].v = elements[i].v_temp + (1.0/6.0*elements[i].E[0][2] + 1.0/3.0*elements[i].E[1][2] + 1.0/3.0*elements[i].E[2][2] + 1.0/6.0*elements[i].E[3][2]);
		elements[i].p = elements[i].p_temp + (1.0/6.0*elements[i].E[0][3] + 1.0/3.0*elements[i].E[1][3] + 1.0/3.0*elements[i].E[2][3] + 1.0/6.0*elements[i].E[3][3]);
		
		if (turbulence != NO_TURB) {
			elements[i].k = elements[i].k_temp + (1.0/6.0*elements[i].E[0][4] + 1.0/3.0*elements[i].E[1][4] + 1.0/3.0*elements[i].E[2][4] + 1.0/6.0*elements[i].E[3][4]);
			elements[i].omega = elements[i].omega_temp + (1.0/6.0*elements[i].E[0][5] + 1.0/3.0*elements[i].E[1][5] + 1.0/3.0*elements[i].E[2][5] + 1.0/6.0*elements[i].E[3][5]);

			if (elements[i].k < 0) elements[i].k = DBL_EPSILON;
			if (elements[i].omega < 0) {
				elements[i].omega = 0.0;
				int n;
				for (int in = 0; in < 3; in++) {
					n = faces[elements[i].face[in]].elementL;
					if (n == i)
						n = faces[elements[i].face[in]].elementR;
					elements[i].omega += elements[n].omega;
				}
				elements[i].omega /= 3.0;
				if (elements[i].omega < 0) elements[i].omega = elements[i].omega_temp;
				if (elements[i].omega < 0) {
					std::cout << global_time << "omega: "<<elements[i].omega << std::endl;
					exit(0);
				}
			}
		}
		computeResNorm(resnorm_laminar,resnorm_turbulent,elements,i);
	}
	file << sqrt(resnorm_laminar)<< " " << sqrt(resnorm_turbulent) << std::endl;	
}

void updateState(std::vector<Point > &points,std::vector<Element > &elements,std::vector<Face > &faces,std::ofstream &file) {
	double delta_t;
	double resnorm_laminar,resnorm_turbulent;
	resnorm_laminar = resnorm_turbulent = 0;
	
	Face face;
	computeResidue(elements,faces);
	calLocalTimeStep(elements,faces);

	for (int i = 0; i < N_elem; i++) {
		delta_t = elements[i].delta_t;
		//UPDATE
		elements[i].rho_temp = elements[i].rho;
		elements[i].u_temp = elements[i].u;
		elements[i].v_temp = elements[i].v;
		elements[i].p_temp = elements[i].p;
		elements[i].k_temp = elements[i].k;
		elements[i].omega_temp = elements[i].omega;

		
		if (other != PI) {
			elements[i].rho += (-1.0*delta_t/elements[i].volume) *elements[i].residue[0];
			elements[i].u += (-1.0*delta_t/elements[i].volume) * (-1.0*elements[i].u_temp*elements[i].residue[0]/elements[i].rho_temp + elements[i].residue[1]/elements[i].rho_temp);
			elements[i].v += (-1.0*delta_t/elements[i].volume) * (-1.0*elements[i].v_temp*elements[i].residue[0]/elements[i].rho_temp + elements[i].residue[2]/elements[i].rho_temp);
			elements[i].p += (-1.0*delta_t/elements[i].volume) * (0.5*(gam-1)*(elements[i].u_temp*elements[i].u_temp + elements[i].v_temp*elements[i].v_temp)*elements[i].residue[0] - (gam-1)*elements[i].u_temp*elements[i].residue[1] - (gam-1)*elements[i].v_temp*elements[i].residue[2] + (gam-1)*elements[i].residue[3]);
			if (turbulence != NO_TURB) {
				elements[i].k += (-1.0*delta_t/elements[i].volume) * (-1.0*elements[i].k_temp*elements[i].residue[0]/elements[i].rho_temp + elements[i].residue[4]/elements[i].rho_temp);
				elements[i].omega += (-1.0*delta_t/elements[i].volume) * (-1.0*elements[i].omega_temp*elements[i].residue[0]/elements[i].rho_temp + elements[i].residue[5]/elements[i].rho_temp);
				//Source terms
				elements[i].k += delta_t*elements[i].source1/elements[i].rho_temp;
				elements[i].omega += delta_t*elements[i].source2/elements[i].rho_temp;
			}
		}
		else {
			elements[i].rho += (-1.0*delta_t/elements[i].volume) * elements[i].residue[0];
			elements[i].u += (-1.0*delta_t/elements[i].volume) * (-1.0*elements[i].u_temp*elements[i].residue[0]/elements[i].rho_temp + elements[i].residue[1]/elements[i].rho_temp);
			elements[i].v += (-1.0*delta_t/elements[i].volume) * (-1.0*elements[i].v_temp*elements[i].residue[0]/elements[i].rho_temp + elements[i].residue[2]/elements[i].rho_temp);
			elements[i].p += (-1.0*delta_t/elements[i].volume) * (0.5*(gam-1)*(elements[i].u_temp*elements[i].u_temp + elements[i].v_temp*elements[i].v_temp)*elements[i].residue[0] - (gam-1)*elements[i].u_temp*elements[i].residue[1] - (gam-1)*elements[i].v_temp*elements[i].residue[2] + (gam-1)*elements[i].residue[3]);	
			if (turbulence != NO_TURB) {
				elements[i].k += (-1.0*delta_t/elements[i].volume) * (-1.0*elements[i].k_temp*elements[i].residue[0]/elements[i].rho_temp + 1.0/(1.0 + betastar*elements[i].omega_temp*delta_t) * (elements[i].residue[4] - elements[i].source1*elements[i].volume)/elements[i].rho_temp);
				elements[i].omega += (-1.0*delta_t/elements[i].volume) * (-1.0*elements[i].omega_temp*elements[i].residue[0]/elements[i].rho_temp + 1.0/(1.0 + 2.0*elements[i].beta*elements[i].omega_temp*delta_t) * (elements[i].residue[5] - elements[i].source2*elements[i].volume)/elements[i].rho_temp);

				if (elements[i].k < 0) elements[i].k = DBL_EPSILON;
				if (elements[i].omega < 0) {
					elements[i].omega = 0.0;
					int n;
					for (int in = 0; in < 3; in++) {
						n = faces[elements[i].face[in]].elementL;
						if (n == i)
							n = faces[elements[i].face[in]].elementR;
						elements[i].omega += elements[n].omega;
					}
					elements[i].omega /= 3.0;
					if (elements[i].omega < 0) elements[i].omega = elements[i].omega_temp;
					if (elements[i].omega < 0) {
						std::cout << global_time << "omega: "<<elements[i].omega << std::endl;
						exit(0);
					}
				}
			}
		
		}

		computeResNorm(resnorm_laminar,resnorm_turbulent,elements,i);

	}
	file << sqrt(resnorm_laminar) << " " << sqrt(resnorm_turbulent) << std::endl;

}

	
#endif


