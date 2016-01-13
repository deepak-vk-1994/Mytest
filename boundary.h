#ifndef _BOUNDARY_H
#define _BOUNDARY_H

#include<math.h>
#include "globals.h"
#include "grid.h"
#include "fluxcal.h"

void applyBC(std::vector<Point > &points, std::vector<Element > &elements, std::vector<Face > &faces) {
	int vertex1, vertex2;
	
	for (std::vector<Face >::iterator it = faces.begin(); it != faces.end(); it++) {
		if ((*it).marker == INTERIOR) continue;
		//INLET and OUTLET
		if (mach_switch == SUBSONIC) {
		
			if ((*it).marker == TOP || (*it).marker == BOTTOM) {
				elements[(*it).elementL].rho = elements[(*it).elementR].rho;
				elements[(*it).elementL].p = elements[(*it).elementR].p;
				#ifdef VISCOUS
					elements[(*it).elementL].u = -elements[(*it).elementR].u;
					elements[(*it).elementL].v = -elements[(*it).elementR].v;
				#else
					elements[(*it).elementL].u = elements[(*it).elementR].u - 2.0 *( elements[(*it).elementR].u*((*it).nx)*((*it).nx) + elements[(*it).elementR].v*((*it).nx)*((*it).ny) );
					elements[(*it).elementL].v = elements[(*it).elementR].v - 2.0 *( elements[(*it).elementR].u*((*it).nx)*((*it).ny) + elements[(*it).elementR].v*((*it).ny)*((*it).ny) );
				#endif
			}
			else if ((*it).marker == OTHER) {
				elements[(*it).elementL].rho = elements[(*it).elementR].rho;
				elements[(*it).elementL].p = elements[(*it).elementR].p;
				#ifdef VISCOUS
					elements[(*it).elementL].u = -elements[(*it).elementR].u;
					elements[(*it).elementL].v = -elements[(*it).elementR].v;
				#else
					elements[(*it).elementL].u = elements[(*it).elementR].u - 2.0 *( elements[(*it).elementR].u*((*it).nx)*((*it).nx) + elements[(*it).elementR].v*((*it).nx)*((*it).ny) );
					elements[(*it).elementL].v = elements[(*it).elementR].v - 2.0 *( elements[(*it).elementR].u*((*it).nx)*((*it).ny) + elements[(*it).elementR].v*((*it).ny)*((*it).ny) );
				#endif
			}

			else if ((*it).marker == LEFT) {
				elements[(*it).elementL].p = elements[(*it).elementR].p;
				elements[(*it).elementL].u = u_inf;
				elements[(*it).elementL].v = v_inf;
				elements[(*it).elementL].w = w_inf;
				elements[(*it).elementL].rho = rho_inf;
			}
			else if ((*it).marker == RIGHT) {
				elements[(*it).elementL].p = p_inf;
				elements[(*it).elementL].u = elements[(*it).elementR].u;
				elements[(*it).elementL].v = elements[(*it).elementR].v;
				elements[(*it).elementL].w = elements[(*it).elementR].w;
				elements[(*it).elementL].rho = elements[(*it).elementR].rho;
			}

		}

		else if (mach_switch == SUPERSONIC) {
			
			if ((*it).marker == TOP || (*it).marker == BOTTOM) {
				elements[(*it).elementL].rho = elements[(*it).elementR].rho;
				elements[(*it).elementL].p = elements[(*it).elementR].p;
				#ifdef VISCOUS
					elements[(*it).elementL].u = -elements[(*it).elementR].u;
					elements[(*it).elementL].v = -elements[(*it).elementR].v;
				#else
					elements[(*it).elementL].u = elements[(*it).elementR].u - 2.0 *( elements[(*it).elementR].u*((*it).nx)*((*it).nx) + elements[(*it).elementR].v*((*it).nx)*((*it).ny) );
					elements[(*it).elementL].v = elements[(*it).elementR].v - 2.0 *( elements[(*it).elementR].u*((*it).nx)*((*it).ny) + elements[(*it).elementR].v*((*it).ny)*((*it).ny) );
				#endif
			}
			else if ((*it).marker == OTHER) {
				elements[(*it).elementL].rho = elements[(*it).elementR].rho;
				elements[(*it).elementL].p = elements[(*it).elementR].p;
				#ifdef VISCOUS
					elements[(*it).elementL].u = -elements[(*it).elementR].u;
					elements[(*it).elementL].v = -elements[(*it).elementR].v;
				#else
					elements[(*it).elementL].u = elements[(*it).elementR].u - 2.0 *( elements[(*it).elementR].u*((*it).nx)*((*it).nx) + elements[(*it).elementR].v*((*it).nx)*((*it).ny) );
					elements[(*it).elementL].v = elements[(*it).elementR].v - 2.0 *( elements[(*it).elementR].u*((*it).nx)*((*it).ny) + elements[(*it).elementR].v*((*it).ny)*((*it).ny) );
				#endif
			}

			else if ((*it).marker == LEFT) {
				elements[(*it).elementL].p = p_inf;
				elements[(*it).elementL].u = u_inf;
				elements[(*it).elementL].v = v_inf;
				elements[(*it).elementL].w = w_inf;
				elements[(*it).elementL].rho = rho_inf;
			}
			else if ((*it).marker == RIGHT) {
				elements[(*it).elementL].p = elements[(*it).elementR].p;
				elements[(*it).elementL].u = elements[(*it).elementR].u;
				elements[(*it).elementL].v = elements[(*it).elementR].v;
				elements[(*it).elementL].w = elements[(*it).elementR].w;
				elements[(*it).elementL].rho = elements[(*it).elementR].rho;
			}

		}
	}

	#ifdef DEBUG
		ofstream file("./Debug/BC.txt");
		for (std::vector<Element >::iterator it = elements.begin(); it != elements.end(); it++) {
			file << (*it).rho <<" "<<(*it).p<<" "<<(*it).u<<" "<<(*it).v<<std::endl;
		}
	#endif
}

void calLocalTimeStep(double &delta_t,std::vector<Element > &elements,std::vector<Face > &faces,int i) {
	 double lambda,a;
	 Face e; 
	 int EL,ER;

	 lambda = 0;
	 for (int k = 0; k < 3; k++) {
	 	e = faces[elements[i].face[k]];
	 	EL = (e).elementL;
	 	ER = (e).elementR;
	 	a = 0.5 * ( sqrt(gam*(elements[EL].p/elements[EL].rho)) + sqrt(gam*(elements[ER].p/elements[ER].rho)) );
	 	lambda += (fabs(elements[i].u*(e.nx) + elements[i].v*(e.ny)) + a)*(e.area);
	 }
	 delta_t = elements[i].volume*2.0*CFL/lambda;
	 // cout << delta_t <<std::endl;

	 // delta_t = 1e-7;
}

void computeResNorm(double &resnorm,std::vector<Element > &elements,int i) {
	resnorm += (( (elements[i].residue[0]/(rho_inf*u_inf))*(elements[i].residue[0]/(rho_inf*u_inf)) + (elements[i].residue[1]/(rho_inf*u_inf*u_inf))*(elements[i].residue[1]/(rho_inf*u_inf*u_inf)) + (elements[i].residue[2]/(rho_inf*u_inf*u_inf))*(elements[i].residue[2]/(rho_inf*u_inf*u_inf)) + (elements[i].residue[3]/(rho_inf*rho_inf*u_inf*u_inf*u_inf))*(elements[i].residue[3]/(rho_inf*rho_inf*u_inf*u_inf*u_inf)) ));
}

void computeResNormSingle(double &resnorm,std::vector<Element > &elements, double residue[4]) {
	resnorm += ( (residue[0]/(rho_inf*u_inf))*(residue[0]/(rho_inf*u_inf)) + (residue[1]/(rho_inf*u_inf*u_inf))*(residue[1]/(rho_inf*u_inf*u_inf)) + (residue[2]/(rho_inf*u_inf*u_inf))*(residue[2]/(rho_inf*u_inf*u_inf)) + (residue[3]/(rho_inf*rho_inf*u_inf*u_inf*u_inf))*(residue[3]/(rho_inf*rho_inf*u_inf*u_inf*u_inf)) );
}

void getE(std::vector<Element > &elements, int index) {
	for (int i = 0; i < N_elem; i++) {
		if (index == 0) {
			elements[i].rho_temp = elements[i].rho;
			elements[i].u_temp = elements[i].u;
			elements[i].v_temp = elements[i].v;
			elements[i].p_temp = elements[i].p;
		}
		elements[i].E[index][0] = elements[i].residue[0];
		elements[i].E[index][1] = (-1.0*elements[i].u*elements[i].residue[0]/elements[i].rho + elements[i].residue[1]/elements[i].rho);
		elements[i].E[index][2] = (-1.0*elements[i].v*elements[i].residue[0]/elements[i].rho + elements[i].residue[2]/elements[i].rho);
		elements[i].E[index][3] = (0.5*(gam-1)*(elements[i].u*elements[i].u + elements[i].v*elements[i].v)*elements[i].residue[0] - (gam-1)*elements[i].u*elements[i].residue[1] - (gam-1)*elements[i].v*elements[i].residue[2] + (gam-1)*elements[i].residue[3]);
	}
}

void computeResidue(std::vector<Element > &elements,std::vector<Face > &faces) {
	for (int i = 0; i < N_elem; i++) {
		for (int k = 0; k < 4; k++) elements[i].residue[k] = 0.0;
		for (int j = 0; j < 3; j++) {
			if (faces[elements[i].face[j]].elementR == i) {
				elements[i].residue[0] += faces[elements[i].face[j]].mass_flux;
				elements[i].residue[1] += faces[elements[i].face[j]].mom_x_flux;
				elements[i].residue[2] += faces[elements[i].face[j]].mom_y_flux;
				elements[i].residue[3] += faces[elements[i].face[j]].energy_flux;
			}
			else if (faces[elements[i].face[j]].elementL == i) {
				elements[i].residue[0] -= faces[elements[i].face[j]].mass_flux;
				elements[i].residue[1] -= faces[elements[i].face[j]].mom_x_flux;
				elements[i].residue[2] -= faces[elements[i].face[j]].mom_y_flux;
				elements[i].residue[3] -= faces[elements[i].face[j]].energy_flux;
			}
			else 
				cout<<"ERROR in face-triangle indexing in UPDATE"<<endl;
		}
	}
}

void update(std::vector<Element > &elements,std::vector<Face > &faces,int stage) {
	double delta_t = 0.0;
	for (int i = 0; i < N_elem; i++) {
		calLocalTimeStep(delta_t,elements,faces,i);	
		if (stage != 2) delta_t *= 0.5;
		elements[i].rho = elements[i].rho_temp + (-1.0*delta_t/elements[i].volume) * elements[i].E[stage][0];
		elements[i].u = elements[i].u_temp + (-1.0*delta_t/elements[i].volume) * elements[i].E[stage][1];
		elements[i].v = elements[i].v_temp + (-1.0*delta_t/elements[i].volume) * elements[i].E[stage][2];
		elements[i].p = elements[i].p_temp + (-1.0*delta_t/elements[i].volume) * elements[i].E[stage][3];
	}
}

#ifdef RK4
void updateState(std::vector<Point > &points,std::vector<Element > &elements,std::vector<Face > &faces,std::ofstream &file) {
	double delta_t,resnorm;
	resnorm = delta_t = 0.0;
	
	//Stage 1
	computeResidue(elements,faces);
	getE(elements,0);

	//Stage 2
	update(elements,faces,0);

	//Stage 3
	applyBC(points,elements,faces);
	calFlux(points,elements,faces);
	computeResidue(elements,faces);
	getE(elements,1);
	update(elements,faces,1);

	//Stage 4
	applyBC(points,elements,faces);
	calFlux(points,elements,faces);
	computeResidue(elements,faces);
	getE(elements,2);
	update(elements,faces,2);

	//Final update
	applyBC(points,elements,faces);
	calFlux(points,elements,faces);
	computeResidue(elements,faces);
	getE(elements,3);

	for (int i = 0; i < N_elem; i++) {
		calLocalTimeStep(delta_t,elements,faces,i);	
		elements[i].rho = elements[i].rho_temp + (-1.0*delta_t/elements[i].volume) * (1.0/6.0*elements[i].E[0][0] + 1.0/3.0*elements[i].E[1][0] + 1.0/3.0*elements[i].E[2][0] + 1.0/6.0*elements[i].E[3][0]);
		elements[i].u = elements[i].u_temp + (-1.0*delta_t/elements[i].volume) * (1.0/6.0*elements[i].E[0][1] + 1.0/3.0*elements[i].E[1][1] + 1.0/3.0*elements[i].E[2][1] + 1.0/6.0*elements[i].E[3][1]);
		elements[i].v = elements[i].v_temp + (-1.0*delta_t/elements[i].volume) * (1.0/6.0*elements[i].E[0][2] + 1.0/3.0*elements[i].E[1][2] + 1.0/3.0*elements[i].E[2][2] + 1.0/6.0*elements[i].E[3][2]);
		elements[i].p = elements[i].p_temp + (-1.0*delta_t/elements[i].volume) * (1.0/6.0*elements[i].E[0][3] + 1.0/3.0*elements[i].E[1][3] + 1.0/3.0*elements[i].E[2][3] + 1.0/6.0*elements[i].E[3][3]);
		computeResNorm(resnorm,elements,i);
	}
	file << sqrt(resnorm) <<std::endl;
	
}

#else
void updateState(std::vector<Point > &points,std::vector<Element > &elements,std::vector<Face > &faces,std::ofstream &file) {
	double delta_t,s,a,b,c;
	double residue[4],resnorm;
	resnorm = 0;
	
	Face face;
	for (int i = 0; i < N_elem; i++) {

		for (int k = 0; k < 4; k++) residue[k] = 0.0;

		for (int j = 0; j < 3; j++) {
			if (faces[elements[i].face[j]].elementR == i) {
				residue[0] += faces[elements[i].face[j]].mass_flux;
				residue[1] += faces[elements[i].face[j]].mom_x_flux;
				residue[2] += faces[elements[i].face[j]].mom_y_flux;
				residue[3] += faces[elements[i].face[j]].energy_flux;
			}
			else if (faces[elements[i].face[j]].elementL == i) {
				residue[0] -= faces[elements[i].face[j]].mass_flux;
				residue[1] -= faces[elements[i].face[j]].mom_x_flux;
				residue[2] -= faces[elements[i].face[j]].mom_y_flux;
				residue[3] -= faces[elements[i].face[j]].energy_flux;
			}
			else 
				cout<<"ERROR in face-triangle indexing in UPDATE"<<endl;
		}

		calLocalTimeStep(delta_t,elements,faces,i);
		
		//UPDATE
		
		elements[i].rho += (-1.0*delta_t/elements[i].volume) * residue[0];
		elements[i].u += (-1.0*delta_t/elements[i].volume) * (-1.0*elements[i].u*residue[0]/elements[i].rho + residue[1]/elements[i].rho);
		elements[i].v += (-1.0*delta_t/elements[i].volume) * (-1.0*elements[i].v*residue[0]/elements[i].rho + residue[2]/elements[i].rho);
		elements[i].p += (-1.0*delta_t/elements[i].volume) * (0.5*(gam-1)*(elements[i].u*elements[i].u + elements[i].v*elements[i].v)*residue[0] - (gam-1)*elements[i].u*residue[1] - (gam-1)*elements[i].v*residue[2] + (gam-1)*residue[3]);

		computeResNormSingle(resnorm,elements,residue);

	}
	file << sqrt(resnorm) <<std::endl;
}
#endif
	
#endif
