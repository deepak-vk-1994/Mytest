#include<iostream>
#include<vector>
#include<fstream>

#include "globals.h"
#include "grid.h"
#include "initialize.h"
#include "fluxcal.h"
#include "boundary.h"
#include "print.h"

using namespace std;

int main() {
	setInputVariables();
	vector<Point> points;  points.reserve(N_pts);
	vector<Element> elements; elements.reserve(N_elem);
	vector<Face> faces;
	cout <<p_inf<<" "<<rho_inf<<" "<<u_inf<<endl;
	#ifdef SECONDORDER
	ofstream fileres("./Results/ResnormHigh.txt");
	#else	
	ofstream fileres("./Results/Resnorm.txt");
	#endif

	populateFromSTL(points,elements);
	storeNeighbours(points,elements,faces);
	initializeGeometery(points,elements,faces);	
	initializeGhostCells(points,elements,faces);
	
	for (int t = 0; t < SIMULATION_TIME; t++) {
		applyBC(points,elements,faces);
		calFlux(points,elements,faces);
		updateState(points,elements,faces,fileres);

		if (t%1000 == 0) {
			printFlow(points,elements,t);
			printSurfaceP(points,elements,faces,t);
		}
	}
}
