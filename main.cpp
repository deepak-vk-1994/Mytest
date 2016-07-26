#include <iostream>
#include <vector>
#include <fstream>

#include "globals.h"
#include "grid.h"
#include "initialize.h"
#include "fluxcal.h"
#include "boundary.h"
#include "print.h"

using namespace std;

int main() {
	setInputVariables();
	vector<Point> points; 
	vector<Element> elements; 
	vector<Face> faces;

	ofstream fileres,fileLD,fileLDSurf;
	if (ord_accuracy == HO) {
		fileres.open("./Results/ResnormHigh.txt");
		fileLD.open("./Results/LDHigh.txt");
		fileLDSurf.open("./Results/LDSurfHigh.txt");
	}
	else {	
		fileres.open("./Results/Resnorm.txt");
		fileLD.open("./Results/LD.txt");
		fileLDSurf.open("./Results/LDSurf.txt");
	}

	if (restart == 0) {
		populateFromSTL(points,elements);
		storeNeighbours(points,elements,faces);
		initializeGeometeryFaces(points,elements,faces);	
		initializeGhostCells(points,elements,faces);
		initializeGeometeryElems(points,elements,faces);
	}
	else {
		populateGrid(points,elements,faces);
		storeNeighbours(points,elements,faces);
		initializeGeometeryFaces(points,elements,faces);	
		initializeGhostCells(points,elements,faces);
		initializeGeometeryElems(points,elements,faces);
	}


	if (BCT == SLIP || BCT == NO_SLIP)  bcs.push_back(TOP);
	if (BCB == SLIP || BCB == NO_SLIP)  bcs.push_back(BOTTOM);
	if (BCO == SLIP || BCO == NO_SLIP)  bcs.push_back(OTHER);

	if (BCT == NO_SLIP)  bcns.push_back(TOP);
	if (BCB == NO_SLIP)  bcns.push_back(BOTTOM);
	if (BCO == NO_SLIP)  bcns.push_back(OTHER);
	
	for (int t = 0; t < SIMULATION_TIME; t++) {	
		global_time = t;
		applyBC(points,elements,faces);
		calFlux(points,elements,faces);
		if (time_int == RK4)
			updateStateRK4(points,elements,faces,fileres);
		else
			updateState(points,elements,faces,fileres);
		if (t%1000 == 0) {
			printFlow(points,elements,t);
			// printBLData(points,elements,faces,t);
//			printCP(points,elements,faces,t);
//			extractFromLine(points,elements,faces);
//			printYPlus(points,elements,faces,t);
		}

//		if (t%1000 == 0) {
//			calLiftAndDrag(points,elements,faces,fileLD);
//			calLiftAndDragUsingSurf(points,elements,faces,fileLDSurf);
//		}
	}
}
