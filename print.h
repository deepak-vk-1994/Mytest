#ifndef _PRINT_H
#define _PRINT_H

#include<math.h>
#include "globals.h"
#include "grid.h"


void printFlow(std::vector<Point > &points, std::vector<Element > &elements, int t) {

	char name[50];
	ofstream file;
	#ifdef SECONDORDER
	sprintf(name,"./Results/BGhigh-%d.vtk",t);
	#else	
	sprintf(name,"./Results/BG-%d.vtk",t);
	#endif

	file.open(name);
	file<<"# vtk DataFile Version 3.0"<<std::endl;
	file<<"Unstructured Grid for Triangles"<<std::endl;
	file<<"ASCII"<<std::endl;
	file<<"DATASET UNSTRUCTURED_GRID"<<std::endl;
	file<<"POINTS "<<N_pts<<" DOUBLE"<<std::endl;
	for (int i = 0; i < N_pts; i++) {
		file<<points[i].x<<" "<<points[i].y<<" "<<points[i].z<<std::endl;
	}
	file<<std::endl;
	file<<"CELLS "<<N_elem<<" "<<N_elem*3.0 + N_elem<<std::endl;
	for (int i = 0; i < N_elem; i++) {
		file <<"3 "<<elements[i].vertex[0]<<" "<<elements[i].vertex[1]<<" "<<elements[i].vertex[2]<<std::endl;
	}
	file<<std::endl;

	file<<"CELL_TYPES "<<N_elem<<std::endl;
	for (int i = 0; i < N_elem; i++) {
		file <<"5"<<std::endl;
	}
	file<<std::endl;
	//DATA
	file<<"CELL_DATA "<<N_elem<<std::endl;
/*	file<<"SCALARS pressure DOUBLE"<<std::endl;
	file<<"LOOKUP_TABLE default"<<std::endl;
	for(int i = 0; i < N_elem; i++) {
		file << elements[i].p<<std::endl;
////		file << sqrt(elements[i].u*elements[i].u + elements[i].v*elements[i].v)/sqrt(gam*R*elements[i].calTemp())<<std::endl;
	}
	file<<std::endl;
	file<<"SCALARS density DOUBLE"<<std::endl;
	file<<"LOOKUP_TABLE default"<<std::endl;
	for(int i = 0; i < N_elem; i++) {
		file << elements[i].rho<<std::endl;
	}
	file<<std::endl;
*/
	file<<"VECTORS velocity DOUBLE"<<std::endl;
	for(int i = 0; i < N_elem; i++) {
		file << elements[i].u<<" "<<elements[i].v<<" "<<elements[i].w<<std::endl;
	}

	file.close();
}

void printSurfaceP(std::vector<Point > &points, std::vector<Element > &elements, std::vector<Face > &faces, int t) {

	char name[50];
	ofstream file;

	#ifdef SECONDORDER
	sprintf(name,"./Results/SPHigh%d.txt",t);
	#else	
	sprintf(name,"./Results/SP%d.txt",t);
	#endif

	file.open(name);
	int elem;
	std::map<double,double> surfp;
	std::map<double,double>::iterator it;
	for (std::vector<Face >::iterator it = faces.begin(); it != faces.end(); it++) {
		if ((*it).marker == BOTTOM) {
			elem = (*it).elementR;
			surfp[elements[elem].xc] = elements[elem].p;
		}
	}
	for (it = surfp.begin(); it != surfp.end(); it++) {
		file << (*it).first << " " << (*it).second << std::endl;
	}
}		


#endif
