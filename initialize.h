#ifndef _INITIALIZE_H
#define _INITIALIZE_H

#include<fstream>
#include<vector>
#include<map>
#include<math.h>
#include<utility>
#include<tuple>
#include<assert.h>
#include<string>
#include "globals.h"
#include "grid.h"

using namespace std;

void setInputVariables() {
	ifstream file("./Input.txt");
	string line;
	int i,j;
	i = j = 0;
	double variables[11],p0,rho0,M;
	
	while(file) {
		std::getline(file,line);
		if (i%2 == 0) {
			i++;
			continue;
		}
		variables[j] = stod(line);
		i++;
		j++;
	}
	p0 = variables[0];
	M = variables[1];
	rho0 = variables[2];
	CFL = variables[3];
	T_ref = variables[4];
	mu_ref = variables[5];
	Pr = variables[6];
	N_elem = (int)variables[7];
	N_pts = (int)variables[8];
	SIMULATION_TIME = (int)variables[9];
	mach_switch = (int)variables[10];

	p_inf = 1013.25; //p0/(pow((1 + (gam-1)*0.5*M*M),(gam/(gam-1))));
	rho_inf = 0.011699;//rho0/(pow((1 + (gam-1)*0.5*M*M),(1/(gam-1))));
	u_inf = 139.28;//M*sqrt(gam*p_inf/rho_inf);
	
	double T_inf = p_inf/(rho_inf*R);
	double visc_inf = ( mu_ref * pow((T_inf/T_ref),1.5) * (T_ref + S1)/(T_inf + S1) );
	// double visc_inf = mu_ref;
	Re = rho_inf*u_inf*0.18/visc_inf;
	cout << "Reynold's number: "<<Re<<"  Vel: "<<u_inf<<"   Rho: "<<rho_inf<<"D"<<D<<"  p_inf: "<<p_inf<<" T_inf: "<< T_inf <<"Visc:"<<visc_inf<< std::endl;
	file.close();
}


void populateFromSTL(std::vector<Point > &points, std::vector<Element > &elements) {
	ifstream file("./Grids/check.txt");
	double coord1,coord2,coord3;

	Point pt;
	Element elem;

	std::tuple<double,double,double> nextpoint;
	std::map<std::tuple<double,double,double>,int> pointmap;
	std::map<std::tuple<double,double,double>,int>::iterator it;
	int index_pt = 0;
	int index_elem = 0;
	int index = 0;
	int i = 0;
	int line = 0;
	while(file >> coord1 >> coord2 >> coord3) {
		line++;
		nextpoint = std::make_tuple(coord1,coord2,coord3); //Point read
		it = pointmap.find(nextpoint);
		if (it != pointmap.end()) {	//Point already indexed -- just add to triangle
			index = (*it).second;
			elem.vertex[i] = index;			
			i++;
			if (i%3 == 0) i = 0;	
		}
		else {   //New point -- index it
			pt.x = coord1; 
			pt.y = coord2;
			pt.z = coord3; 
			points.push_back(pt);
			elem.vertex[i] = index_pt;			
			i++;
			if (i%3 == 0) i = 0;

			pointmap[nextpoint] = index_pt;
			index_pt++;
		}
		if (line%3 == 0) elements.push_back(elem); //Every three points is a triangle
	}

	N_elem = elements.size();
	N_pts = points.size();


	#ifdef DEBUG
		ofstream ofile("./Debug/elements.txt");
		for (int i = 0; i < N_elem; i++) {
			ofile << points[elements[i].vertex[0]].x<<" "<<points[elements[i].vertex[0]].y<<" "<<points[elements[i].vertex[0]].z<<endl;
			ofile << points[elements[i].vertex[1]].x<<" "<<points[elements[i].vertex[1]].y<<" "<<points[elements[i].vertex[1]].z<<endl;
			ofile << points[elements[i].vertex[2]].x<<" "<<points[elements[i].vertex[2]].y<<" "<<points[elements[i].vertex[2]].z<<endl;
		}
	#endif
	file.close();
//	ofile.close();
}

void populateFromNastran(std::vector<Point > &points, std::vector<Element > &elements) {
	//TODO
}

void storeNeighbours(std::vector<Point > &points, std::vector<Element > &elements, std::vector<Face > &faces) {
	std::pair<int,int> face,edge_reverse,index_elem;
	//face pair, reverse of the face pair, and index_elem keeps track of the triangle index and the vertex index

	std::map<std::pair<double,double>,std::pair<int,int>> neighbourmap;
	std::map<std::pair<double,double>,std::pair<int,int>>::iterator it;

	Face F;
	int faceindex = 0;
	
	for (int i = 0; i < N_elem; i++) {  //i is the triangle number
		//EDGE 1
		face.first = elements[i].vertex[0];
		face.second = elements[i].vertex[1];
		edge_reverse.first = face.second;
		edge_reverse.second = face.first;

		index_elem.first = i;
		index_elem.second = 0;

		neighbourmap[face] = index_elem; 
		
		it = neighbourmap.find(edge_reverse);
		if (it != neighbourmap.end()) {
			//neighbour of the reverse face is the current triangle
			F.vertex1 = face.first;
			F.vertex2 = face.second;
			F.elementL = (*it).second.first;
			F.elementR = i;
			faces.push_back(F);	
			elements[i].face[0] = faceindex;
			elements[(*it).second.first].face[(*it).second.second] = faceindex;
			faceindex++;
		}
	
		//EDGE 2
		face.first = elements[i].vertex[1];
		face.second = elements[i].vertex[2];
		edge_reverse.first = face.second;
		edge_reverse.second = face.first;

		index_elem.first = i;
		index_elem.second = 1;

		neighbourmap[face] = index_elem; 
		
		it = neighbourmap.find(edge_reverse);
		if (it != neighbourmap.end()) {
			F.vertex1 = face.first;
			F.vertex2 = face.second;
			F.elementL = (*it).second.first;
			F.elementR = i;
			faces.push_back(F);	
			elements[i].face[1] = faceindex;
			elements[(*it).second.first].face[(*it).second.second] = faceindex;
			faceindex++;
		}

		//EDGE 3
		face.first = elements[i].vertex[2];
		face.second = elements[i].vertex[0];
		edge_reverse.first = face.second;
		edge_reverse.second = face.first;

		index_elem.first = i;
		index_elem.second = 2;

		neighbourmap[face] = index_elem; 
		
		it = neighbourmap.find(edge_reverse);
		if (it != neighbourmap.end()) {
			F.vertex1 = face.first;
			F.vertex2 = face.second;
			F.elementL = (*it).second.first;
			F.elementR = i;
			faces.push_back(F);	
			elements[i].face[2] = faceindex;
			elements[(*it).second.first].face[(*it).second.second] = faceindex;
			faceindex++;
		}
	}

	int vertex1, vertex2;

	for (int i = 0; i < N_elem; i++) { 
		for (int j = 0; j < 3; j++) {
			if (elements[i].face[j] == -1) {
				vertex1 = elements[i].vertex[j];
				vertex2 = elements[i].vertex[(j == 2 ? 0 : j+1)];

				F.vertex1 = vertex1;
				F.vertex2 = vertex2;
				F.elementL = -1;
				F.elementR = i;
				
				if(points[vertex1].x == xleft && points[vertex2].x == xleft) 
					F.marker = LEFT;
				else if(points[vertex1].x == xright && points[vertex2].x == xright) 
					F.marker = RIGHT;
				else if(points[vertex1].y == ytop && points[vertex2].y == ytop) 
					F.marker = TOP;
//				else if(points[vertex1].y == ybottom && points[vertex2].y == ybottom)
//					F.marker = BOTTOM;
				else 
					F.marker = BOTTOM;

				faces.push_back(F);
				elements[i].face[j] = faceindex;
				faceindex++;		
			}
		}
	}


	
	
	#ifdef DEBUG
		for (int i = 0; i < N_elem; i++) { 
			if (i != faces[elements[i].face[0]].elementR && i != faces[elements[i].face[0]].elementL)
				cout << "ERROR in face-triangle indexing -1";
			if (i != faces[elements[i].face[1]].elementR && i != faces[elements[i].face[1]].elementL)
				cout << "ERROR in face-triangle indexing -2";
			if (i != faces[elements[i].face[2]].elementR && i != faces[elements[i].face[2]].elementL)
				cout << "ERROR in face-triangle indexing -3";

		}

		ofstream file2("./Debug/faces.txt");
		for (std::vector<Face >::iterator it = faces.begin(); it != faces.end(); it++) {
//			file2 << (*it).vertex1 <<" "<<(*it).vertex2<<" "<<(*it).elementL<<" "<<(*it).elementR<<" "<<(*it).marker<<std::endl;
		//	if ((*it).marker == LEFT || (*it).marker == RIGHT)
		//	if ((*it).marker == TOP || (*it).marker == BOTTOM || (*it).marker == LEFT || (*it).marker == RIGHT)
			if ((*it).marker == RIGHT)
				file2 << points[(*it).vertex1].x<<" "<< points[(*it).vertex1].y <<" "<<(*it).elementL<<" "<<(*it).elementR<<" "<<(*it).marker<<std::endl;
		}

		
	#endif
}
	
void initializeGhostCells(std::vector<Point > &points, std::vector<Element > &elements, std::vector<Face > &faces) {
//	Element elem;
	int ghost_index = elements.size(); //Tracks the index of ghost elements, so that the boundary neighbours can get updated
	int pt_index = points.size();

	int i = 0;
	int j;
	int ER,V;
	double x0,y0,z0,nx,ny;
	for (std::vector<Face >::iterator it = faces.begin(); it != faces.end(); it++) {
		if ((*it).marker == INTERIOR) continue;

		i = (*it).elementR;


		if ((*it).marker == LEFT || (*it).marker == RIGHT || (*it).marker == TOP || (*it).marker == BOTTOM || (*it).marker == OTHER) {
			Element elem;
			Point pt;

			(*it).elementL = ghost_index;
			elem.vertex[0] = (*it).vertex2;
			elem.vertex[1] = (*it).vertex1;
			ER = (*it).elementR;
			for (int k = 0; k < 3; k++) {
				if (elements[ER].vertex[k] != (*it).vertex1 && elements[ER].vertex[k] != (*it).vertex2) {
					V = k;
					break;
				}
			}
			x0 = points[elements[ER].vertex[V]].x - points[(*it).vertex1].x;
			y0 = points[elements[ER].vertex[V]].y - points[(*it).vertex1].y;
			nx = (*it).nx; ny = (*it).ny;
			pt.x = points[(*it).vertex1].x + x0 - 2.0*(x0*nx + y0*ny)*nx;
			pt.y = points[(*it).vertex1].y + y0 - 2.0*(x0*nx + y0*ny)*ny;
			pt.z = 0.0;
			points.push_back(pt);
			elem.vertex[2] = pt_index;
			elements.push_back(elem);

			pt_index++;
			ghost_index++;
		}
		
	}

	#ifdef DEBUG
		assert(ghost_index == elements.size());
		ofstream file2("./Debug/edgesmodified.txt");
		for (std::vector<Face >::iterator it = faces.begin(); it != faces.end(); it++) {
			file2 << (*it).vertex1 <<" "<<(*it).vertex2<<" "<<(*it).elementL<<" "<<(*it).elementR<<" "<<(*it).marker<<std::endl;
		}
	#endif			
}


void initializeGeometery(std::vector<Point > &points, std::vector<Element > &elements, std::vector<Face > &faces) {
	int vertex1, vertex2, vertex3;
	double x1,y1,z1,x2,y2,z2,x3,y3,z3,a,b,c,s;
	for (std::vector<Face >::iterator it = faces.begin(); it != faces.end(); it++) {
		vertex1 = (*it).vertex1;
		vertex2 = (*it).vertex2;

		x1 = points[vertex1].x; y1 = points[vertex1].y; z1 = points[vertex1].z;
		x2 = points[vertex2].x; y2 = points[vertex2].y; z2 = points[vertex2].z;

		(*it).area = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1));
		//facet normal : 0 0 1
		(*it).nx = -(y1 - y2)/sqrt((y1-y2)*(y1-y2) + (x1-x2)*(x1-x2));
		(*it).ny = (x1 - x2)/sqrt((y1-y2)*(y1-y2) + (x1-x2)*(x1-x2));
		//facet normal : 0 0 -1
		// (*it).nx = (y1 - y2)/sqrt((y1-y2)*(y1-y2) + (x1-x2)*(x1-x2));
		// (*it).ny = -(x1 - x2)/sqrt((y1-y2)*(y1-y2) + (x1-x2)*(x1-x2));

		(*it).xmp = (x1 + x2)/2.0;
		(*it).ymp = (y1 + y2)/2.0;
	}


	for (int i = 0; i < N_elem; i++) {
		vertex1 = elements[i].vertex[0]; vertex2 = elements[i].vertex[1]; vertex3 = elements[i].vertex[2];
		x1 = points[vertex1].x; y1 = points[vertex1].y; z1 = points[vertex1].z;
		x2 = points[vertex2].x; y2 = points[vertex2].y; z2 = points[vertex2].z;
		x3 = points[vertex3].x; y3 = points[vertex3].y; z3 = points[vertex3].z;
		a = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1));
		b = sqrt((x2-x3)*(x2-x3) + (y2-y3)*(y2-y3) + (z2-z3)*(z2-z3));
		c = sqrt((x3-x1)*(x3-x1) + (y3-y1)*(y3-y1) + (z3-z1)*(z3-z1));

		s = 0.5*(a+b+c);
		elements[i].volume = sqrt(s*(s-a)*(s-b)*(s-c));

		elements[i].xc = (x1 + x2 + x3)/3.0;
		elements[i].yc = (y1 + y2 + y3)/3.0;

	}
	#ifdef DEBUG
		double vol = 0;
		for (int i = 0; i < N_elem; i++) 
			vol += elements[i].volume;
		cout << "Volume of the domain: "<< vol << endl;

		ofstream file("./Debug/Geometery.txt");
		for (std::vector<Face >::iterator it = faces.begin(); it != faces.end(); it++) {
			file << (*it).nx<<" "<<(*it).ny<<" "<<(*it).area<<std::endl;
		}
		file.close();
		
	#endif
}			
#endif


