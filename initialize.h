#ifndef _INITIALIZE_H
#define _INITIALIZE_H

#include <fstream>
#include <vector>
#include <algorithm>
#include <map>
#include <math.h>
#include <utility>
#include <tuple>
#include <assert.h>
#include <string>
#include <float.h>

using namespace std;

void populateFromSTL(std::vector<Point > &points, std::vector<Element > &elements) {
	ifstream file(gridfile);
	double coord1,coord2,coord3;

	Point pt;
	Element elem;

	std::tuple<double,double,double> nextpoint;
	std::map<std::tuple<double,double,double>,int> pointmap;
	std::map<std::tuple<double,double,double>,int>::iterator it;
	int index_pt = 0;
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
	if (debug == 1) {
		ofstream ofile("./Debug/elements.txt");
		for (int i = 0; i < N_elem; i++) {
			ofile << points[elements[i].vertex[0]].x<<" "<<points[elements[i].vertex[0]].y<<" "<<points[elements[i].vertex[0]].z<<endl;
			ofile << points[elements[i].vertex[1]].x<<" "<<points[elements[i].vertex[1]].y<<" "<<points[elements[i].vertex[1]].z<<endl;
			ofile << points[elements[i].vertex[2]].x<<" "<<points[elements[i].vertex[2]].y<<" "<<points[elements[i].vertex[2]].z<<endl;
		}
	}
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
				
				if( fabs(points[vertex1].x-xleft) < DBL_EPSILON && fabs(points[vertex2].x-xleft) < DBL_EPSILON) 
					F.marker = LEFT;
				else if( fabs(points[vertex1].x-xright) < DBL_EPSILON && fabs(points[vertex2].x-xright) < DBL_EPSILON ) 
					F.marker = RIGHT;
				else if( fabs(points[vertex1].y-ytop) < DBL_EPSILON && fabs(points[vertex2].y-ytop) < DBL_EPSILON) 
					F.marker = TOP;
				// else if ((points[vertex1].y) > 0 && (points[vertex2].y) > 0)
				// 	F.marker = TOP;
				// else if( fabs(points[vertex1].y-ybottom) < DBL_EPSILON && fabs(points[vertex2].y-ybottom) < DBL_EPSILON && (points[vertex1].x < 0.0 || points[vertex2].x < 0.0) )
				// 	F.marker = OTHER;
				else if( fabs(points[vertex1].y-ybottom) < DBL_EPSILON && fabs(points[vertex2].y-ybottom) < DBL_EPSILON)
					F.marker = BOTTOM;
				else 
					F.marker = OTHER;

				faces.push_back(F);
				elements[i].face[j] = faceindex;
				faceindex++;		
			}
		}
	}


	
	
	if (debug == 1) {
		for (int i = 0; i < N_elem; i++) { 
			if (i != faces[elements[i].face[0]].elementR && i != faces[elements[i].face[0]].elementL)
				cout << "ERROR in face-triangle indexing -1";
			if (i != faces[elements[i].face[1]].elementR && i != faces[elements[i].face[1]].elementL)
				cout << "ERROR in face-triangle indexing -2";
			if (i != faces[elements[i].face[2]].elementR && i != faces[elements[i].face[2]].elementL)
				cout << "ERROR in face-triangle indexing -3";

		}

		ofstream file2L("./Debug/facesL.txt");
		ofstream file2R("./Debug/facesR.txt");
		ofstream file2T("./Debug/facesT.txt");
		ofstream file2B("./Debug/facesB.txt");
		ofstream file2O("./Debug/facesO.txt");
		for (std::vector<Face >::iterator it = faces.begin(); it != faces.end(); it++) {
			if ((*it).marker == LEFT)
				file2L << points[(*it).vertex1].x<<" "<< points[(*it).vertex1].y <<" "<<(*it).elementL<<" "<<(*it).elementR<<" "<<(*it).marker<<std::endl;
			if ((*it).marker == RIGHT)
				file2R << points[(*it).vertex1].x<<" "<< points[(*it).vertex1].y <<" "<<(*it).elementL<<" "<<(*it).elementR<<" "<<(*it).marker<<std::endl;
			if ((*it).marker == TOP)
				file2T << points[(*it).vertex1].x<<" "<< points[(*it).vertex1].y <<" "<<(*it).elementL<<" "<<(*it).elementR<<" "<<(*it).marker<<std::endl;
			if ((*it).marker == BOTTOM)
				file2B << points[(*it).vertex1].x<<" "<< points[(*it).vertex1].y <<" "<<(*it).elementL<<" "<<(*it).elementR<<" "<<(*it).marker<<std::endl;
			if ((*it).marker == OTHER)
				file2O << points[(*it).vertex1].x<<" "<< points[(*it).vertex1].y <<" "<<(*it).elementL<<" "<<(*it).elementR<<" "<<(*it).marker<<std::endl;
		}

		
	}
}
	
void initializeGhostCells(std::vector<Point > &points, std::vector<Element > &elements, std::vector<Face > &faces) {
	int ghost_index = elements.size(); //Tracks the index of ghost elements, so that the boundary neighbours can get updated
	int pt_index = points.size();
	int flag = 0;
	int ER,V;
	V = 0;
	double x0,y0,nx,ny;
	for (std::vector<Face >::iterator it = faces.begin(); it != faces.end(); it++) {
		flag = 0;
		if ((*it).marker == INTERIOR) continue;

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
					flag = 1;
					break;
				}
			}
			assert(flag == 1);
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

	if (debug == 1) {
		assert(ghost_index == elements.size());
		ofstream file2("./Debug/edgesmodified.txt");
		for (std::vector<Face >::iterator it = faces.begin(); it != faces.end(); it++) {
			file2 << (*it).vertex1 <<" "<<(*it).vertex2<<" "<<(*it).elementL<<" "<<(*it).elementR<<" "<<(*it).marker<<std::endl;
		}
	}	
}


void initializeGeometeryFaces(std::vector<Point > &points, std::vector<Element > &elements, std::vector<Face > &faces) {
	int vertex1, vertex2;
	double x1,y1,z1,x2,y2,z2;
	for (std::vector<Face >::iterator it = faces.begin(); it != faces.end(); it++) {
		vertex1 = (*it).vertex1;
		vertex2 = (*it).vertex2;

		x1 = points[vertex1].x; y1 = points[vertex1].y; z1 = points[vertex1].z;
		x2 = points[vertex2].x; y2 = points[vertex2].y; z2 = points[vertex2].z;

		(*it).area = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1));
		if (n_facet == 1) {
			(*it).nx = -(y1 - y2)/sqrt((y1-y2)*(y1-y2) + (x1-x2)*(x1-x2));
			(*it).ny = (x1 - x2)/sqrt((y1-y2)*(y1-y2) + (x1-x2)*(x1-x2));
		}
		else {
			(*it).nx = (y1 - y2)/sqrt((y1-y2)*(y1-y2) + (x1-x2)*(x1-x2));
			(*it).ny = -(x1 - x2)/sqrt((y1-y2)*(y1-y2) + (x1-x2)*(x1-x2));
		}

		(*it).xmp = (x1 + x2)/2.0;
		(*it).ymp = (y1 + y2)/2.0;
	}
}

void initializeGeometeryElems(std::vector<Point > &points, std::vector<Element > &elements, std::vector<Face > &faces) {
	int i = 0;
	int vertex1, vertex2, vertex3;
	double x1,y1,z1,x2,y2,z2,x3,y3,z3,a,b,c,s,distance;
	double x1f,y1f,x2f,y2f;
	std::vector<double> distnode1,distnode2,distnode3;
	for (std::vector<Element >::iterator it = elements.begin(); it != elements.end(); it++) {
		i = it - elements.begin();
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

		if (turbulence != NO_TURB) {
			for (std::vector<Face >::iterator it_f = faces.begin(); it_f != faces.end(); it_f++) {
				for (bcns_it = bcns.begin(); bcns_it != bcns.end(); bcns_it++) {
					if ((*it_f).marker == (*bcns_it)) {
						x1f = points[(*it_f).vertex1].x;
						y1f = points[(*it_f).vertex1].y;
						x2f = points[(*it_f).vertex2].x;
						y2f = points[(*it_f).vertex2].y;
						//node1
						distance = sqrt( (x1-x1f)*(x1-x1f) + (y1-y1f)*(y1-y1f) );
						distnode1.push_back(distance);
						distance = sqrt( (x1-x2f)*(x1-x2f) + (y1-y2f)*(y1-y2f) );
						distnode1.push_back(distance);
						distance = sqrt( (x1-(*it_f).xmp)*(x1-(*it_f).xmp) + (y1-(*it_f).ymp)*(y1-(*it_f).ymp) );
						distnode1.push_back(distance);
						//node2
						distance = sqrt( (x2-x1f)*(x2-x1f) + (y2-y1f)*(y2-y1f) );
						distnode2.push_back(distance);
						distance = sqrt( (x2-x2f)*(x2-x2f) + (y2-y2f)*(y2-y2f) );
						distnode2.push_back(distance);
						distance = sqrt( (x2-(*it_f).xmp)*(x2-(*it_f).xmp) + (y2-(*it_f).ymp)*(y2-(*it_f).ymp) );
						distnode2.push_back(distance);
						//node3
						distance = sqrt( (x3-x1f)*(x3-x1f) + (y3-y1f)*(y3-y1f) );
						distnode3.push_back(distance);
						distance = sqrt( (x3-x2f)*(x3-x2f) + (y3-y2f)*(y3-y2f) );
						distnode3.push_back(distance);
						distance = sqrt( (x3-(*it_f).xmp)*(x3-(*it_f).xmp) + (y3-(*it_f).ymp)*(y3-(*it_f).ymp) );
						distnode3.push_back(distance);
					}
				}
			}
			distance = 0.0;
			distance += *(std::min_element(distnode1.begin(),distnode1.end()));
			distance += *(std::min_element(distnode2.begin(),distnode2.end()));
			distance += *(std::min_element(distnode3.begin(),distnode3.end()));
			elements[i].dist = distance/3.0;
			distnode1.clear();
			distnode2.clear();
			distnode3.clear();
		}
	}
	if (debug == 1) {
		double vol = 0;
		for (int i = 0; i < N_elem; i++) 
			vol += elements[i].volume;
		cout << "Volume of the domain: "<< vol << endl;
	}
}			
#endif


