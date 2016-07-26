#ifndef _PRINT_H
#define _PRINT_H

#include<math.h>


void printFlow(std::vector<Point > &points, std::vector<Element > &elements, int t) {

	char name[250];
	double machno;
	ofstream file;
	double printelem,printpts;

	if (ghostoutput == 1) {
		printelem = elements.size();
		printpts = points.size();
	}
	else {
		printelem = N_elem;
		printpts = N_pts;
	}

	if (ord_accuracy == HO) 
		sprintf(name,"./Results/BGhigh-%d.vtk",t);
	else	
		sprintf(name,"./Results/BG-%d.vtk",t);

	file.open(name);
	file<<"# vtk DataFile Version 3.0"<<std::endl;
	file<<"Unstructured Grid for Triangles"<<std::endl;
	file<<"ASCII"<<std::endl;
	file<<"DATASET UNSTRUCTURED_GRID"<<std::endl;
	file<<"POINTS "<<printpts<<" DOUBLE"<<std::endl;
	for (int i = 0; i < printpts; i++) {
		file<<points[i].x<<" "<<points[i].y<<" "<<points[i].z<<std::endl;
	}
	file<<std::endl;
	file<<"CELLS "<<printelem<<" "<<printelem*3.0 + printelem<<std::endl;
	for (int i = 0; i < printelem; i++) {
		file <<"3 "<<elements[i].vertex[0]<<" "<<elements[i].vertex[1]<<" "<<elements[i].vertex[2]<<std::endl;
	}
	file<<std::endl;
	file<<"CELL_TYPES "<<printelem<<std::endl;
	for (int i = 0; i < printelem; i++) {
		file <<"5"<<std::endl;
	}
	file<<std::endl;
	//DATA
	file<<"CELL_DATA "<<printelem<<std::endl;
	file<<"SCALARS mut/mu DOUBLE"<<std::endl;
	file<<"LOOKUP_TABLE default"<<std::endl;
	for(int i = 0; i < printelem; i++) {
		file << elements[i].mu_t/elements[i].calViscosity()<<std::endl;
	}
	file<<std::endl;
	file<<"SCALARS k DOUBLE"<<std::endl;
	file<<"LOOKUP_TABLE default"<<std::endl;
	for(int i = 0; i < printelem; i++) {
		file << elements[i].k<<std::endl;
	}
	file<<std::endl;
	file<<"SCALARS omega DOUBLE"<<std::endl;
	file<<"LOOKUP_TABLE default"<<std::endl;
	for(int i = 0; i < printelem; i++) {
		file << elements[i].omega<<std::endl;
	}
	file<<std::endl;
	file<<"SCALARS Ydist DOUBLE"<<std::endl;
	file<<"LOOKUP_TABLE default"<<std::endl;
	for(int i = 0; i < printelem; i++) {
		file << elements[i].dist<<std::endl;
	}
	file<<std::endl;
	file<<"SCALARS F1 DOUBLE"<<std::endl;
	file<<"LOOKUP_TABLE default"<<std::endl;
	for(int i = 0; i < printelem; i++) {
		file << elements[i].F1<<std::endl;
	}
	file<<std::endl;
	file<<"SCALARS vorticity DOUBLE"<<std::endl;
	file<<"LOOKUP_TABLE default"<<std::endl;
	for(int i = 0; i < printelem; i++) {
		file << elements[i].calVorticity()<<std::endl;
	}
	file<<std::endl;
	file<<"SCALARS pressure DOUBLE"<<std::endl;
	file<<"LOOKUP_TABLE default"<<std::endl;
	for(int i = 0; i < printelem; i++) {
		file << elements[i].p<<std::endl;
	}
	file<<std::endl;
	file<<"SCALARS density DOUBLE"<<std::endl;
	file<<"LOOKUP_TABLE default"<<std::endl;
	for(int i = 0; i < printelem; i++) {
		file << elements[i].rho<<std::endl;
	}
	file<<std::endl;
	file<<"SCALARS Machno DOUBLE"<<std::endl;
	file<<"LOOKUP_TABLE default"<<std::endl;
	for(int i = 0; i < printelem; i++) {
		file << sqrt(elements[i].u*elements[i].u + elements[i].v*elements[i].v)/sqrt(gam*R*elements[i].calTemp())<<std::endl;
	}
	file<<std::endl;
	file<<"SCALARS TotalP DOUBLE"<<std::endl;
	file<<"LOOKUP_TABLE default"<<std::endl;
	for(int i = 0; i < printelem; i++) {
		machno = sqrt(elements[i].u*elements[i].u + elements[i].v*elements[i].v)/sqrt(gam*R*elements[i].calTemp());
		file << elements[i].p * (pow((1 + (gam-1)*0.5*machno*machno),(gam/(gam-1)))) <<std::endl;
	}
	file<<std::endl;
	file<<"SCALARS check DOUBLE"<<std::endl;
	file<<"LOOKUP_TABLE default"<<std::endl;
	for(int i = 0; i < printelem; i++) {
		file << elements[i].check<<std::endl;
	}
	file<<std::endl;
	file<<"VECTORS velocity DOUBLE"<<std::endl;
	for(int i = 0; i < printelem; i++) {
		//if (elements[i].v < 1e-15) file << elements[i].u<<" "<<"0.0"<<" "<<elements[i].w<<std::endl;
		file << elements[i].u<<" "<<elements[i].v<<" "<<elements[i].w<<std::endl;
	}
	file<<std::endl;
	file<<"VECTORS gradp DOUBLE"<<std::endl;
	for(int i = 0; i < printelem; i++) {
		file << elements[i].gradp[0]<<" "<<elements[i].gradp[1]<<" "<<"0.0"<<std::endl;
	}
	file<<std::endl;
	file<<"VECTORS gradentropy DOUBLE"<<std::endl;
	for(int i = 0; i < printelem; i++) {
		file << elements[i].gradentropy[0]<<" "<<elements[i].gradentropy[1]<<" "<<"0.0"<<std::endl;
	}
	file.close();
}

void printBLData(std::vector<Point > &points, std::vector<Element > &elements, std::vector<Face > &faces, int t) {

	char nameu[50], nameCF[50], namev[50];
	ofstream fileu, fileCF, filev;
	double tw; //wall shear stress

	if (ord_accuracy == HO) {
		sprintf(nameu,"./Results/UHigh%d.txt",t);
		sprintf(namev,"./Results/VHigh%d.txt",t);
		sprintf(nameCF,"./Results/CFHigh%d.txt",t);
	}
	else {	
		sprintf(nameu,"./Results/U%d.txt",t);
		sprintf(namev,"./Results/V%d.txt",t);
		sprintf(nameCF,"./Results/CF%d.txt",t);
	}

	fileu.open(nameu);
	filev.open(namev);
	fileCF.open(nameCF);
	int elem;
	std::map<double,double> Cf;
	std::map<double,double> BLu;
	std::map<double,double> BLv;
	std::map<double,double>::iterator it;
	for (std::vector<Face >::iterator it = faces.begin(); it != faces.end(); it++) {
		if ((*it).marker == BOTTOM) {
			elem = (*it).elementR;
			tw = elements[elem].calViscosity()*elements[elem].gradu[1];
			Cf[elements[elem].xc] = tw/(0.5*elements[elem].rho*u_inf*u_inf);
		}
	}
	for (it = Cf.begin(); it != Cf.end(); it++) {
		fileCF << (*it).first << " " << (*it).second << std::endl;
	}

	for (int i = 0; i < N_elem; i++) {
		if (elements[i].xc > 0.52 && elements[i].xc < 0.54) {
			if (elements[i].yc < 0.0156)
				BLu[elements[i].yc] = elements[i].u;
		}
	}
	for (it = BLu.begin(); it != BLu.end(); it++) {
		fileu << (*it).first << " " << (*it).second << std::endl;
	}

	for (int i = 0; i < N_elem; i++) {
		if (elements[i].xc > 0.52 && elements[i].xc < 0.54) {
			if (elements[i].yc < 0.0156)
				BLv[elements[i].yc] = elements[i].v;
		}
	}
	for (it = BLv.begin(); it != BLv.end(); it++) {
		filev << (*it).first << " " << (*it).second << std::endl;
	}

}		

void printCP(std::vector<Point > &points, std::vector<Element > &elements, std::vector<Face > &faces, int t) {

	char nameu[250],namel[250];
	ofstream fileu,filel;

	if (ord_accuracy == HO) {
		sprintf(nameu,"./Results/CPupperHigh-%d.txt",t);
		sprintf(namel,"./Results/CPlowerHigh-%d.txt",t);
	}
	else {
		sprintf(nameu,"./Results/CPupper%d.txt",t);
		sprintf(namel,"./Results/CPlower%d.txt",t);
	}

	fileu.open(nameu);
	filel.open(namel);

	int elem;
	std::map<double,double> CPu;
	std::map<double,double> CPl;
	
	std::map<double,double>::iterator it;
	for (std::vector<Face >::iterator it = faces.begin(); it != faces.end(); it++) {
		if ((*it).marker == OTHER) {
			elem = (*it).elementR;
			if ((*it).ny > 0.0)
				CPl[elements[elem].xc] = (elements[elem].p - p_inf)/(0.5*rho_inf*u_inf*u_inf);
			else
				CPu[elements[elem].xc] = (elements[elem].p - p_inf)/(0.5*rho_inf*u_inf*u_inf);
		}
	}
	for (it = CPu.begin(); it != CPu.end(); it++) {
		fileu << (*it).first/chord << " " << (*it).second << std::endl;
	}
	for (it = CPl.begin(); it != CPl.end(); it++) {
		filel << (*it).first/chord << " " << (*it).second << std::endl;
	}
}	

void printYPlus(std::vector<Point > &points, std::vector<Element > &elements, std::vector<Face > &faces, int t) {
	char nameu[250],namel[250];
	ofstream fileu,filel;

	if (ord_accuracy == HO) {
		sprintf(nameu,"./Results/YPupperHigh-%d.txt",t);
		sprintf(namel,"./Results/YPlowerHigh-%d.txt",t);
	}
	else {
		sprintf(nameu,"./Results/YPupper%d.txt",t);
		sprintf(namel,"./Results/YPlower%d.txt",t);
	}

	fileu.open(nameu);
	filel.open(namel);

	std::map<double,double> YPu;
	std::map<double,double> YPl;
	int elem;
	double ustar,tw;
	std::map<double,double>::iterator it;
	for (std::vector<Face >::iterator it = faces.begin(); it != faces.end(); it++) {
		if ((*it).marker == OTHER) {
			elem = (*it).elementR;
			tw = elements[elem].calViscosity()*elements[elem].gradu[1];
			ustar = sqrt(fabs(tw/elements[elem].rho));
			if ((*it).ny > 0.0)
				YPl[elements[elem].xc] = ustar*elements[elem].dist*elements[elem].rho/elements[elem].calViscosity();
			else
				YPu[elements[elem].xc] = ustar*elements[elem].dist*elements[elem].rho/elements[elem].calViscosity();
		}
	}
	for (it = YPu.begin(); it != YPu.end(); it++) {
		fileu << (*it).first/chord << " " << (*it).second << std::endl;
	}
	for (it = YPl.begin(); it != YPl.end(); it++) {
		filel << (*it).first/chord << " " << (*it).second << std::endl;
	}

}

void calLiftAndDrag(std::vector<Point > &points, std::vector<Element > &elements, std::vector<Face > &faces, std::ofstream &file) {	
	double length,lift,drag,x1,x2,y1,y2,dot,elem;
	lift = 0.0;
	drag = 0.0;
	for (std::vector<Face >::iterator it = faces.begin(); it != faces.end(); it++) {
		if ((*it).marker == TOP || (*it).marker == BOTTOM || (*it).marker == LEFT || (*it).marker == RIGHT) {
			elem = (*it).elementR;
			x1 = points[(*it).vertex1].x;
			x2 = points[(*it).vertex2].x;
			y1 = points[(*it).vertex1].y;
			y2 = points[(*it).vertex2].y;
			length = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) );
			dot = (elements[elem].u*(*it).nx + elements[elem].v*(*it).ny)*length;
			lift += -elements[elem].rho*dot*elements[elem].v - elements[elem].p*(*it).ny*length;
			drag += -elements[elem].rho*dot*elements[elem].u - elements[elem].p*(*it).nx*length;
		}
	}
	file << lift <<" "<< drag <<" "<<lift/(0.5*rho_inf*u_inf*u_inf*chord)<<" "<<drag/(0.5*rho_inf*u_inf*u_inf*chord)<< std::endl;
}

void calLiftAndDragUsingSurf(std::vector<Point > &points, std::vector<Element > &elements, std::vector<Face > &faces, std::ofstream &file) {	
	double length,lift,drag,x1,x2,y1,y2,elem,tau;
	lift = 0.0;
	drag = 0.0;
	for (std::vector<Face >::iterator it = faces.begin(); it != faces.end(); it++) {
		if ((*it).marker == OTHER) {
			elem = (*it).elementR;
			x1 = points[(*it).vertex1].x;
			x2 = points[(*it).vertex2].x;
			y1 = points[(*it).vertex1].y;
			y2 = points[(*it).vertex2].y;
			length = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) );
			if ((*it).ny < 0.0) {
				tau = fabs(elements[elem].calViscosity()*elements[elem].gradu[1]);
				lift += -1.0*(-elements[elem].p*length*(*it).ny - tau*length*(*it).nx);
				drag += -1.0*(-elements[elem].p*length*(*it).nx + tau*length*(*it).ny);
			}
			else {
				tau = fabs(elements[elem].calViscosity()*elements[elem].gradu[1]);
				lift += -1.0*(-elements[elem].p*length*(*it).ny + tau*length*(*it).nx);
				drag += -1.0*(-elements[elem].p*length*(*it).nx - tau*length*(*it).ny);
			}
			
		}
	}
	file << lift <<" "<< drag <<" "<<lift/(0.5*rho_inf*u_inf*u_inf*chord)<<" "<<drag/(0.5*rho_inf*u_inf*u_inf*chord)<< std::endl;
}

void extractFromLine(std::vector<Point > &points, std::vector<Element > &elements, std::vector<Face > &faces) {	
	ifstream ifile("./Grids/mid.txt");
	ofstream ofile("./Results/lineextractP.txt");
	ofstream ofile2("./Results/lineextractU.txt");

	double coord1,coord2,x1,y1,x2,y2,xcheck,ycheck,value;
	int flag;
	std::vector<std::pair<double,double>> line;
	std::pair<double,double> coord;
	std::map<double,double> extractp;
	std::map<double,double> extractu;

	std::map<double,double>::iterator it;
	while (ifile >> coord1 >> coord2) {
		coord.first = coord1*chord;
		coord.second = coord2*chord;
		line.push_back(coord);
	}
	for (std::vector<std::pair<double,double>>::iterator it = line.begin(); it != line.end()-1; it++) {
		x1 = (*it).first;
		y1 = (*it).second;
		x2 = (*(it+1)).first;
		y2 = (*(it+1)).second;
		for (int i = 0; i < N_elem; i++) {
			for (int j = 0; j < 3; j++) {
				if (points[elements[i].vertex[j]].x >= x1 && points[elements[i].vertex[j]].x <= x2) {
					for (int k = 0; k < 3; k++) {
						xcheck = points[elements[i].vertex[k]].x;
						ycheck = points[elements[i].vertex[k]].y;
						value = (y2-y1)*xcheck - (x2-x1)*ycheck + (x2*y1 - x1*y2); //Equation of line
						if (k == 0) flag = (value > 0) - (value < 0);
						if (((value > 0) - (value < 0)) != flag) {
							extractp[elements[i].xc/chord] = (elements[i].p - p_inf)/(0.5*rho_inf*u_inf*u_inf);
							extractu[elements[i].xc/chord] = elements[i].u;
							break;
						}
					}
					break;
				}
			}
		}
	}
	for (it = extractp.begin(); it != extractp.end(); it++) {
		ofile << (*it).first << " " << (*it).second << std::endl;
	}
	for (it = extractu.begin(); it != extractu.end(); it++) {
		ofile2 << (*it).first << " " << (*it).second << std::endl;
	}
}

#endif