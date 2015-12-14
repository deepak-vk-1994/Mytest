#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<ctype.h>
#include<cstring>
#include <sstream>
#include <stdexcept>

#define STL
//#define NASTRAN
using namespace std;
double scale = 1;
int main() {
	#ifdef STL
		string line;
		ofstream out("inputvijay.txt");
		int flag = 1;
		int check = 1;
		ifstream file("inputvijay.grd");
		while(file) {
			std::getline(file,line);
			istringstream iss(line); //split the string and store in string streams
			string word;
			flag = 1;
			while(iss >> word) {	
				check = 1;
				if (word == "facet" || word == "outer" || word == "endloop" || word == "endfacet" || word == "solid") {
					flag = 0;
					break;
				}
				try {
					stod(word);   //string to double
				}
				catch(const std::invalid_argument& e) {
					check = 0;
				}
				if (check != 0) {
					double d = stod(word);
					out << d*scale << " ";
				}
			}
			if (flag != 0) out << endl;
		}
		file.close();
		out.close();
	#endif
//	#ifdef NASTRAN
//		string line;
//		ofstream out("outNastran.txt");
//		int flag = 1;
//		int check = 1;
//		ifstream file("bumpgrid_unstruct_nastran.grd");
//		while(file) {
//			std::getline(file,line);
//			istringstream iss(line); //split the string and store in string streams
//			string word;
//			flag = 1;
//			while(iss >> word) {	
//				check = 1;
//				if (word == "MAT1" || word == "PSHELL"){
//					flag = 0;
//					break;
//				}
//				
////				try {
////					stod(word);   //string to double
////				}
////				catch(const std::invalid_argument& e) {
////					check = 0;
////				}
////				if (check != 0) {
////					double d = stod(word);
////					out << d << " ";
////				}
//				out << word <<"  Y  " ;
//			}
//			if (flag != 0) out << endl;
//		}
//		file.close();
//		out.close();
//	#endif
}

