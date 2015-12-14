#include<iostream>
#include<vector>
#include<math.h>
using namespace std;


int main() {
//	vector<int> vec;
//	for(int i = 0; i < 10; i++) {
//		vec.push_back(i*24);
//	}
//	std::vector<int>::iterator it = vec.end()-1;
//	cout <<vec[9] <<" "<<vec[it-vec.begin()]<<endl;

//	double M = -2.1;
//	cout << (M > 0) - (M < 0);

    double a,b,c,d;
    a = 8;
    b = 27.89991007;;
    c = 66.7899;
    d = 88.73;
    cout<<std::max(std::max(a,b),std::max(c,d))<<endl;
    cout <<"cbrt "<<cbrt(a);
}
