#include <iostream>
#include "entropy/fullSearch.h"
#include <string>
#include <fstream>
//using namespace std;
#include <algorithm>
#include <ctime>
#include <cmath>
#include <vector>
#include <map>

using namespace std;


int main() {


	int nsamples = 469;
	int nvars = 101;
	int c = 0;
	std::vector<int> d = readData("/Users/Omar/Documents/Year4/M4R/fullSearch/data/prjebintTEST.txt",
	469, 101);
//
	const int* p = d.data();
//
	//ofstream myfile ("/Users/Omar/Documents/Year4/M4R/fullSearch/output/prjebintTESToutput/output.txt");
	bool v[101] = { false };
	v[100] = true;
	double Hcl = entropy(p,nsamples,nvars,c,v);
	v[100] = false;

	for (int i = 0; i < 100-2;i++){
		v[i] = true;
		double Hp1 = entropy(p,nsamples,nvars,c,v); // this is 0??
		v[100] = true;
		double Hp1cl = entropy(p,nsamples,nvars,c,v);

		v[100] = false;
		v[i] = false;

		for (int j = i+1;j < 100-1;j++){

			v[j] = true;
			double Hp2 = entropy(p,nsamples,nvars,c,v);
			v[j] = false;
			for(int k = j+1; k < 100; k++){
				double *measureArray = calculateMeasures(i, Hp1, j, Hp2, k, 100, p,
						nsamples, nvars, c, Hcl, Hp1cl);

				cout << "(" << i << "," << j << ","<< k << "): " << *(measureArray+0) <<
						"	"<< *(measureArray+1) <<"	"<< *(measureArray+3) <<
						"	"<<*(measureArray+4) << "\n";
			}
		}
	}

//	myfile.close();
//

	return 0;

}
