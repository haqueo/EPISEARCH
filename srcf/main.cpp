#include <iostream>

#include "entropy/fullSearch.h"
using namespace std;

int main() {
	cout << "!!!Hello World!!!" << endl; // prints !!!Hello World!!!
	// double entropy(const int *d, int nsamples, int nvars, int c, bool *v) {
	const int d[12] = {0,1,1,1,1,0,0,0,0,1,0,1};
	int nsamples = 6;
	int nvars = 2;
	int c= 0;
	bool v[2] = {true,true};

	double H = entropy(d,nsamples,nvars,c,v);

	cout << "Entropy is " << H;
	return 0;

}
