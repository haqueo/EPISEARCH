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
#include <sstream>

using namespace std;


int main(int argc, char* argv[]) {

	// get command line arguments, warn user if something is wrong.
	if (argc < 5){
		std::cerr << "Usage: " << argv[0] << "-inputfile -outputfile -nsamples -nvars" <<
				std::endl;
		return 1;
	}
	std::string inputfile = "";
	std::string outputfile = "";
	int nsamples = -1;
	int nvars = -1;

	cout << "argc is " << argc << std::endl;

	for (int i = 1; i < argc; i++){
		if (std::string(argv[i]) == "-inputfile"){
			if (i+1 < argc){
				inputfile = argv[i+1];
				cout << "this worked, inputfile" << std::endl;
				i++;
			} else {
				std::cerr << "-input option requires one argument" << std::endl;
				return 1;
			}
		} else if (std::string(argv[i]) == "-outputfile"){
			if (i+1 < argc){
				outputfile = argv[i+1];
				cout << "this worked, outputfile" << std::endl;
				i++;
			} else {
				std::cerr << "-outputfile option requires one argument" << std::endl;
				return 1;
			}
		} else if (std::string(argv[i]) == "-nsamples"){
			if (i+1 < argc){
				std::istringstream iss(argv[i+1]);
				iss >> nsamples;
				i++;
				cout << "this worked, nsamples" << std::endl;


			} else {
				std::cerr << "-nsamples option requires one argument" << std::endl;
				return 1;
			}
		}else if (std::string(argv[i]) == "-nvars"){
			if (i+1 < argc){
				std::istringstream iss(argv[i+1]);
				iss >> nvars;
				i++;
				cout << "this worked, nvars" << std::endl;

			} else {
				std::cerr << "-input option requires one argument" << std::endl;
				return 1;
			}
		} else{
			// error, none of the above.

			std::cerr << "input not in valid format";
			return 1;
		}
	}

	cout << "The values I got for the command line were" << std::endl;
	cout << "inputfile = " << inputfile << std::endl;
	cout << "outputfile = " << outputfile << std::endl;
	cout << "nsamples = " << nsamples << std::endl;
	cout << "nvars = " << nvars << std::endl;
	// run full search

	// end

	runFullSearch(inputfile, outputfile, nsamples, nvars, 0);
	//std::string outputtestersubset = "/Users/Omar/Documents/Year4/M4R/fullSearch/output/outputtestersubset.txt";
	//runFullSearchIndexes(inputfile,outputtestersubset,nsamples,nvars,0,12,61,30,75,41,97);

	return 0;

}
