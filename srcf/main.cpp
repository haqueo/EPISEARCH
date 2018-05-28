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
#include <time.h>

using namespace std;


int main(int argc, char* argv[]) {

	// get command line arguments, warn user if something is wrong.
	if (argc < 6){
		std::cerr << "Usage: " << argv[0] << "-inputfile -outputfile -nsamples "
				"-nvars -nindex -indexfile -printall -assoclevel" <<
				std::endl;
		return 1;
	}
	std::string inputfile = "";
	std::string outputfile = "";
	int nsamples = -1;
	int nvars = -1;
	int nindex = -1000;
	std::string indexfile = "";
	std::string printallDec = "";
	bool printall = false;
	double assoclevel = 0.05;

	cout << "argc is " << argc << std::endl;

	for (int i = 1; i < argc; i++){
		if (std::string(argv[i]) == "-inputfile"){
			if (i+1 < argc){
				inputfile = argv[i+1];
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
		}else if (std::string(argv[i]) == "-nindex"){
			if (i+1 < argc){
				std::istringstream iss(argv[i+1]);
				iss >> nindex;
				i++;


			} else {
				std::cerr << "-nindex option requires one argument" << std::endl;
				return 1;
			}
		}else if (std::string(argv[i]) == "-indexfile"){
			if (i+1 < argc){

				indexfile = argv[i+1];
				i++;

			} else {
				std::cerr << "-indexfile option requires one argument" << std::endl;
				return 1;
			}
		}else if (std::string(argv[i]) == "-printall"){
			if (i+1 < argc){

				printallDec = argv[i+1];
				i++;

				if (printallDec == "y"){
					printall=true;
				} else if (printallDec == "n"){
					printall= false;
				} else{
					std::cerr << "input to -printall is y or n" << std::endl;
					return 1;
				}


			} else {
				std::cerr << "-printall option requires one argument" << std::endl;
				return 1;
			}
		} else if (std::string(argv[i]) == "-assoclevel"){
			if (i+1 < argc){
				std::istringstream iss(argv[i+1]);
				iss >> assoclevel;
				i++;

			} else {
				std::cerr << "-assoclevel option requires one argument" << std::endl;
				return 1;
			}}  else{
			// error, none of the above.

			std::cerr << "input not in valid format";
			return 1;
		}
	}

	clock_t t1,t2;
	t1=clock();
	    //code goes here


	if (nindex == -1000){
		// no index pair given. Running full search.
		runFullSearch(inputfile, outputfile, nsamples, nvars, 0,printall,assoclevel);
	} else {
		// use nindex to find the 6 indices
		std::vector<int> limits = readIndices(indexfile,nindex);
		cout << "limits vector is" << std::endl;
		for(int i=0; i<6; ++i)
		  std::cout << limits[i] << ' ';


		// then run full search indexes
		runFullSearchIndexes(inputfile,outputfile,nsamples,nvars,0,limits[0],
				limits[3],limits[1],limits[4],limits[2],limits[5],printall,assoclevel);
	}
	t2=clock();
	float diff ((float)t2-(float)t1);
	float seconds = diff / CLOCKS_PER_SEC;

	ofstream myfile ("statisticsAn.txt");

	if (myfile.is_open()){
		myfile << "time taken: "<< seconds << " s";
		myfile.close();
	}

	return 0;

}
