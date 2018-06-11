#include <iostream>
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
#include <stdio.h>
#include <stdlib.h>

#include <stdio.h>
#include "searchTools/searchTools.h"

#include "utilities/utilities.h"

using namespace std;


int main(int argc, char* argv[]) {

	// get command line arguments, warn user if something is wrong.
	if (argc < 6){
		std::cerr << "Usage: " << argv[0] << "-inputfile -outputfile -nsamples "
				"-nvars -nindex -indexfile -printall -assoclevel -VIfilter -VIfilterfile -epsilon"
				" -clusterfilter -clusterfilterIDsFile -clusterfilterSizeFile -numClusters -numThread" <<
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
	std::string VIfilterDec = "";
	bool VIfilter = false;
	std::string VIfilterfile = "";
	double epsilon = 0.0;
	std::string clusterFilterDec = "";
	bool clusterFilter = false;
	std::string clusterFilterIDsFile = "";
	std::string clusterFilterSizeFile = "";
	int numClusters = -1;
	int numThread = 1;

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
			}}  else if (std::string(argv[i]) == "-VIfilter"){
			if (i+1 < argc){

						VIfilterDec = argv[i+1];
						i++;

						if (VIfilterDec == "y"){
							VIfilter=true;
						} else if (VIfilterDec == "n"){
							VIfilter= false;
						} else{
							std::cerr << "input to -VIfilter is y or n" << std::endl;
							return 1;
						}


					} else {
						std::cerr << "-VIfilter option requires one argument" << std::endl;
						return 1;
					}
			} else if (std::string(argv[i]) == "-VIfilterfile"){

				if (i+1 < argc){
							VIfilterfile = argv[i+1];
							i++;
						} else {
							std::cerr << "-VIfilterfile option requires one argument" << std::endl;
							return 1;
						}
			} else if (std::string(argv[i]) == "-epsilon"){
				if (i+1 < argc){
					std::istringstream iss(argv[i+1]);
					iss >> epsilon;
					i++;
				} else {
					std::cerr << "-epsilon option requires one argument" << std::endl;
					return 1;
				}
			} else if (std::string(argv[i]) == "-clusterfilter"){
				if (i+1 < argc){

					clusterFilterDec = argv[i+1];
					i++;
					if(clusterFilterDec == "y"){
						clusterFilter=true;
					} else if (clusterFilterDec == "n"){
						clusterFilter = false;
					} else{
						std::cerr << "input to -clusterfilter is y or n" << std::endl;
						return 1;
					}
				} else {
					std::cerr << "input to -clusterfilter option requires one argument" << std::endl;
					return 1;
				}
			}else if (std::string(argv[i]) == "-clusterfilterIDsFile"){

				if(i+1 < argc){
					clusterFilterIDsFile = argv[i+1];
					i++;
				} else {
					std::cerr << "-clusterfilterIDsFile option requires one argument";
				}
			}else if (std::string(argv[i]) == "-clusterfilterSizeFile"){

				if (i+1 < argc){
					clusterFilterSizeFile = argv[i+1];
					i++;
				} else {
					std::cerr << "-clusterfilterSizeFile option requires one argument";
				}
			}else if (std::string(argv[i]) == "-numClusters"){
				if (i+1 < argc){
					std::istringstream iss(argv[i+1]);
					iss >> numClusters;
					i++;
				} else {
					std::cerr << "-numClusters option requires one argument" << std::endl;
					return 1;
				}
			} else if (std::string(argv[i]) == "-numThread"){
				if (i+1 < argc){
					std::istringstream iss(argv[i+1]);
					iss >> numThread;
					i++;

				} else {
					std::cerr << "-numThread option requires one argument" << std::endl;
					return 1;
				}
			} else {
			// error, none of the above.
				std::cerr << std::string(argv[i]) << std::endl;
			std::cerr << "input not in valid format";
			return 1;
		}
	}

	clock_t t1,t2;
	t1=clock();

	std::vector<int> limits;
	if (nindex != -1000){
		limits = readIndices(indexfile,nindex);
	}

	if (nindex != -1000){
		// use indexed functions
		if (VIfilter){
			std::vector<double> viDists = readviDistances(VIfilterfile,nvars);
			runSearchVIFilterIndexes(inputfile, outputfile, nsamples, nvars, 0, limits[0],limits[3],limits[1],limits[4],limits[2],limits[5],epsilon,printall, assoclevel, viDists);

		} else if (clusterFilter){
			std::vector<int> clusterIDs = readClusterIDs(clusterFilterIDsFile,nvars);
			std::vector<int> clusterSizes = readClusterSizes(clusterFilterSizeFile,numClusters);
			runSearchClusterFilterIndexes(inputfile,outputfile,nsamples,nvars,0,limits[0],limits[3],limits[1],limits[4],
									limits[2],limits[5],printall,assoclevel,clusterIDs,clusterSizes);
		} else {
			runFullSearchIndexes(inputfile,outputfile,nsamples,nvars,0,limits[0],
							limits[3],limits[1],limits[4],limits[2],limits[5],printall,assoclevel);
		}
	} else {
		// use parallelised functions

		if (VIfilter){
			std::vector<double> viDists = readviDistances(VIfilterfile,nvars);
			runSearchVIFilter(inputfile,outputfile,nsamples,nvars,0,epsilon,printall,assoclevel,viDists,numThread);
		} else if (clusterFilter){
			std::vector<int> clusterIDs = readClusterIDs(clusterFilterIDsFile,nvars);

			std::vector<int> clusterSizes = readClusterSizes(clusterFilterSizeFile,numClusters);
			runSearchClusterFilter(inputfile,outputfile,nsamples,nvars,0,printall,assoclevel,clusterIDs,clusterSizes,numThread);
		} else {
			runFullSearch(inputfile,outputfile,nsamples,nvars,0,printall,assoclevel,numThread);
		}

	}

	t2=clock();
	float diff ((float)t2-(float)t1);
	float seconds = diff / CLOCKS_PER_SEC;

	ofstream myfile ("statistics.txt");

	if (myfile.is_open()){
		myfile << "time taken: "<< seconds << " s";
		myfile.close();
	}


	return 0;

}
