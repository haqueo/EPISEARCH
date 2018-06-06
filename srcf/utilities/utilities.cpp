/*
 * utilities.cpp
 *
 *  Created on: 6 Jun 2018
 *      Author: Omar
 */
#include <iostream>
//using namespace std;
#include <algorithm>
#include <ctime>
#include <cmath>
#include <vector>
#include <map>
#include <fstream>
using namespace std;
#include <string>
#include <sstream>
#include <iterator>
#include <math.h>
#include "../entropy/entropy.h"

/**
 * Reads in a genetic data file of the form nrows x nvar
 *
 * @param filename The File to be read in
 * @param nrows The number of rows in the .txt
 * @param nvars The number of columns in the.txt
 *
 * @return Pointer to the data in a 1D format, column first.
 */
std::vector<int> readData(std::string filename, int nrows, int nvars){

	const int datasize = nvars*nrows;
	static std::vector<int> data(datasize);
	static std::vector<int> dataReal(datasize);

    std::fstream myfile(filename.c_str(), std::ios_base::in);
    int a;

    // read in the data as is.
    for(int i =0; i< datasize;i++){
    	myfile >> a;
    	data[i] = a;
    }

    // need to transform the data so it's column first.
    for(int s = 0; s < nrows; s++){
    	for (int i = 0; i<nvars;i++){
    		dataReal[s+i*nrows] = data[i+s*nvars];
    	}
    }
	return dataReal;

}



/**
 * Read in the n'th index pair from indexes.txt
 *
 * @param filename The location of indexes.txt
 * @param n How many lines down to get indices from
 *
 * @return The n'th line from indexes.txt
 */
std::vector<int> readIndices(std::string filename, int n){

	std::ifstream infile(filename.c_str());
	std::string s;
	std::vector<int> indexes;

   //skip N lines
   for(int i = 0; i < n; ++i){
	   std::getline(infile, s);
   }

   std::replace( s.begin(), s.end(), ',', ' ' );
   stringstream ss(s);
   copy(istream_iterator<int>(ss), istream_iterator<int>(), back_inserter(indexes));

   return indexes;
}


std::vector<double> readviDistances(std::string filename){
	std::ifstream ifile(filename.c_str(), std::ios::in);
	std::vector<double> viDists;

    double num = 0.0;
    //keep storing values from the text file so long as data exists:
    while (ifile >> num) {
        viDists.push_back(num);
    }

    return viDists;

}


std::vector<int> readClusterIDs(std::string filename){

	std::ifstream ifile(filename.c_str(),std::ios::in);
	std::vector<int> clusterIDs;

	int num = 0;
	while(ifile >> num){
		clusterIDs.push_back(num);
	}
	cout << "finished reading clusterIDs" <<std::endl;
	return clusterIDs;
}

std::vector<int> readClusterSizes(std::string filename){
	std::ifstream ifile(filename.c_str(),std::ios::in);
	std::vector<int> clusterSizes;

	int num = 0;
	while(ifile >> num){
		clusterSizes.push_back(num);
	}
	cout << "finished reading clusterSizes" << std::endl;


	return clusterSizes;
}


// Get the index of the first element in the next cluster. If we're currently in the last
// cluster, return a number greater than the amount of SNPs. The sanity check
// in the "Done" variable will prevent further iterations.
int getIndexNextCluster(int thisClusterID, std::vector<int> clusterIDs, std::vector<int> clusterSizes,
		int numClusters, int nvars){


	if(thisClusterID == numClusters-1){
		return nvars+1;

	}
	int sumCounter = 0;

	for(int i = 0; i <= thisClusterID; i++){
		sumCounter += clusterSizes[i];
	}

	return sumCounter;

}

