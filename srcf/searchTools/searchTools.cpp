//============================================================================
// Name        : searchTools.cpp
// Author      : Omar Haque
// Version     : 2.0
// Copyright   : ECL-2.0
// Description : The main tools and search algorithms.
//============================================================================

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
#include <mpi.h>
#include <stdio.h>
#include "../entropy/entropy.h"
#include "../utilities/utilities.h"



/**
    Returns an array of information theoretic measures

    @param p1 The index of the first site
    @param Hp1 The entropy of the first site
    @param p2 The index of the second site
    @param Hp2 The entropy of the second site
    @param p3 The index of the third site
    @param cl The index of the class
    @param p2 The index of the second site
    @param p3 The index of the third site
    @param *d Pointer to the first element of the data
    @param nsamples The number of samples (rows)
    @param nvars The number of variables (columns) including Class
    @param Hcl The entropy of the class column
    @param Hp1cl The joint entropy of the first site and Class column

    @return Array containing [IGStrict,IGAlt,IGNew,VIp,VIc]
*/
double *calculateMeasures(int p1, double Hp1, int p2, double Hp2, int p3, int cl, const int *d,
		int nsamples, int nvars, int c,double Hcl, double Hp1cl){
	static double final[5];
	int selTempNew[4] = {-1,-1,-1,-1};

	/*** calculation of base entropies and joint entropies ***/
	// H(P2,C)
	selTempNew[0] = p2;
	selTempNew[1] = cl;
	double Hp2cl = entropyFast(d,nsamples,nvars,0,selTempNew);

	// H(P3,C)
	selTempNew[0] = p3;
	double Hp3cl = entropyFast(d,nsamples,nvars,0,selTempNew);

	// H(P3)
	selTempNew[1] = -1;
	double Hp3 = entropyFast(d,nsamples, nvars,0,selTempNew);

	// H(P1,P2)
	selTempNew[0] = p1;
	selTempNew[1] = p2;
	double Hp1p2 = entropyFast(d,nsamples, nvars,0,selTempNew);

	//H(P1,P2,C)
	selTempNew[2] = cl;
	double Hp1p2cl = entropyFast(d,nsamples,nvars,0,selTempNew);

	// H(P1,P3)
	selTempNew[2] = -1;
	selTempNew[1] = p3;
	double Hp1p3 = entropyFast(d,nsamples,nvars,0,selTempNew);

	// H(P1,P3,C)
	selTempNew[2] = cl;
	double Hp1p3cl = entropyFast(d,nsamples,nvars,0,selTempNew);

	// H(P2,P3)
	selTempNew[0] = p2;
	selTempNew[2] = -1;
	double Hp2p3 = entropyFast(d,nsamples,nvars,0,selTempNew);

	// H(P2,P3,cl)
	selTempNew[2] = cl;
	double Hp2p3cl = entropyFast(d,nsamples,nvars,0,selTempNew);

	// H(P1,P2,P3)
	selTempNew[2] = p1;
	double Hp1p2p3 = entropyFast(d,nsamples,nvars,0,selTempNew);

	// H(P1,P2,P3,C)
	selTempNew[3] = cl;
	double Hp1p2p3cl = entropyFast(d,nsamples,nvars,0,selTempNew);

	///////////// calculation of compound measures
	double Ip1p2p3 = Hcl - Hp1p2p3cl + Hp1p2p3;
	double Ip1 = Hcl - Hp1cl + Hp1;
	double Ip2 = Hcl - Hp2cl + Hp2;
	double Ip3 = Hcl - Hp3cl + Hp3;
	double Ip1p2 = Hcl - Hp1p2cl + Hp1p2;
	double Ip1p3 = Hcl - Hp1p3cl + Hp1p3;
	double Ip2p3 = Hcl - Hp2p3cl + Hp2p3;
	double IGp1p2 = Ip1p2 - Ip1 - Ip2;
	double IGp1p3 = Ip1p3 - Ip1 - Ip3;
	double IGp2p3 = Ip2p3 - Ip2 - Ip3;

	///////////// calculation of VI
	double VIp1p2 = 2*Hp1p2 - Hp1 - Hp2;
	double VIp1p3 = 2*Hp1p3 - Hp1 - Hp3;
	double VIp2p3 = 2*Hp2p3 - Hp2 - Hp3;
	double VIp1cl = 2*Hp1cl - Hp1 - Hcl;
	double VIp2cl = 2*Hp2cl - Hp2 - Hcl;
	double VIp3cl = 2*Hp3cl - Hp3 - Hcl;

	// calculation of the measures
	double IG_strict = Ip1p2p3 - max(IGp1p2,0.0) - max(IGp1p3,0.0) - max(IGp2p3,0.0) - Ip1 - Ip2 - Ip3;
	double IG_alt = Ip1p2p3 - IGp1p2 - IGp1p3 - IGp2p3 - Ip1 - Ip2 - Ip3;
	double IG_easy = Ip1p2p3 - IGp1p2 - IGp1p3 - IGp2p3;
	double VIp = VIp1p2 + VIp1p3 + VIp2p3;
	double VIc = VIp1cl + VIp2cl + VIp3cl;

	final[0] = IG_strict;
	final[1] = IG_alt;
	final[2] = IG_easy;
	final[3] = VIp;
	final[4] = VIc;
	return final;
}

/**
 * Run through all of the triples of the dataset, calculating measures for each triple.
 *
 * @param filename The filename for the genetic dataset
 * @param outputFilename The filename for the printed measures to be outputted
 * @param nsamples Number of samples in the genetic dataset (number of rows)
 * @param nvars Number of variables in the genetic dataset (number of columns) including Class
 * @param c The estimator to use, 0 uses empirical entropy. See entropy.cpp for details
 * @param printall Indicates whether to print measures for every triple
 * @param assocLevel If printall is false, measures are printed if IGStrict > assocLevel*H(CL)
 */
void runFullSearch(std::string filename, std::string outputFilename, int nsamples, int nvars,
		int c, bool printall, double assocLevel, int numThread = 1){

	// initialise variables
	std::vector<int> d = readData(filename, nsamples, nvars);
	const int* p = d.data();
	int v[4] = {-1,-1,-1,-1};
	v[0] = nvars-1;
	double Hcl = entropyFast(p,nsamples,nvars,c,v);
	v[0] = -1;
	double measure0sum = 0.0;
	double measure0squaredsum = 0.0;
	double measure1sum = 0.0;
	double measure1squaredsum = 0.0;
	double measure2sum = 0.0;
	double measure2squaredsum = 0.0;
	double measure3sum = 0.0;
	double measure3squaredsum = 0.0;
	double measure4sum = 0.0;
	double measure4squaredsum = 0.0;
	int iterations = 0;


	ofstream myfile (outputFilename.c_str());

	if (myfile.is_open()){

	#pragma omp parallel for reduction(+:measure0sum,measure0squaredsum,measure1sum,measure1squaredsum,measure2sum,measure2squaredsum,measure3sum,measure3squaredsum,measure4sum,measure4squaredsum) private(v) num_threads(numThread)
	for (int i = 0; i < nvars-3;i++){
		v[0] = i;
		double Hp1 = entropyFast(p,nsamples,nvars,c,v);
		v[1] = nvars-1;
		double Hp1cl = entropyFast(p,nsamples,nvars,c,v);

		v[0] = -1;
		v[1] = -1;

		for (int j = i+1;j < nvars-2;j++){

			v[0] = j;
			double Hp2 = entropyFast(p,nsamples,nvars,c,v);
			v[0] = -1;

			for(int k = j+1; k < nvars-1; k++){
				double *measureArray = calculateMeasures(i, Hp1, j, Hp2, k, nvars-1, p,
						nsamples, nvars, c, Hcl, Hp1cl);

				if (printall) {
					myfile << "(" << i << "," << j << "," << k << "): "
							<< *(measureArray + 0) << "	" << *(measureArray + 1) << "	"
							<< *(measureArray + 2) << "	" << *(measureArray + 3) << "	"
							<< *(measureArray + 4) << "\n";
				} else if (*(measureArray + 0) > assocLevel * Hcl) {
					myfile << "(" << i << "," << j << "," << k << "): "
							<< *(measureArray + 0) << "	" << *(measureArray + 1) << "	"
							<< *(measureArray + 2) << "	" << *(measureArray + 3) << "	"
							<< *(measureArray + 4) << "\n";
				}


				measure0sum+=*(measureArray+0);
				measure0squaredsum+=pow(*(measureArray+0),2);
				measure1sum+=*(measureArray+1);
				measure1squaredsum+=pow(*(measureArray+1),2);
				measure2sum+=*(measureArray+2);
				measure2squaredsum+=pow(*(measureArray+2),2);
				measure3sum+=*(measureArray+3);
				measure3squaredsum+=pow(*(measureArray+3),2);
				measure4sum+=*(measureArray+4);
				measure4squaredsum+=pow(*(measureArray+4),2);
				iterations++;
			}
		}
	}

	myfile.close();
	} else {
		cout << "Unable to open file";
	}
	writeStatistics(measure0sum, measure0squaredsum, measure1sum,
			measure1squaredsum, measure2sum, measure2squaredsum, measure3sum,
			measure3squaredsum, measure4sum, measure4squaredsum, iterations);
}

/**
 * Run through specified indexes in the dataset, calculating measures for each triple.
 *
 * @param filename The filename for the genetic dataset
 * @param outputFilename The filename for the printed measures to be outputted
 * @param nsamples Number of samples in the genetic dataset (number of rows)
 * @param nvars Number of variables in the genetic dataset (number of columns) including Class
 * @param c The estimator to use, 0 uses empirical entropy. See entropy.cpp for details
 * @param istart The index i to start at, i.e. (istart,jstart,kstart) - inclusive
 * @param iend The index i to end at, i.e. (iend,jend,kend) - inclusive
 * @param jstart The index j to start at, i.e. (istart,jstart,kstart) - inclusive
 * @param jend The index j to end at, i.e. (iend,jend,kend) - inclusive
 * @param kstart The index k to start at, i.e. (istart,jstart,kstart) - inclusive
 * @param kend The index k to end at, i.e. (iend,jend,kend) - inclusive
 * @param printall Indicates whether to print measures for every triple
 * @param assocLevel If printall is false, measures are printed if IGStrict > assocLevel*H(CL)
 */
void runFullSearchIndexes(std::string filename, std::string outputFilename, int nsamples, int nvars,
		int c, int istart, int iend, int jstart, int jend, int kstart, int kend,
		bool printall, double assocLevel){


	/* (istart,jstart,jend) to (iend,jend,kend) INCLUSIVE BOTH SIDES */


	double measure0sum = 0.0;
	double measure0squaredsum = 0.0;
	double measure1sum = 0.0;
	double measure1squaredsum = 0.0;
	double measure2sum = 0.0;
	double measure2squaredsum = 0.0;
	double measure3sum = 0.0;
	double measure3squaredsum = 0.0;
	double measure4sum = 0.0;
	double measure4squaredsum = 0.0;
	int iterations = 0;


	std::vector<int> d = readData(filename, nsamples, nvars);
	// read in the array of indexes through filenameIndexes

	const int* p = d.data();

	int v[4] = {-1,-1,-1,-1};
	v[0] = nvars-1;

	double Hcl = entropyFast(p,nsamples,nvars,c,v);
	v[0] = -1;

	ofstream myfile (outputFilename.c_str());

	bool Done = false;

	if (myfile.is_open()){

	int jstartreal;
	int kstartreal;

	for (int i = istart; !Done; i++){
		v[0] = i;
		double Hp1 = entropyFast(p,nsamples,nvars,c,v);
		v[1] = nvars-1;
		double Hp1cl = entropyFast(p,nsamples,nvars,c,v);

		v[0] = -1;
		v[1] = -1;

		if (i == istart){
			jstartreal = jstart;
		} else {
			jstartreal = i+1;
		}

		for (int j = jstartreal; (j < nvars-2) && !Done; j++){

			v[0] = j;
			double Hp2 = entropyFast(p,nsamples,nvars,c,v);
			v[0] = -1;

			if (j == jstart){
				kstartreal = kstart;
			} else {
				kstartreal = j+1;
			}

			for(int k = kstartreal; (k < nvars-1) && !Done; k++){
				double *measureArray = calculateMeasures(i, Hp1, j, Hp2, k, nvars-1, p,
						nsamples, nvars, c, Hcl, Hp1cl);

				if (printall) {
					myfile << "(" << i << "," << j << "," << k << "): "
							<< *(measureArray + 0) << "	" << *(measureArray + 1) << "	"
							<< *(measureArray + 2) << "	" << *(measureArray + 3) << "	"
							<< *(measureArray + 4) << "\n";
				} else if (*(measureArray + 0) > assocLevel * Hcl) {
					myfile << "(" << i << "," << j << "," << k << "): "
							<< *(measureArray + 0) << "	" << *(measureArray + 1) << "	"
							<< *(measureArray + 2) << "	" << *(measureArray + 3) << "	"
							<< *(measureArray + 4) << "\n";
				}


				// TODO: add overflow protection
				measure0sum+=*(measureArray+0);
				measure0squaredsum+=pow(*(measureArray+0),2);
				measure1sum+=*(measureArray+1);
				measure1squaredsum+=pow(*(measureArray+1),2);
				measure2sum+=*(measureArray+2);
				measure2squaredsum+=pow(*(measureArray+2),2);
				measure3sum+=*(measureArray+3);
				measure3squaredsum+=pow(*(measureArray+3),2);
				measure4sum+=*(measureArray+4);
				measure4squaredsum+=pow(*(measureArray+4),2);
				iterations++;

				Done = (i > iend) && (j > jend) && (k > kend);
			}
		}
	}

	myfile.close();
	} else {
		cout << "Unable to open file";
	}
	writeStatistics(measure0sum, measure0squaredsum, measure1sum,
			measure1squaredsum, measure2sum, measure2squaredsum, measure3sum,
			measure3squaredsum, measure4sum, measure4squaredsum, iterations);
}

/**
 * Compute epistatic measures using the VI filter. Only compute measures for triples such that
 * VI[i] + VI[j] + VI[k] > epsilon
 *
 * @param filename The input file
 * @param outputFilename The output file name
 * @param nsamples number of samples (rows) in the data
 * @param nvars number of variables (columns) INCLUDING phenotype. nvars = #snps+1
 * @param c estimator used. 0 is empirical entropy
 * @param printall denote whether to print all the snps you come across
 * @param assocLevel if printall is false, only print if I{*} > H(cl) * assocLevel
 * @param viDists vector of distances to the phenotype, assumed to be sorted increasing
 * @param numThreads to parallelise over i
 */
void runSearchVIFilter(std::string filename, std::string outputFilename, int nsamples, int nvars,
		int c,double epsilon,bool printall, double assocLevel, std::vector<double> &viDists,
		int numThread = 1){

// initialise variables
	double measure0sum = 0.0;
	double measure0squaredsum = 0.0;
	double measure1sum = 0.0;
	double measure1squaredsum = 0.0;
	double measure2sum = 0.0;
	double measure2squaredsum = 0.0;
	double measure3sum = 0.0;
	double measure3squaredsum = 0.0;
	double measure4sum = 0.0;
	double measure4squaredsum = 0.0;
	int iterations = 0;

	std::vector<double>::iterator bin;



	std::vector<int> d = readData(filename, nsamples, nvars);
	// read in the array of indexes through filenameIndexes

	const int* p = d.data();
	// sorted dataframe
	int v[4] = {-1,-1,-1,-1};
	v[0] = nvars-1;

	double Hcl = entropyFast(p,nsamples,nvars,c,v);
	v[0] = -1;

	ofstream myfile (outputFilename.c_str());



	if (myfile.is_open()){
	int kstartreal;

	#pragma omp parallel for reduction(+:measure0sum,measure0squaredsum,measure1sum,measure1squaredsum,measure2sum,measure2squaredsum,measure3sum,measure3squaredsum,measure4sum,measure4squaredsum) private(v,bin,kstartreal) num_threads(numThread)
	for (int i = 0; i < nvars-3; i++){
		v[0] = i;
		double Hp1 = entropyFast(p,nsamples,nvars,c,v);
		v[1] = nvars-1;
		double Hp1cl = entropyFast(p,nsamples,nvars,c,v);

		v[0] = -1;
		v[1] = -1;

		for (int j = i+1; j < nvars-2; j++){

			v[0] = j;
			double Hp2 = entropyFast(p,nsamples,nvars,c,v);
			v[0] = -1;


			bin = std::upper_bound(viDists.begin() + j+1, viDists.end(), epsilon - viDists[i] - viDists[j]);
			kstartreal = bin - viDists.begin();


		for(int k = kstartreal; k < nvars-1; k++){

						double *measureArray = calculateMeasures(i, Hp1, j, Hp2, k, nvars-1, p,
								nsamples, nvars, c, Hcl, Hp1cl);

						if (printall) {
							myfile << "(" << i << "," << j << "," << k << "): "
									<< *(measureArray + 0) << "	" << *(measureArray + 1) << "	"
									<< *(measureArray + 2) << "	" << *(measureArray + 3) << "	"
									<< *(measureArray + 4) << "\n";
						} else if (*(measureArray + 0) > assocLevel * Hcl) {
							myfile << "(" << i << "," << j << "," << k << "): "
									<< *(measureArray + 0) << "	" << *(measureArray + 1) << "	"
									<< *(measureArray + 2) << "	" << *(measureArray + 3) << "	"
									<< *(measureArray + 4) << "\n";
						}


						// TODO: add overflow protection
						measure0sum+=*(measureArray+0);
						measure0squaredsum+=pow(*(measureArray+0),2);
						measure1sum+=*(measureArray+1);
						measure1squaredsum+=pow(*(measureArray+1),2);
						measure2sum+=*(measureArray+2);
						measure2squaredsum+=pow(*(measureArray+2),2);
						measure3sum+=*(measureArray+3);
						measure3squaredsum+=pow(*(measureArray+3),2);
						measure4sum+=*(measureArray+4);
						measure4squaredsum+=pow(*(measureArray+4),2);
						iterations++;



				}
			}
		}
		myfile.close();
		}
			else {
				cout << "Unable to open file";
			}
		writeStatistics(measure0sum, measure0squaredsum, measure1sum,
				measure1squaredsum, measure2sum, measure2squaredsum, measure3sum,
				measure3squaredsum, measure4sum, measure4squaredsum, iterations);
	}

/**
 * Compute epistatic measures using the VI filter. Only compute measures for triples such that
 * VI[i] + VI[j] + VI[k] > epsilon. Only compute over specified range of indexes
 *
 * @param filename The input file
 * @param outputFilename The output file name
 * @param nsamples number of samples (rows) in the data
 * @param nvars number of variables (columns) INCLUDING phenotype. nvars = #snps+1
 * @param c estimator used. 0 is empirical entropy
 * @param istart The index i to start at, i.e. (istart,jstart,kstart) - inclusive
 * @param iend The index i to end at, i.e. (iend,jend,kend) - inclusive
 * @param jstart The index j to start at, i.e. (istart,jstart,kstart) - inclusive
 * @param jend The index j to end at, i.e. (iend,jend,kend) - inclusive
 * @param kstart The index k to start at, i.e. (istart,jstart,kstart) - inclusive
 * @param kend The index k to end at, i.e. (iend,jend,kend) - inclusive
 * @param printall denote whether to print all the snps you come across
 * @param assocLevel if printall is false, only print if I{*} > H(cl) * assocLevel
 * @param viDists vector of distances to the phenotype, assumed to be sorted increasing
 * @param numThreads to parallelise over i
 */
void runSearchVIFilterIndexes(std::string filename, std::string outputFilename, int nsamples, int nvars,
		int c, int istart, int iend, int jstart, int jend, int kstart,int kend, double epsilon,
		bool printall, double assocLevel, std::vector<double> &viDists){

	// initialise variables
	double measure0sum = 0.0;
	double measure0squaredsum = 0.0;
	double measure1sum = 0.0;
	double measure1squaredsum = 0.0;
	double measure2sum = 0.0;
	double measure2squaredsum = 0.0;
	double measure3sum = 0.0;
	double measure3squaredsum = 0.0;
	double measure4sum = 0.0;
	double measure4squaredsum = 0.0;
	int iterations = 0;

	std::vector<double>::iterator bin;



	std::vector<int> d = readData(filename, nsamples, nvars);
	// read in the array of indexes through filenameIndexes

	const int* p = d.data();
	// sorted dataframe
	int v[4] = {-1,-1,-1,-1};
	v[0] = nvars-1;

	double Hcl = entropyFast(p,nsamples,nvars,c,v);
	v[0] = -1;

	ofstream myfile (outputFilename.c_str());

	bool Done = false;

	if (myfile.is_open()){

	int jstartreal;
	int kstartcandidate;
	int kstartreal;

	for (int i = istart; (i <= iend) && !Done; i++){
		v[0] = i;
		double Hp1 = entropyFast(p,nsamples,nvars,c,v);
		v[1] = nvars-1;
		double Hp1cl = entropyFast(p,nsamples,nvars,c,v);

		v[0] = -1;
		v[1] = -1;

		if (i == istart){
			jstartreal = jstart;
		} else {
			jstartreal = i+1;
		}



		for (int j = jstartreal; (j < nvars-2) && !Done; j++){

			v[0] = j;
			double Hp2 = entropyFast(p,nsamples,nvars,c,v);
			v[0] = -1;


			if (j == jstart){
				kstartcandidate = kstart;
			} else {
				kstartcandidate = j+1;
			}

			viDists.begin();
			bin = std::upper_bound(viDists.begin() + kstartcandidate, viDists.end(), epsilon - viDists[i] - viDists[j]);
			kstartreal = bin - viDists.begin();




		for(int k = kstartreal; (k < nvars-1) && !Done; k++){


						double *measureArray = calculateMeasures(i, Hp1, j, Hp2, k, nvars-1, p,
								nsamples, nvars, c, Hcl, Hp1cl);

						if (printall) {
							myfile << "(" << i << "," << j << "," << k << "): "
									<< *(measureArray + 0) << "	" << *(measureArray + 1) << "	"
									<< *(measureArray + 2) << "	" << *(measureArray + 3) << "	"
									<< *(measureArray + 4) << "\n";
						} else if (*(measureArray + 0) > assocLevel * Hcl) {
							myfile << "(" << i << "," << j << "," << k << "): "
									<< *(measureArray + 0) << "	" << *(measureArray + 1) << "	"
									<< *(measureArray + 2) << "	" << *(measureArray + 3) << "	"
									<< *(measureArray + 4) << "\n";
						}


						// TODO: add overflow protection
						measure0sum+=*(measureArray+0);
						measure0squaredsum+=pow(*(measureArray+0),2);
						measure1sum+=*(measureArray+1);
						measure1squaredsum+=pow(*(measureArray+1),2);
						measure2sum+=*(measureArray+2);
						measure2squaredsum+=pow(*(measureArray+2),2);
						measure3sum+=*(measureArray+3);
						measure3squaredsum+=pow(*(measureArray+3),2);
						measure4sum+=*(measureArray+4);
						measure4squaredsum+=pow(*(measureArray+4),2);
						iterations++;

						Done = (i > iend) && (j > jend) && (k > kend);

				}
			}
		}
	myfile.close();
	}
		else {
			cout << "Unable to open file";
		}
	writeStatistics(measure0sum, measure0squaredsum, measure1sum,
			measure1squaredsum, measure2sum, measure2squaredsum, measure3sum,
			measure3squaredsum, measure4sum, measure4squaredsum, iterations);
}



/**
 * Compute epistatic measures using the Cluster filter. Only compute measures for triples such that
 * clusterID[i]!=clusterID[j]!=clusterID[k]
 *
 * @param filename The input file
 * @param outputFilename The output file name
 * @param nsamples Number of samples (rows) in the data
 * @param nvars Number of variables (columns) INCLUDING phenotype. nvars = #snps+1
 * @param c Estimator used. 0 is empirical entropy
 * @param printall Denote whether to print all the snps you come across
 * @param assocLevel If printall is false, only print if I{*} > H(cl) * assocLevel
 * @param clusterIDs Vector containing the cluster IDs for all the SNPs
 * @param clusterSizes Vector of length numClusters containing the size of all clusters
 * @param numThread The number of threads to parallelise over. Defaults to 1.
 */
void runSearchClusterFilter(std::string filename, std::string outputFilename, int nsamples, int nvars,
		int c, bool printall, double assocLevel, const std::vector<int> &clusterIDs,
		const std::vector<int> &clusterSizes, int numThread){

	// initialise variables
	double measure0sum = 0.0;
	double measure0squaredsum = 0.0;
	double measure1sum = 0.0;
	double measure1squaredsum = 0.0;
	double measure2sum = 0.0;
	double measure2squaredsum = 0.0;
	double measure3sum = 0.0;
	double measure3squaredsum = 0.0;
	double measure4sum = 0.0;
	double measure4squaredsum = 0.0;
	int iterations = 0;



	std::vector<int> d = readData(filename, nsamples, nvars);
	// read in the array of indexes through filenameIndexes

	const int* p = d.data();

	int v[4] = {-1,-1,-1,-1};
	v[0] = nvars-1;

	double Hcl = entropyFast(p,nsamples,nvars,c,v);
	v[0] = -1;

	ofstream myfile (outputFilename.c_str());

	if (myfile.is_open()){

	int jstartreal;
	int kstartreal;


	std::vector<int> nextClusterIndices(clusterSizes.size());
	int sum = 0;
	for (int i = 0; i < nextClusterIndices.size(); i++){
		sum+=clusterSizes[i];
		nextClusterIndices[i] = sum;

	}


	#pragma omp parallel for reduction(+:measure0sum,measure0squaredsum,measure1sum,measure1squaredsum,measure2sum,measure2squaredsum,measure3sum,measure3squaredsum,measure4sum,measure4squaredsum) private(v,jstartreal,kstartreal) num_threads(numThread)
	for (int i = 0; i < nvars-3; i++){
		v[0] = i;
		double Hp1 = entropyFast(p,nsamples,nvars,c,v);
		v[1] = nvars-1;
		double Hp1cl = entropyFast(p,nsamples,nvars,c,v);

		v[0] = -1;
		v[1] = -1;


		jstartreal = nextClusterIndices[clusterIDs[i]];

		for (int j = jstartreal; j < nvars-2; j++){

					v[0] = j;
					double Hp2 = entropyFast(p,nsamples,nvars,c,v);
					v[0] = -1;



			kstartreal = nextClusterIndices[clusterIDs[j]];




		for(int k = kstartreal; k < nvars-1; k++){




						double *measureArray = calculateMeasures(i, Hp1, j, Hp2, k, nvars-1, p,
								nsamples, nvars, c, Hcl, Hp1cl);

						if (printall) {
							myfile << "(" << i << "," << j << "," << k << "): "
									<< *(measureArray + 0) << "	" << *(measureArray + 1) << "	"
									<< *(measureArray + 2) << "	" << *(measureArray + 3) << "	"
									<< *(measureArray + 4) << "\n";
						} else if (*(measureArray + 0) > assocLevel * Hcl) {
							myfile << "(" << i << "," << j << "," << k << "): "
									<< *(measureArray + 0) << "	" << *(measureArray + 1) << "	"
									<< *(measureArray + 2) << "	" << *(measureArray + 3) << "	"
									<< *(measureArray + 4) << "\n";
						}


						// TODO: add overflow protection
						measure0sum+=*(measureArray+0);
						measure0squaredsum+=pow(*(measureArray+0),2);
						measure1sum+=*(measureArray+1);
						measure1squaredsum+=pow(*(measureArray+1),2);
						measure2sum+=*(measureArray+2);
						measure2squaredsum+=pow(*(measureArray+2),2);
						measure3sum+=*(measureArray+3);
						measure3squaredsum+=pow(*(measureArray+3),2);
						measure4sum+=*(measureArray+4);
						measure4squaredsum+=pow(*(measureArray+4),2);
						iterations++;

		}
	}
	}
		myfile.close();
} else {
		cout << "Unable to open file";
	}
	writeStatistics(measure0sum, measure0squaredsum, measure1sum,
			measure1squaredsum, measure2sum, measure2squaredsum, measure3sum,
			measure3squaredsum, measure4sum, measure4squaredsum, iterations);

	}



/**
 * Compute epistatic measures using the Cluster filter. Only compute measures for triples such that
 * clusterID[i]!=clusterID[j]!=clusterID[k]
 *
 * @param filename The input file
 * @param outputFilename The output file name
 * @param nsamples Number of samples (rows) in the data
 * @param nvars Number of variables (columns) INCLUDING phenotype. nvars = #snps+1
 * @param c Estimator used. 0 is empirical entropy
 * @param istart The index i to start at, i.e. (istart,jstart,kstart) - inclusive
 * @param iend The index i to end at, i.e. (iend,jend,kend) - inclusive
 * @param jstart The index j to start at, i.e. (istart,jstart,kstart) - inclusive
 * @param jend The index j to end at, i.e. (iend,jend,kend) - inclusive
 * @param kstart The index k to start at, i.e. (istart,jstart,kstart) - inclusive
 * @param kend The index k to end at, i.e. (iend,jend,kend) - inclusive
 *
 * @param printall Denote whether to print all the snps you come across
 * @param assocLevel If printall is false, only print if I{*} > H(cl) * assocLevel
 * @param clusterIDs Vector containing the cluster IDs for all the SNPs
 * @param clusterSizes Vector of length numClusters containing the size of all clusters
 * @param numThread The number of threads to parallelise over. Defaults to 1.
 */
void runSearchClusterFilterIndexes(std::string filename, std::string outputFilename, int nsamples, int nvars,
		int c, int istart, int iend, int jstart, int jend, int kstart,int kend,
		bool printall, double assocLevel, const std::vector<int> &clusterIDs, const std::vector<int> &clusterSizes){

	// initialise variables
	double measure0sum = 0.0;
	double measure0squaredsum = 0.0;
	double measure1sum = 0.0;
	double measure1squaredsum = 0.0;
	double measure2sum = 0.0;
	double measure2squaredsum = 0.0;
	double measure3sum = 0.0;
	double measure3squaredsum = 0.0;
	double measure4sum = 0.0;
	double measure4squaredsum = 0.0;
	int iterations = 0;

	int numClusters = (clusterSizes).size();
	int iUpperLim = getIndexNextCluster(numClusters - 3,clusterIDs,clusterSizes,numClusters,nvars); // upper limit i exclusive
	int jUpperLim = getIndexNextCluster(numClusters - 2,clusterIDs,clusterSizes,numClusters,nvars); // upper limit j exclusive
	int iTrueUpperLim = min(iUpperLim,iend);
	int jTrueUpperLim = min(jUpperLim,jend);
	int kTrueUpperLim = min(nvars-2,kend);


	std::vector<int> nextClusterIndices(clusterSizes.size());
	int sum = 0;
	for (int i = 0; i < nextClusterIndices.size(); i++){
		sum+=clusterSizes[i];
		nextClusterIndices[i] = sum;

	}

	std::vector<int> d = readData(filename, nsamples, nvars);
	// read in the array of indexes through filenameIndexes

	const int* p = d.data();

	int v[4] = {-1,-1,-1,-1};
	v[0] = nvars-1;

	double Hcl = entropyFast(p,nsamples,nvars,c,v);
	v[0] = -1;

	ofstream myfile (outputFilename.c_str());

	bool Done = false;

	if (myfile.is_open()){

	int jstartreal;
	int kstartreal;

	for (int i = istart; (i <= iend) && !Done; i++){
		v[0] = i;
		double Hp1 = entropyFast(p,nsamples,nvars,c,v);
		v[1] = nvars-1;
		double Hp1cl = entropyFast(p,nsamples,nvars,c,v);

		v[0] = -1;
		v[1] = -1;

		if (i == istart){

		jstartreal = max(nextClusterIndices[clusterIDs[i]],jstart);


		} else {

			jstartreal = nextClusterIndices[clusterIDs[i]];


		}


		for (int j = jstartreal; (j < jend) && !Done; j++){

					v[0] = j;
					double Hp2 = entropyFast(p,nsamples,nvars,c,v);
					v[0] = -1;

					if (j == jstart){

						kstartreal = max(nextClusterIndices[j],kstart);

					} else {
						kstartreal = nextClusterIndices[j];


					}

		for(int k = kstartreal; (k < nvars-1) && !Done; k++){




						double *measureArray = calculateMeasures(i, Hp1, j, Hp2, k, nvars-1, p,
								nsamples, nvars, c, Hcl, Hp1cl);

						if (printall) {
							myfile << "(" << i << "," << j << "," << k << "): "
									<< *(measureArray + 0) << "	" << *(measureArray + 1) << "	"
									<< *(measureArray + 2) << "	" << *(measureArray + 3) << "	"
									<< *(measureArray + 4) << "\n";
						} else if (*(measureArray + 0) > assocLevel * Hcl) {
							myfile << "(" << i << "," << j << "," << k << "): "
									<< *(measureArray + 0) << "	" << *(measureArray + 1) << "	"
									<< *(measureArray + 2) << "	" << *(measureArray + 3) << "	"
									<< *(measureArray + 4) << "\n";
						}


						// TODO: add overflow protection
						measure0sum+=*(measureArray+0);
						measure0squaredsum+=pow(*(measureArray+0),2);
						measure1sum+=*(measureArray+1);
						measure1squaredsum+=pow(*(measureArray+1),2);
						measure2sum+=*(measureArray+2);
						measure2squaredsum+=pow(*(measureArray+2),2);
						measure3sum+=*(measureArray+3);
						measure3squaredsum+=pow(*(measureArray+3),2);
						measure4sum+=*(measureArray+4);
						measure4squaredsum+=pow(*(measureArray+4),2);
						iterations++;



						Done = (i >= iTrueUpperLim) && (j >= jTrueUpperLim) && (k >= kTrueUpperLim);

		}
	}
	}
		myfile.close();
} else {
		cout << "Unable to open file";
	}
	writeStatistics(measure0sum, measure0squaredsum, measure1sum,
			measure1squaredsum, measure2sum, measure2squaredsum, measure3sum,
			measure3squaredsum, measure4sum, measure4squaredsum, iterations);

	}



